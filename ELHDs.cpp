//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 30.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>
#include <vector>

// Inclusion from Inside Loop library
#include <curses.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/norm.h>
#include <il/math.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "Dilatancy.h"
#include "ELHDs.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
#include "Mesh.h"

namespace hfp2d {

/* This function solves the 2D problem of fluid injection into a frictional
 * weakening dilatant fault for one time increment.
 *
 * Inputs -> look at "TimeIncr" function
 *
 * Output:
 *
 *   Structure called 'Result', which is composed of:
 *
 *     - friction -> vector of friction coefficient at each time step
 *     - dilatancy -> vector of dilatancy at each time step
 *     - tot_stress_state -> array (matrix) of total stress at each time step
 *     - incr_d ->  vectr of increment of slip
 *     - P ->  vector of pore pressure profile at nodal points for each time
 *             step
 *     - iter -> number of iterations for each time step to satisfy MC
 *     - final_slip_points -> vector of final slip points per each time step
 *     - t -> current time
 *     - dt -> current time step
 *     - slippagezone -> length slippage zone for each time step
 *
 */

void elhds(Mesh mesh, int p, double Cohes, const il::Array2D<double> &kmat,
           double Incr_dil, double d_wd, il::Array2D<double> rho,
           double Init_dil, double CompressFluid, double TimeStep, double Visc,
           il::Array<double> S, int InjPoint, int dof_dim, double Peak_fric,
           double Resid_fric, double d_wf, il::io_t, Result &res) {

  // Vector of friction at the beginning of each time step
  il::Array<double> fric{res.friction.size(), 0};
  fric = res.friction;

  // Vector of dilatancy at the beginning of each time step
  il::Array<double> dil{res.dilatancy.size()};
  dil = res.dilatancy;

  // Matrix of total stress at the beginning of each time step
  // Size -> Nelts x 2
  il::Array2D<double> sigma_tot{res.tot_stress_state.size(0),
                                res.tot_stress_state.size(1), 0};
  sigma_tot = res.tot_stress_state;

  // Vector of total slip at nodal points at the beginning of each time step
  il::Array<double> incr_slip{res.incr_d.size(), 0};
  incr_slip = res.incr_d;

  // Length of slippage zone at the beginning of each time step
  double sl;
  sl = res.slippagezone;

  // Vector of pressure profile at nodes at the beginning of each time step
  // Size -> Nelts + 1 -> {P_n1, P_n2, P_n3 ..}
  il::Array<double> press_prof{res.P.size(), 0};
  press_prof = res.P;

  // Create the matrix of pressure at collocation points {{0,p_c1},{0,p_c2}..}
  // Size -> NCollPoints x 2
  il::Array2D<double> Pcm{2 * mesh.nelts(), 2, 0};
  for (il::int_t l = 0, k = 1; l < Pcm.size(0); ++l) {

    Pcm(l, 1) = il::dot(
        hfp2d::from_edge_to_col_cg(
            dof_dim,
            hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io),
            hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io), il::io),
        press_prof)[k];
    k = k + 2;
  }

  // Create matrix of effective stress for each collocation point {tau,sigma_n}
  // Size -> NCollPoints x 2
  il::Array2D<double> sigma_eff{2 * mesh.nelts(), 2, 0};
  for (int m = 0; m < sigma_eff.size(0); ++m) {
    for (int i = 0; i < sigma_eff.size(1); ++i) {

      sigma_eff(m, i) = sigma_tot(m, i) - Pcm(m, i);
    }
  }

  // Initial effective stress state at the beginning of each time step
  il::Array2D<double> sigma_eff_old{2 * mesh.nelts(), 2, 0};
  sigma_eff_old = sigma_eff;

  // Vector of shear stress at the beginning of each time step & vector of
  // effective normal stress at the beginning of each time step
  il::Array<double> tau{2 * mesh.nelts(), 0};
  il::Array<double> sigma_n{2 * mesh.nelts(), 0};
  for (il::int_t n = 0; n < 2 * mesh.nelts(); ++n) {

    tau[n] = sigma_eff(n, 0);
    sigma_n[n] = sigma_eff(n, 1);
  }

  // Auxiliary vector to flag slipping elements
  il::Array<il::int_t> kk{0};

  // Slipping collocation points
  il::Array<il::int_t> q{0};

  // Auxiliary vector to flag DOFs of slipping elements
  il::Array<int> Dof_slip_elem{0};

  // Declaration of variables
  double SL;
  double DTnew;
  il::Array<il::int_t> slipp_elem;
  il::Status status;
  il::Norm norm;
  norm = il::Norm::L1;

  // Initialization of the while loop
  int iter = 0;
  double itermax = 10e10;

  // Check is the criterion that has to be satisfied in the while loop:
  // tau must be lower (or equal) than (fric * sigma_n + Cohes) at each
  // collocation point
  int check;
  check = boole_mc(tau, sigma_n, fric, Cohes, il::io);

  il::Array<double> Dd{incr_slip.size(), 0};
  il::Array<double> DiffDd{Dd.size(), 0};
  il::Array<double> Ddold = Dd;
  il::Array<double> Dp{press_prof.size(), 0};
  il::Array<double> Dpold{press_prof.size(), 0};
  il::Array<double> DiffDp{Dp.size(), 0};
  il::Array<double> pnodes = press_prof;

  // Iterate until each point falls below the Mohr-Coulomb failure line
  while (check != 2 * mesh.nelts() && iter <= itermax) {

    // Calculate the excess of shear stress that has to be release
    il::Array<double> T{2 * mesh.nelts(), 0};
    for (il::int_t j = 0, k = 0, l = 0; j < T.size(); ++j) {

      T[j] = (Cohes + fric[j] * sigma_eff(j, 1)) - sigma_eff(j, 0);

      // If the jth component of T is negative, then flag slipping elements (the
      // ones in which the corresponding collocation point do not satisfy the
      // M-C criterion)
      if (T[j] <= 0) {

        kk.resize(k + 1);
        q.resize(k + 1);
        q[k] = j;
        kk[k] = hfp2d::position_2d_array(
            hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io), j, il::io)(0,
                                                                           0);
        k = k + 1;
      }
    }

    // Calculate the excess of shear stress that has to be release
    il::Array<double> TT{2 * mesh.nelts(), 0};
    for (il::int_t j = 0, k = 0, l = 0; j < T.size(); ++j) {

      TT[j] = (Cohes + fric[j] * sigma_eff_old(j, 1)) - sigma_eff_old(j, 0);
    }

    // Since we have 2 collocation points per element,
    // we have to delete the duplicate elements
    slipp_elem = delete_duplicates2(kk, il::io);

    // Find the corresponding DOFs of the slipping elements
    for (il::int_t i = 0, k = 0; i < slipp_elem.size(); ++i) {

      Dof_slip_elem.resize(k + 4);

      Dof_slip_elem[k] = hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p,
                                                    il::io)(slipp_elem[i], 0);
      Dof_slip_elem[k + 1] = hfp2d::dofhandle_dg_full2d(
          dof_dim, mesh.nelts(), p, il::io)(slipp_elem[i], 1);
      Dof_slip_elem[k + 2] = hfp2d::dofhandle_dg_full2d(
          dof_dim, mesh.nelts(), p, il::io)(slipp_elem[i], 2);
      Dof_slip_elem[k + 3] = hfp2d::dofhandle_dg_full2d(
          dof_dim, mesh.nelts(), p, il::io)(slipp_elem[i], 3);
      k = k + 4;
    }

    //// SOLUTION OF THE COUPLED PROBLEM ////
    // Initialization of the system BigA*BigX = BigB
    // Just for shear DDs!!
    il::Array2D<double> BigA{(Dof_slip_elem.size() / 2) + mesh.nelts() + 1,
                             (Dof_slip_elem.size() / 2) + mesh.nelts() + 1, 0};
    il::Array<double> BigB{(Dof_slip_elem.size() / 2) + mesh.nelts() + 1, 0};

    // Select the elasticity matrix for just shear DDs of slipping DOFs
    il::Array2D<double> kmatd{(Dof_slip_elem.size() / 2),
                              (Dof_slip_elem.size() / 2), 0};
    for (int m = 0, k = 0; m < kmatd.size(0); ++m, k = k + 2) {
      for (int i = 0, q = 0; i < kmatd.size(1); ++i, q = q + 2) {

        kmatd(m, i) = kmat(Dof_slip_elem[k], Dof_slip_elem[q]);
      }
    }

    // Assembling just elasticity part of BigA
    hfp2d::set_submatrix_non_linear_system(BigA, 0, 0, kmatd);

    // Implicit solution of the non-linear system of equations by using fixed
    // point iterations

    // Initialization of Finite Volume matrices
    il::Array2D<double> Vd;
    il::Array2D<double> Vp;
    il::Array2D<double> L{mesh.nelts() + 1, mesh.nelts() + 1, 0};
    il::Array2D<double> incrdk{mesh.nelts(), 2, 0};

    double betarela = 0.6;
    double errDd = 2;
    double errDp = 2;
    double tolerance = 10e-6;
    int j = 0;

    il::Array<double> BigX{BigA.size(0), 0};
    il::Array<double> Ddd{Dd.size(), 0};

    while (j < itermax && (errDd > tolerance || errDp > tolerance)) {
      ++j;

      for (il::int_t i = 0; i < incrdk.size(0); ++i) {
        for (il::int_t l = 0; l < incrdk.size(1); ++l) {

          incrdk(i, l) =
              incr_slip[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io)(i,
                                                                           l)] +
              Dd[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io)(i, l)];
        }
      }

      // Mass matrix "Vd"
      Vd = hfp2d::build_vd_matrix_p1(
          mesh, Incr_dil, d_wd,
          hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io), rho,
          incrdk, il::io);

      // Pressure matrix "P"
      Vp = hfp2d::build_vp_matrix_p1(mesh, Incr_dil, Init_dil, CompressFluid,
                                     incrdk, d_wd, il::io);
      // Finite Difference matrix "L"
      L = hfp2d::build_l_matrix(mesh, incrdk, rho, Visc, Incr_dil, d_wd,
                                Init_dil, TimeStep, il::io);

      /// Assembling the system ///

      for (il::int_t k = 0; k < Vd.size(0); ++k) {

        for (il::int_t i = 0, q = 0; i < (Dof_slip_elem.size() / 2);
             ++i, q = q + 2) {

          BigA((Dof_slip_elem.size() / 2) + k, i) = Vd(k, Dof_slip_elem[q]);
        }
      }

      for (il::int_t m = 0; m < Vp.size(0); ++m) {

        for (il::int_t i = 0; i < Vp.size(1); ++i) {

          BigA((Dof_slip_elem.size() / 2) + m, (Dof_slip_elem.size() / 2) + i) =
              Vp(m, i) - L(m, i);
        }
      }

      for (il::int_t n = 0, q = 0; n < slipp_elem.size(); ++n, q = q + 2) {

        BigB[q] =
            TT[hfp2d::dofhandle_dg(2, mesh.nelts(), il::io)(slipp_elem[n], 0)];
        BigB[q + 1] =
            TT[hfp2d::dofhandle_dg(2, mesh.nelts(), il::io)(slipp_elem[n], 1)];
      }

      for (il::int_t i1 = 0; i1 < L.size(0); ++i1) {

        BigB[(Dof_slip_elem.size() / 2) + i1] =
            il::dot(L, pnodes)[i1] + (TimeStep * S[i1]);
      }

      /// Boundary conditions ///
      for (il::int_t k1 = 0; k1 < BigA.size(1); ++k1) {
        BigA((Dof_slip_elem.size() / 2) + InjPoint, k1) = 0;
      }

      for (il::int_t l1 = 0; l1 < BigA.size(0); ++l1) {
        BigA(l1, (Dof_slip_elem.size() / 2) + InjPoint) = 0;
      }

      BigA((Dof_slip_elem.size() / 2) + InjPoint,
           (Dof_slip_elem.size() / 2) + InjPoint) = 1;

      BigB[(Dof_slip_elem.size() / 2) + InjPoint] = 0;

      /// Solve the system ///
      BigX = il::linear_solve(BigA, BigB, il::io, status);
      status.abort_on_error();

      // Under relaxation technique
      for (il::int_t m1 = 0, q = 0; m1 < slipp_elem.size(); ++m1, q = q + 2) {
        Ddd[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io)(slipp_elem[m1],
                                                               0)] = BigX[q];
        Ddd[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(),
                                il::io)(slipp_elem[m1], 1)] = BigX[q + 1];
      }

      for (il::int_t m1 = 0; m1 < Dd.size(); ++m1) {
        Dd[m1] = (1 - betarela) * Ddold[m1] + betarela * Ddd[m1];
      }

      for (il::int_t n1 = 0; n1 < Dp.size(); ++n1) {

        Dp[n1] = (1 - betarela) * Dpold[n1] +
                 betarela * BigX[(Dof_slip_elem.size() / 2) + n1];
      }

      // Error on increments
      for (il::int_t i2 = 0; i2 < DiffDd.size(); ++i2) {
        DiffDd[i2] = Dd[i2] - Ddold[i2];
      }

      for (il::int_t j2 = 0; j2 < DiffDp.size(); ++j2) {
        DiffDp[j2] = Dp[j2] - Dpold[j2];
      }

      errDd = (il::norm(DiffDd, norm) / il::norm(Dd, norm));
      errDp = (il::norm(DiffDp, norm)) / (il::norm(Dp, norm));

      // Update
      Ddold = Dd;
      Dpold = Dp;
    }

    // Calculate the increment of shear stress due to the increment of slip
    il::Array2D<double> IncrT{2 * mesh.nelts(), 2, 0};
    for (il::int_t k2 = 0, k = 0; k2 < IncrT.size(0); ++k2, k = k + 2) {
      for (il::int_t i = 0, q = 0; i < 2 * mesh.nelts(); ++i, q = q + 2) {
        IncrT(k2, 0) = IncrT(k2, 0) + kmat(k, q) * Dd[i];
      }
    }

    // Force the slipping nodes at the first iteration to stay in the MC
    // line in the next iterations
    //    il::Array2D<double> MC{q[q.size() - 1] - q[0] + 1, 2, 0};
    il::Array2D<double> MC{Dof_slip_elem.size() / 2, 2, 0};
    for (il::int_t m2 = 0, z2 = 0; m2 < MC.size(0); ++m2, z2 = z2 + 2) {
      MC(m2, 0) = fabs(
          Cohes + fric[Dof_slip_elem[z2] / 2] * sigma_n[Dof_slip_elem[z2] / 2]);
      MC(m2, 1) = fabs(sigma_n[Dof_slip_elem[z2] / 2]);
    }

    // Update the effective stress state
    for (il::int_t n2 = 0; n2 < sigma_eff.size(0); ++n2) {
      for (il::int_t i = 0; i < sigma_eff.size(1); ++i) {
        sigma_eff(n2, i) = sigma_eff_old(n2, i) + IncrT(n2, i);
      }
    }

    for (il::int_t k3 = 0; k3 < MC.size(0); ++k3) {
      for (il::int_t i = 0; i < sigma_eff.size(1); ++i) {
        sigma_eff(Dof_slip_elem[0] / 2 + k3, i) = MC(k3, i);
      }
    }

    // Update the vector of shear stress & vector of normal stress
    for (il::int_t i3 = 0; i3 < sigma_eff.size(0); ++i3) {
      tau[i3] = sigma_eff(i3, 0);
      sigma_n[i3] = sigma_eff(i3, 1);
    }

    // Update the pore pressure profile at nodal points
    for (il::int_t l2 = 0; l2 < press_prof.size(); ++l2) {
      press_prof[l2] = press_prof[l2] + Dp[l2];
    }

    // Total stress state at collocation points for next time increment
    for (il::int_t j3 = 0; j3 < sigma_tot.size(0); ++j3) {
      for (il::int_t i = 0; i < sigma_tot.size(1); ++i) {
        sigma_tot(j3, i) = sigma_eff(j3, i) + Pcm(j3, i);
      }
    }

    // Number of slipping elements
    int NSlippingElem;
    NSlippingElem = slipp_elem.size();

    // Number of collocation points that do not satisfy M-C criterion
    int NSlipCollPoints;
    NSlipCollPoints = q.size();

    // Calculate slippage length
    SL = euclidean_distance(mesh.Coor(slipp_elem[0], 0),
                            mesh.Coor(slipp_elem[slipp_elem.size() - 1] + 1, 0),
                            mesh.Coor(slipp_elem[0], 1),
                            mesh.Coor(slipp_elem[slipp_elem.size() - 1] + 1, 1),
                            il::io);

    // Adaptive time stepping
    double psi = 1.;
    if (SL - sl > psi) {
      DTnew = TimeStep / 10;
    } else {
      DTnew = TimeStep;
    }

    // Update the iterations
    ++iter;

    /*
            // Cumulative slip at nodal points
            for (il::int_t k3 = 0; k3 < Dd.size(); ++k3) {
              Dd[k3] = Dd[k3] + incr_slip[k3];
            }

            // Cumulative slip at collocation points
            il::Array<double> s{2 * mesh.nelts(), 0};
            s = il::dot(hfp2d::from_edge_to_col_dg(
                            dof_dim, hfp2d::dofhandle_dg(dof_dim, mesh.nelts(),
                            il::io),
                            il::io),
                        Dd);

            // Update the friction weakening coefficient for next time step
            fric = hfp2d::exp_friction(Peak_fric, Resid_fric, d_wf, s, il::io);

            // Update the dilatancy weakening coefficient for next time step
            dil = hfp2d::dilatancy(Init_dil, Incr_dil, d_wd, s, il::io);
    */

    // Update the check criterion
    check = boole_mc(tau, sigma_n, fric, Cohes, il::io);
  }

  // Cumulative slip at nodal points
  for (il::int_t k3 = 0; k3 < Dd.size(); ++k3) {
    Dd[k3] = Dd[k3] + incr_slip[k3];
  }

  // Cumulative slip at collocation points
  il::Array<double> s{2 * mesh.nelts(), 0};
  s = il::dot(
      hfp2d::from_edge_to_col_dg(
          dof_dim, hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io), il::io),
      Dd);

  // Update the friction weakening coefficient for next time step
  fric = hfp2d::exp_friction(Peak_fric, Resid_fric, d_wf, s, il::io);

  // Update the dilatancy weakening coefficient for next time step
  dil = hfp2d::dilatancy(Init_dil, Incr_dil, d_wd, s, il::io);

  // Assign the desired outputs to structure's members
  res.friction = fric;
  res.dilatancy = dil;
  res.dt = DTnew;
  res.iter = iter;
  res.incr_d = Dd;
  res.tot_stress_state = sigma_tot;
  res.slippagezone = SL;
  res.final_slip_elements = slipp_elem;
  res.P = press_prof;
}

///// OTHER UTILITIES /////

/// 1
// This routine turns a 2D array of double precision values into a vector of
// double precision values.
// Input:
//   - Arr -> 2D array that we want to "flat"
il::Array<double> flatten1(const il::Array2D<double> &Arr, il::io_t) {

  il::Array2D<int> A{Arr.size(0), 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * 2;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  il::Array<double> out{2 * Arr.size(0), 0};

  for (il::int_t i = 0; i < Arr.size(0); ++i) {

    for (il::int_t j = 0; j < 2; ++j) {

      out[A(i, j)] = Arr(i, j);
    }
  }

  return out;
}

/// 2
// This routine turns a 2D array of integer values into a vector of
// integer values.
// Input:
//   - Arr -> 2D array that we want to "flat"
il::Array<int> flatten2(const il::Array2D<int> &Arr, il::io_t) {

  il::Array2D<int> A{Arr.size(0), 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * 2;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  il::Array<int> out{2 * Arr.size(0), 0};

  for (il::int_t i = 0; i < Arr.size(0); ++i) {

    for (il::int_t j = 0; j < 2; ++j) {

      out[A(i, j)] = Arr(i, j);
    }
  }

  return out;
}

/// 3
// This function returns an integer which represents the number of collocation
// points that satisfy the Mohr-Coulomb criterion
// Inputs:
//  - tau -> vector of shear stress at each collocation point
//  - sigma_n -> vector of normal stress at each collocation point
//  - fric -> vector of friction coefficient at each collocation point
//  - Cohes -> cohesion (double precision value)
//  - NCollPoints -> number of collocation points
int boole_mc(il::Array<double> tau, il::Array<double> sigma_n,
             il::Array<double> fric, double Cohes, il::io_t) {

  il::Array<int> T{tau.size(), 0};
  int out = 0.;

  for (il::int_t i = 0; i < T.size(); ++i) {

    if (tau[i] <= (Cohes + (fric[i] * sigma_n[i]))) {

      T[i] = 1;
    }
  }

  for (il::int_t m = 0; m < T.size(); ++m) {

    if (T[m] == 1) {

      out = out + T[m];
    }
  }

  return out;
}

/// 4
// This function extracts the integer values greater (or equal) to 0 in a
// vector composed of integer values
// Input:
//  - arr -> vector in which we want to extract the values >=0
// It is general in a sense that the output can be either a vector or an
// integer
il::Array<int> select(const il::Array<int> &arr, il::io_t) {

  il::Array<int> ans{0};
  il::int_t k = 0;

  for (il::int_t i = 0; i < arr.size(); ++i) {

    if (arr[i] >= 0) {

      ans.resize(k + 1);
      ans[k] = arr[i];

      ++k;
    }
  }

  return ans;
};

/// 5
// This function allows us to set a submatrix of the A matrix as B matrix
// Inputs:
//  - A -> matrix in which we want to set the submatrix
//  - i0 -> initial row index in the A matrix where we want to set the submatrix
//  - i1 -> initial column index in the A matrix where we want to set
//          the submatrix
//  - B -> submatrix
void set_submatrix_non_linear_system(il::Array2D<double> &A, int i0, int i1,
                                     const il::Array2D<double> &B) {

  IL_EXPECT_FAST(i0 + B.size(0) <= A.size(0));
  IL_EXPECT_FAST(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
}

/// 6
// This function allows us to set a subvector of the A vector as B vector
// Inputs:
//  - A -> vector in which we want to set the subvector
//  - i0 -> initial index in the A vector where we want to set the subvector
//  - B -> subvector
void set_subvector_non_linear_system(il::Array<double> &A, int i0,
                                     const il::Array<double> &B) {

  IL_EXPECT_FAST(i0 + B.size() <= A.size());

  for (int j0 = 0; j0 < B.size(); ++j0) {
    A[i0 + j0] = B[j0];
  }
}

/// 7
// This function allows us to extract a submatrix from the A matrix
// Inputs:
//  - i0 -> initial row index where we want to extract the submatrix
//  - i1 -> final row index where we want to extract the submatrix
//  - j0 -> initial column index where we want to extract the submatrix
//  - j1 -> final column index where we want to extract the submatrix
//  - A -> matrix in which we want to extract the submatrix
il::Array2D<double>
take_submatrix_non_linear_system(int i0, int i1, int j0, int j1,
                                 const il::Array2D<double> &A) {

  il::Array2D<double> sub{i1 - i0 + 1, j1 - j0 + 1, 0};

  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      sub(i - i0, j - j0) = A(i, j);
    }
  }
  return sub;
}

/// 8
// This function allows us to extract a subvector from the A vector
// Inputs:
//  - i0 -> initial index where we want to extract the vector
//  - i1 -> final index where we want to extract the vector
//  - A -> vector in which we want to extract the subvector
il::Array<double> take_subvector_non_linear_system(int i0, int i1,
                                                   const il::Array<double> &A) {

  il::Array<double> sub{i1 - i0 + 1, 0};

  for (int i = i0; i <= i1; ++i) {

    sub[i - i0] = A[i];
  }
  return sub;
}

/// 9
// This funciton calcualtes the 2D euclidean distance between two points
// {x1,y1}, {x2,y2}
// Input -> x- y- coordinates of the two points
// Output -> double precision values that represents the euclidean distance
// between the two aforementioned points
double euclidean_distance(double x1, double x2, double y1, double y2,
                          il::io_t) {

  double dist;

  double x = x1 - x2;
  double y = y1 - y2;

  dist = pow(x, 2) + pow(y, 2);
  dist = sqrt(dist);

  return dist;
}

/// 10
// Function that return the index ( IndRow,IndColumn ) of a given value ("seek")
// in an array2D.
// arr2D -> array2D in which we want to find the index of  given value
// seek -> value for which we want to find out the index
// It return an array that contain the {N.row, N.col} of the seek value

il::Array<int> find_2d(const il::Array2D<double> &arr2D, double_t seek,
                       il::io_t) {

  il::Array<int> outp{2};

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (arr2D(i, j) == seek)
        outp[0] = i, outp[1] = j;
    }
  }

  return outp;
}

/// 11

il::Array<int> delete_duplicates(const il::Array<il::int_t> &arr, il::io_t) {

  il::Array<int> out{0};

  for (il::int_t i = 0, k = 0; i < arr.size(); ++i) {

    for (il::int_t j = i + 1; j < arr.size(); ++j) {

      if (arr[i] == arr[j]) {
        out.resize(k + 1);
        out[k] = arr[i];
        k = k + 1;
      }
    }
  }

  return out;
}

/// 12

il::Array<il::int_t> delete_duplicates2(const il::Array<il::int_t> &arr,
                                        il::io_t) {
  il::Array<il::int_t> out{};

  for (il::int_t i = 0; i < arr.size(); ++i) {
    bool already_there = false;
    for (il::int_t j = 0; j < out.size(); ++j) {
      if (arr[i] == out[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
      out.append(arr[i]);
    }
  }

  return out;
}
}
