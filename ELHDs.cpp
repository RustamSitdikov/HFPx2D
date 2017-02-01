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

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/blas/norm.h>
#include <il/math.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "Dilatancy.h"
#include "ELHDs.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
#include "Mesh.h"

namespace hfp2d {

/* This function solves the 2D problem of fluid injection into a frictional
 * weakening dilatant fault for one increment of time.
 *
 * Inputs -> look at "TimeIncr" function
 *
 * Output:
 *
 *   Structure called 'Result', which is composed of:
 *
 *     - friction -> vector of friction coefficient for each time step
 *     - dilatancy -> vector of dilatancy for each time step
 *     - tot_stress_state -> array (matrix) of total stress for each time step
 *     - incr_d ->  matrix of increment of slip
 *     - P ->  vector of pore pressure profile at nodal points for each time
 * step
 *     - iter -> number of iterations for each time step to satisfy MC
 *     - final_slip_points -> vector of final slip points per each time step
 *     - t -> current time
 *     - dt -> current time step
 *     - slippagezone -> length slippage zone for each time step
 *
 * */

Result elhds(Result SolutionAtTj, Mesh mesh, const int p, const double Cohes,
             const int itermax, il::Array2D<double> &kmat,
             const double tolerance, const double Incr_dil, const double d_wd,
             il::Array2D<int> &Dof, il::Array2D<double> rho,
             const double Init_dil, const double CompressFluid,
             const double TimeStep, const double Visc, il::Array<double> S,
             const int InjPoint, const int dof_dim, const double t,
             const double Peak_fric, const double Resid_fric,
             const double d_wf) {

  Result res;

  // Get the collocation points' information from mesh
  il::Array2D<double> xe{2, 2, 0};
  hfp2d::SegmentCharacteristic mysege;
  // CollPoints -> Matrix of coord. of collocation points
  // {{xcoll_1,ycoll_1},{xcoll_2,ycoll_2}..}
  il::Array2D<double> CollPoints{2 * mesh.nelts(), 2, 0.};

  for (il::int_t e = 0, k = 0; e < mesh.nelts(); ++e) {

    hfp2d::take_submatrix(
        xe, mesh.conn(e, 0), mesh.conn(e, 1), 0, 1,
        mesh.Coor); // take the coordinates of element e from the mesh object

    mysege = hfp2d::get_segment_DD_characteristic(
        xe, p); // get the segment characteristic.

    // assemble the matrix
    for (int j = 0; j < 2; ++j) {

      CollPoints(k, j) = mysege.CollocationPoints(0, j);
      CollPoints(k + 1, j) = mysege.CollocationPoints(1, j);
    }

    k = k + 2;
  }

  // Total numbers of collocation points (i.e 2*Nelts)
  il::int_t NCollPoints = CollPoints.size(0);

  // Vector of friction at the beginning of each time step
  il::Array<double> fric;
  fric = SolutionAtTj.friction;

  // Vector of dilatancy at the beginning of each time step
  il::Array<double> dil;
  dil = SolutionAtTj.dilatancy;

  // Matrix of total stress at the beginning of each time step
  // Size -> Nelts x 2
  il::Array2D<double> sigma_tot;
  sigma_tot = SolutionAtTj.tot_stress_state;

  // Vector of total slip at nodal points at the beginning of each time step
  il::Array2D<double> incr_slip;
  incr_slip = SolutionAtTj.incr_d;

  // Length of slippage zone at the beginning of each time step
  double sl;
  sl = SolutionAtTj.slippagezone;

  // Vector of pressure profile at nodes at the beginning of each time step
  // Size -> 2*Nelts
  // {p1_left,p1_right,p2_left,p2_right..}
  // Remember p1_right = p2_left
  il::Array<double> press_prof;
  press_prof = SolutionAtTj.P;

  // Vectot of pressure profile at nodes at the beginning of each time step
  // Size -> Nnodes
  il::Array<double> pnodes{mesh.nelts() + 1, 0.};
  for (int i = 0, k = 0; i < pnodes.size() - 1; ++i) {

    pnodes[i] = press_prof[k];
    k = k + 2;
  }
  pnodes[pnodes.size()] = press_prof[press_prof.size()];

  // Vector of pressure profile at collocation points at the beginning
  // of each time step. Size -> NCollPoints = 2*Nelts
  il::Array<double> Pc{NCollPoints, 0.};
  Pc = il::dot(from_edge_to_col(mesh.nelts(), 2), press_prof);

  // Create the matrix of pressure at collocation points {{0,pc1},{0,pc2}..}
  // Size -> NCollPoints x 2
  il::Array2D<double> Pcm{Pc.size(), 2, 0.};
  for (il::int_t l = 0; l < Pcm.size(0); ++l) {

    Pcm(l, 1) = Pc[l];
  }

  // Create matrix of effective stress for each collocation point {tau,sigma_n}
  // Size -> NCollPoints x 2
  il::Array2D<double> sigma_eff{NCollPoints, 2, 0.};
  for (int m = 0; m < sigma_eff.size(0); ++m) {
    for (int i = 0; i < sigma_eff.size(1); ++i) {

      sigma_eff(m, i) = sigma_tot(m, i) - Pcm(m, i);
    }
  }

  // Initial effective stress state at the beginning of each time step
  il::Array2D<double> I{NCollPoints, 2, 0.};
  I = sigma_eff;

  // vector of shear stress at the beginning of each time step
  il::Array<double> tau{NCollPoints, 0.};
  for (il::int_t n = 0; n < tau.size(); ++n) {

    tau[n] = sigma_eff(n, 0);
  }

  // vector of effective normal stress at the beginning of each time step
  il::Array<double> sigma_n{NCollPoints, 0.};
  for (il::int_t i1 = 0; i1 < sigma_n.size(); ++i1) {

    sigma_n[i1] = sigma_eff(i1, 1);
  }

  il::Array<int> kk{NCollPoints, -1};      // Auxiliary vector
  il::Array2D<int> qq{NCollPoints, 2, -1}; // Auxiliary matrix

  // Declaration of variables
  il::Array<double> dd{2 * mesh.nelts(), 0.};
  double SL;
  il::Array<int> k;
  double DTnew;

  // Initialization of the while loop
  int iter = 1;
  int check;
  // check is the criterion that has to be satisfied in the while loop:
  // tau must be lower (or equal) than (fric + sigma_n*p) at each
  // collocation point
  check = boole_mc(tau, sigma_n, fric, Cohes, NCollPoints);

  // Iterate until each point falls below the Mohr-Coulomb failure line
  while (check != NCollPoints && iter <= itermax) {

    // Loop over collocation points to check the M-C criterion and to
    // find out the slipping nodes and their degrees of freedom
    for (int i = 0; i < NCollPoints; ++i) {

      if (sigma_eff(i, 0) >= Cohes + fric[i] * sigma_eff(i, 1)) {

        kk[i] = i;
        qq(i, 0) = 2 * i - 1;
        qq(i, 1) = 2 * i;
      }
    }

    il::Array<int> q{2 * NCollPoints, 0};
    q = flatten2(qq);

    // Select the component >= 0, i.e the slipping nodes
    k = select(kk);

    // Select the component >= 0, i.e the DOFs of slipping nodes
    il::Array<int> dofd;
    dofd = select(q);

    // Calculate the excess of shear stress that has to be release
    il::Array<double> T{NCollPoints, 0.};
    for (il::int_t j = 0; j < T.size(); ++j) {

      T[j] = (Cohes + fric[j] * I(j, 1)) - I(j, 0);
    }

    //// SOLUTION OF THE COUPLED PROBLEM ////
    // Initialization of the system BigABigX = BigB
    il::Array2D<double> BigA{k.size() + mesh.nelts() + 1,
                             k.size() + mesh.nelts() + 1, 0.};
    il::Array<double> BigB{k.size() + mesh.nelts() + 1, 0.};

    // Select the elasticity matrix for just shear DDs
    il::Array2D<double> kmatd{2 * mesh.nelts(), 2 * mesh.nelts(), 0.};
    for (int m = 0; m < kmatd.size(0); ++m) {
      for (int i = 0, k = 0, q = 1; i < kmatd.size(1); ++i) {

        kmatd(m, i) = kmat(k, q);
        k = k + 2;
        q = q + 2;
      }
    }

    // Assembling just elasticity part of BigA
    hfp2d::set_submatrix_non_linear_system(
        BigA, 1, 1, hfp2d::take_submatrix_non_linear_system(
                        k[0], k[k.size()], k[0], k[k.size()], kmatd));

    // Implicit solution of the non-linear system of equations
    il::Array2D<double> Dd{mesh.nelts(), 2, 0.};
    il::Array2D<double> Ddold{mesh.nelts(), 2, 0.};

    il::Array<double> Dp{pnodes.size(), 0.};
    il::Array<double> Dpold{pnodes.size(), 0.};

    il::Array2D<double> Vd, Vp, L;
    il::Array2D<double> incrdk{incr_slip.size(0), incr_slip.size(1), 0.};

    double betarela = 0.6;
    double errDd = 2;
    double errDp = 2;
    int j = 0;

    while ((j < itermax && errDd > tolerance) || (errDp > tolerance)) {

      ++j;

      for (il::int_t i = 0; i < incr_slip.size(0); ++i) {

        for (il::int_t l = 0; l < incr_slip.size(1); ++l) {

          incrdk(i, l) = incr_slip(i, l) + Dd(i, l);
        }
      }

      // Mass matrix "Vd"
      Vd = hfp2d::build_vd_matrix_p1(mesh, Incr_dil, d_wd, Dof, rho, incrdk);
      // Pressure matrix "P"
      Vp = hfp2d::build_vp_matrix_p1(mesh, Incr_dil, Init_dil, CompressFluid,
                                     incrdk, d_wd);
      // Finite Difference matrix "L"
      L = hfp2d::build_l_matrix(mesh, incrdk, rho, Visc, Incr_dil, d_wd,
                                Init_dil, TimeStep);

      /// Assembling the system ///
      hfp2d::set_submatrix_non_linear_system(
          BigA, k.size() + 1, 1, hfp2d::take_submatrix_non_linear_system(
                                     1, Vd.size(0), k[0], k[k.size()], Vd));

      il::Array2D<double> C{Vp.size(0), Vp.size(1), 0.};
      for (il::int_t m = 0; m < C.size(0); ++m) {

        for (il::int_t i = 0; i < C.size(1); ++i) {

          C(m, i) = Vp(m, i) - L(m, i);
        }
      }

      hfp2d::set_submatrix_non_linear_system(BigA, k.size() + 1, k.size() + 1,
                                             C);

      hfp2d::set_subvector_non_linear_system(
          BigB, 1,
          hfp2d::take_subvector_non_linear_system(k[0], k[k.size()], T));

      il::Array<double> R{mesh.nelts() + 1, 0.};
      for (il::int_t i1 = 0; i1 < R.size(); ++i1) {

        R[i1] = il::dot(L, pnodes)[i1] + (S[i1] * TimeStep);
      }

      hfp2d::set_subvector_non_linear_system(BigB, k.size() + 1, R);

      /// Boundary conditions ///
      // Dirichlet boundary conditions for solution with constant overpressure
      for (il::int_t k1 = 0; k1 < BigA.size(1); ++k1) {

        BigA(k.size() + InjPoint, k1) = 0.;
      }

      for (il::int_t l1 = 0; l1 < k.size(); ++l1) {

        BigA(l1, k.size() + InjPoint) = 0.;
      }

      BigA(k.size() + InjPoint, k.size() + InjPoint) = 1;
      BigB[k.size() + InjPoint] = 0.;

      /// Solve the system ///
      il::Array<double> BigX{BigA.size(0), 0.};
      il::Status status;
      BigX = il::linear_solve(BigA, BigB, il::io, status);

      // Select the contribution in terms of increment of nodal slip and
      // pressure
      il::Array<double> Dd2;
      il::Array<double> Dp2;
      Dd2 = hfp2d::take_subvector_non_linear_system(1, k.size(), BigX);
      Dp2 = hfp2d::take_subvector_non_linear_system(k.size() + 1, BigX.size(),
                                                    BigX);

      // under relaxation technique
      for (il::int_t m1 = 0; m1 < Dp.size(); ++m1) {

        Dp[m1] = (1 - betarela) * Dpold[m1] + betarela * Dp2[m1];
      }

      for (il::int_t n1 = 0; n1 < k.size(); ++n1) {

        dd[k[n1]] = Dd2[n1];
      }

      il::Array2D<int> Aux{mesh.nelts(), 2, 0};

      for (il::int_t k = 0, j; k < Aux.size(0); ++k) {

        j = k * dof_dim;

        for (il::int_t i = 0; i < Aux.size(1); ++i) {

          Aux(k, i) = i + j;
        }
      }

      for (il::int_t i2 = 0; i2 < Dd.size(0); ++i2) {

        for (il::int_t i = 0; i < Dd.size(1); ++i) {

          Dd(i2, i) = (1 - betarela) * Ddold(i2, i) + betarela * dd[Aux(i2, i)];
        }
      }

      // Error on increments
      double errDd, errDp;

      il::Array2D<double> diffDdDdold{Dd.size(0), Dd.size(1), 0.};
      for (il::int_t j2 = 0; j2 < diffDdDdold.size(0); ++j2) {

        for (il::int_t i = 0; i < diffDdDdold.size(1); ++i) {

          diffDdDdold(j2, i) = Dd(j2, i) - Ddold(j2, i);
        }
      }

      il::Array<double> diffDpDpold{Dp.size(), 0.};
      for (il::int_t k2 = 0; k2 < diffDpDpold.size(); ++k2) {

        diffDpDpold[k2] = Dp[k2] - Dpold[k2];
      }

      il::Norm norm;
      norm = il::Norm::L2;
      errDd = (il::norm(flatten1(diffDdDdold), norm) /
               il::norm(flatten1(Dd), norm));
      errDp = (il::norm(diffDpDpold, norm)) / (il::norm(Dp, norm));

      // Update
      Ddold = Dd;
      Dpold = Dp;
    }

    // Calculate the increment of shear stress due to the increment of slip
    il::Array2D<double> IncrT{NCollPoints, 2, 0.};
    il::Array<double> IncrTau{2 * mesh.nelts(), 0.};
    IncrTau = il::dot(kmatd, dd);

    for (il::int_t l2 = 0; l2 < IncrT.size(0); ++l2) {

      IncrT(l2, 0) = IncrTau[l2];
    }

    // Force the slipping nodes at the first iteration to stay in the MC line in
    // the second iteration
    il::Array2D<double> MC{k[k.size()] - k[0] + 1, 2, 0.};

    for (il::int_t m2 = 0; m2 < MC.size(0); ++m2) {

      MC(m2, 0) = fabs(Cohes + fric[k[m2]] * sigma_n[k[m2]]);
      MC(m2, 1) = fabs(sigma_n[k[m2]]);
    }

    // Update the effective stress state
    for (il::int_t n2 = 0; n2 < sigma_eff.size(0); ++n2) {
      for (il::int_t i = 0; i < sigma_eff.size(1); ++i) {

        sigma_eff(n2, i) = I(n2, i) + IncrT(n2, i);
      }
    }

    set_submatrix_non_linear_system(sigma_eff, k[0], k[0], MC);

    // Update the vector of shear stress & vector of normal stress
    for (il::int_t i3 = 0; i3 < sigma_eff.size(0); ++i3) {

      tau[i3] = sigma_eff(i3, 0);
      sigma_n[i3] = sigma_eff(i3, 1);
    }

    // Update the pore pressure profile at nodal points
    il::Array2D<double> DP2{mesh.nelts(), 2, 0.};
    for (il::int_t j3 = 0; j3 < DP2.size(0); ++j3) {

      for (il::int_t i = 0; i < DP2.size(1); ++i) {

        DP2(j3, i) = Dp[mesh.conn(j3, i)];
      }
    }

    il::Array<double> DP{2 * mesh.nelts(), 0.};
    DP = flatten1(DP2);

    for (il::int_t k3 = 0; k3 < press_prof.size(); ++k3) {

      press_prof[k3] = press_prof[k3] + DP[k3];
    }

    // Total stress state at collocation points
    for (il::int_t l3 = 0; l3 < sigma_tot.size(0); ++l3) {

      for (il::int_t i = 0; i < sigma_tot.size(1); ++i) {

        sigma_tot(l3, i) = sigma_eff(l3, i) + Pcm(l3, i);
      }
    }

    // Number of slip points
    int NSlipPoints;
    NSlipPoints = k.size();

    // Slippage length at each time step
    SL = euclidean_distance(flatten1(CollPoints)[k[0]],
                            flatten1(CollPoints)[k[k.size()]], 0., 0.);

    // Adaptive time stepping
    double psi = 1.;

    if (SL - sl > psi) {
      DTnew = TimeStep / 10;
    } else {
      DTnew = TimeStep;
    }

    // Update the iterations
    ++iter;

    // Update the check criterion
    check = boole_mc(tau, sigma_n, fric, Cohes, NCollPoints);
  }

  // Cumulative slip at nodal points
  for (il::int_t m3 = 0; m3 < dd.size(); ++m3) {

    dd[m3] = dd[m3] + flatten1(incr_slip)[m3];
  }

  // Cumulative slip at collocation points
  il::Array<double> s;
  s = il::dot(from_edge_to_col(mesh.nelts(), 2), dd);

  // Update the friction coefficient for next time step
  fric = exp_friction(Peak_fric, Resid_fric, d_wf, s);

  // Update the dilatancy weakening coefficient for next time step
  dil = dilatancy(Init_dil, Incr_dil, d_wd, s);

  // Assign the desired outputs to structure's members
  res.dilatancy = dil;
  res.friction = fric;
  res.dt = DTnew;
  res.t = t;
  res.iter = iter;
  res.tot_stress_state = sigma_tot;
  res.slippagezone = SL;
  res.final_slip_points = k;
  res.P = press_prof;

  il::Array2D<double> ddout{mesh.nelts(), 2, 0.};

  for (il::int_t n3 = 0; n3 < ddout.size(0); ++n3) {

    for (il::int_t i = 0; i < ddout.size(1); ++i) {

      ddout(n3, i) = dd[mesh.conn(n3, i)];
    }
  }

  res.incr_d = ddout;

  return res;
}

///// OTHER UTILITIES /////

/// 1
il::Array<double> flatten1(il::Array2D<double> &Arr) {

  il::Array2D<int> A{Arr.size(0), 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * 2;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  il::Array<double> out{2 * Arr.size(0), 0.};

  for (il::int_t i = 0; i < Arr.size(0); ++i) {

    for (il::int_t j = 0; j < 2; ++j) {

      out[A(i, j)] = Arr(i, j);
    }
  }

  return out;
}

/// 2
il::Array<int> flatten2(il::Array2D<int> &Arr) {

  il::Array2D<int> A{Arr.size(0), 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * 2;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  il::Array<int> out{2 * Arr.size(0), 0.};

  for (il::int_t i = 0; i < Arr.size(0); ++i) {

    for (il::int_t j = 0; j < 2; ++j) {

      out[A(i, j)] = Arr(i, j);
    }
  }

  return out;
}

/// 3
int boole_mc(il::Array<double> tau, il::Array<double> sigma_n,
             il::Array<double> fric, const double Cohes,
             const il::int_t NCollPoints) {

  il::Array<int> T{tau.size(), 0};
  int out;

  for (il::int_t i = 0; i < NCollPoints; ++i) {

    if (tau[i] <= (Cohes + (fric[i] * sigma_n[i]))) {

      T[i] = 1;
    }
  }

  for (il::int_t m = 0; m < NCollPoints; ++m) {

    if (T[m] == 1) {

      out = out + T[m];
    }
  }

  return out;
}

/// 4
il::Array<int> select(il::Array<int> &arr) {

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
void set_submatrix_non_linear_system(il::Array2D<double> &A, int i0, int i1,
                                     const il::Array2D<double> &B) {

  IL_ASSERT(i0 + B.size(0) <= A.size(0));
  IL_ASSERT(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
}

/// 6
void set_subvector_non_linear_system(il::Array<double> &A, int i0,
                                     const il::Array<double> &B) {

  IL_ASSERT(i0 + B.size() <= A.size());

  for (int j0 = 0; j0 < B.size(); ++j0) {
    A[i0 + j0] = B[j0];
  }
}

/// 7
il::Array2D<double>
take_submatrix_non_linear_system(int i0, int i1, int j0, int j1,
                                 const il::Array2D<double> &A) {

  il::Array2D<double> sub{i1 - i0 + 1, j1 - j0 + 1, 0.};

  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      sub(i - i0, j - j0) = A(i, j);
    }
  }
  return sub;
}

/// 8
il::Array<double> take_subvector_non_linear_system(int i0, int i1,
                                                   const il::Array<double> &A) {

  il::Array<double> sub{i1 - i0 + 1, 0.};

  for (int i = i0; i <= i1; ++i) {

    sub[i - i0] = A[i];
  }
  return sub;
}

/// 9
double euclidean_distance(double x1, double x2, double y1, double y2) {

  double dist;

  double x = x1 - x2;
  double y = y1 - y2;

  dist = pow(x, 2) + pow(y, 2);
  dist = sqrt(dist);

  return dist;
}
}
