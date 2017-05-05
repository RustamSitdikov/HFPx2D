//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 02.03.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>

// Inclusion from Inside Loop library
#include <il/benchmark/tools/timer/Timer.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/norm.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "Dilatancy.h"
#include "EHLDs.h"
#include "FromEdgeToCol.h"

namespace hfp2d {

void MC_criterion(
    Mesh mesh, int p, il::Array<double> cohes, const il::Array2D<double> &kmat,
    LayerParameters1 &layer_parameters1, LayerParameters2 &layer_parameters2,
    LayerParameters3 &layer_parameters3, il::Array<il::int_t> id_layers,
    Parameters_dilatancy &dilat_parameters, Parameters_fluid &fluid_parameters,
    il::Array<double> S, int inj_point, int dof_dim, il::Array<double> XColl,
    il::Array2D<double> &Fetc, il::Array2D<double> Sigma0,
    hfp2d::simulation_parameters simulation_parameters,
    Parameters_permeability &permeab_parameters, il::io_t,
    Results_one_timeincrement &res) {

  // Vector of friction at the beginning of each time step
  il::Array<double> fric{res.friction.size(), 0};
  fric = res.friction;

  // Matrix of TOTAl stress at the beginning of each time step
  // Size -> Nelts x 2
  il::Array2D<double> sigma_tot{res.tot_stress_state.size(0),
                                res.tot_stress_state.size(1), 0};
  sigma_tot = res.tot_stress_state;

  // Vector of total slip at nodal points at the beginning of each time step
  il::Array<double> tot_slip{res.d_tot.size(), 0};
  tot_slip = res.d_tot;

  // Length of slippage zone at the beginning of each time step
  double sl;
  sl = res.slippagezone;

  // Vector of pressure profile at nodes at the beginning of each time step
  // Size -> Nelts + 1 -> {P_n1, P_n2, P_n3 ..}
  il::Array<double> press_prof{res.P.size(), 0};
  press_prof = res.P;

  // Matrix of pressure profile at collocation points at the beginning of each
  // time step. {{0,p_c1},{0,p_c2}..}
  il::Array2D<double> Pcm{res.Pcm.size(0), res.Pcm.size(1), 0};
  Pcm = res.Pcm;

  // Create matrix of EFFECTIVE stress for each collocation point {tau,sigma_n}
  // Size -> NCollPoints x 2
  il::Array2D<double> sigma_eff{2 * mesh.nelts(), 2, 0};
  for (int m = 0; m < sigma_eff.size(0); ++m) {
    for (int i = 0; i < sigma_eff.size(1); ++i) {
      sigma_eff(m, i) = sigma_tot(m, i) - Pcm(m, i);
    }
  }

  // Initialization of vector of DOFs of active set of collocation points
  il::Array<il::int_t> Dof_slip_coll{};

  // Declaration of variables
  double SL;
  il::Status status;
  il::Norm norm;
  norm = il::Norm::L1;

  // Initialization of the while loop
  Results_solution_nonlinearsystem res_nonlinearsystem;
  il::Array2D<double> sigma_tot_new{2 * mesh.nelts(), 2, 0};
  il::Array2D<double> sigma_eff_new = sigma_eff;
  il::Array<double> press_prof_new{press_prof.size(), 0};
  int iter = 1;
  MCcheck check;
  int cvg = 0;
  il::Array<double> s{2 * mesh.nelts(), 0};
  il::Array2D<double> Pcm_new{2 * mesh.nelts(), 2, 0};
  double crack_vel;

  // Iterate until each point falls below the Mohr-Coulomb failure line
  while (cvg != 1 && iter <= simulation_parameters.itermax_MCcriterion) {

    std::cout << "Iter to satisfy M-C criterion = " << iter << "\n";

    // Find the corresponding DOFs of the slipping collocation points
    for (il::int_t l3 = 0, k2 = 0; l3 < res.active_set_collpoints.size();
         ++l3) {
      Dof_slip_coll.resize(k2 + 2);
      Dof_slip_coll[k2] = 2 * res.active_set_collpoints[l3];
      Dof_slip_coll[k2 + 1] = 2 * res.active_set_collpoints[l3] + 1;
      k2 = k2 + 2;
    }

    // Select the elasticity matrix for just shear DDs of slipping DOFs
    il::Array2D<double> kmatd{res.active_set_collpoints.size(),
                              res.active_set_collpoints.size(), 0};
    for (il::int_t m = 0, k = 0; m < kmatd.size(0); ++m, k = k + 2) {
      for (il::int_t i = 0, qq = 0; i < kmatd.size(1); ++i, qq = qq + 2) {
        kmatd(m, i) = kmat(Dof_slip_coll[k], Dof_slip_coll[qq]);
      }
    }

    // Select just the "slipping" part of the matrix to switch from nodal
    // points to collocation points
    il::Array2D<double> Npc{res.active_set_collpoints.size(), press_prof.size(),
                            0};
    for (il::int_t l2 = 0, q3 = 1; l2 < Npc.size(0); ++l2, q3 = q3 + 2) {
      for (il::int_t i = 0; i < Npc.size(1); ++i) {
        Npc(l2, i) = Fetc(Dof_slip_coll[q3], i);
      }
    }

    // Call the Elasto-Hydrodynamic Lubrication Dilatant solver for the solution
    // of the fully implicit coupled problem
    res_nonlinearsystem = hfp2d::EHLDs(
        mesh, kmatd, Npc, layer_parameters1, layer_parameters2,
        layer_parameters3, id_layers, dilat_parameters, fluid_parameters, res,
        press_prof, tot_slip, dof_dim, p, cohes, status, norm, inj_point, S,
        Dof_slip_coll, Sigma0, sigma_tot, simulation_parameters,
        permeab_parameters, il::io);

    std::cout << "Total N. of iterations for solving non-linear system of "
                 "equations = "
              << res_nonlinearsystem.Niterations << " || "
              << " err on Dd: " << res_nonlinearsystem.errDd << " & "
              << " err on Dp: " << res_nonlinearsystem.errDp << "\n"
              << "\n";

    // Update the pore pressure profile at nodal points
    for (il::int_t l2 = 0; l2 < press_prof.size(); ++l2) {
      press_prof_new[l2] = press_prof[l2] + res_nonlinearsystem.Dp[l2];
    }

    // Update the pore pressure profile at collocation points
    auto press_prof_coll_new = il::dot(
        hfp2d::from_edge_to_col_cg(
            dof_dim, hfp2d::dofhandle_dp(dof_dim, mesh.nelts(), p, il::io),
            hfp2d::dofhandle_cg(dof_dim, mesh.nelts(), il::io), il::io),
        press_prof_new);

    for (il::int_t l = 0, k = 1; l < Pcm.size(0); ++l) {
      Pcm_new(l, 1) = press_prof_coll_new[k];
      k = k + 2;
    }

    // Calculate the increment of shear stress due to the increment of slip
    il::Array2D<double> IncrT{2 * mesh.nelts(), 2, 0};
    for (il::int_t k2 = 0, k = 0; k2 < IncrT.size(0); ++k2, k = k + 2) {
      for (il::int_t i = 0, q5 = 0; i < 2 * mesh.nelts(); ++i, q5 = q5 + 2) {
        IncrT(k2, 0) = IncrT(k2, 0) + kmat(k, q5) * res_nonlinearsystem.Dd[i];
      }
    }

    // Update the CURRENT effective stress state
    for (il::int_t n2 = 0; n2 < sigma_eff.size(0); ++n2) {
      for (il::int_t i = 0; i < sigma_eff.size(1); ++i) {
        sigma_eff_new(n2, i) = sigma_tot(n2, i) - Pcm_new(n2, i) + IncrT(n2, i);
      }
    }

    // Cumulative slip at collocation points
    s = il::dot(hfp2d::from_edge_to_col_dg(
                    dof_dim, hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io),
                    il::io),
                flatten1(res_nonlinearsystem.new_tot_slip, il::io));

    // Update friction coefficient
    fric = hfp2d::lin_friction(
        layer_parameters1, layer_parameters2, layer_parameters3, id_layers,
        hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io), s, il::io);

    // Force the slipping nodes at the first iteration to stay in the MC
    // line in the next iterations (because all the slipping nodes in one time
    // increment displace simultaneously)
    il::Array2D<double> MC{res.active_set_collpoints.size(), 2, 0};
    for (il::int_t m2 = 0; m2 < MC.size(0); ++m2) {
      MC(m2, 0) = fabs(cohes[res.active_set_collpoints[m2]] +
                       fric[res.active_set_collpoints[m2]] *
                           (sigma_tot(res.active_set_collpoints[m2], 1) -
                            Pcm_new(res.active_set_collpoints[m2], 1)));
      MC(m2, 1) = fabs(sigma_tot(res.active_set_collpoints[m2], 1) -
                       Pcm_new(res.active_set_collpoints[m2], 1));
    }

    for (il::int_t k3 = 0; k3 < MC.size(0); ++k3) {
      for (il::int_t i = 0; i < sigma_eff.size(1); ++i) {
        sigma_eff_new(res.active_set_collpoints[k3], i) = MC(k3, i);
      }
    }

    // Total stress state (current) at collocation points for next time
    // increment
    for (il::int_t j3 = 0; j3 < sigma_tot.size(0); ++j3) {
      for (il::int_t i = 0; i < sigma_tot.size(1); ++i) {
        sigma_tot_new(j3, i) = sigma_eff_new(j3, i) + Pcm_new(j3, i);
      }
    }

    // Update the check criterion of the external while loop with current state
    check = boole_mc(sigma_eff_new, fric, cohes, il::io);
    // If the number of collocation points the satisfy the MC criterion is equal
    // to the total number of collocation points, then converged is reached (cvg
    // = 1). Else, append the new active collocation points to the past active
    // set of collocation points.

    // the capacity of the vector active_set_collpoints is at maximum equal to
    // the number of collocation points
    res.active_set_collpoints.reserve(2 * mesh.nelts());

    if (check.Ncollpoint_satisfMC == 2 * mesh.nelts()) {

      cvg = 1;

    } else {

      for (il::int_t i = 0; i < check.CollPoint_notsatisfMC.size(); ++i) {

        res.active_set_collpoints.append(check.CollPoint_notsatisfMC[i]);
      }

      // Sort the current vector of active collocation points
      sort_ascending_order(res.active_set_collpoints, il::io);
    }

    // Calculate current slippage length
    if (res.active_set_collpoints.size() == 0) {
      SL = 0;
    } else {
      SL = euclidean_distance(
          XColl[res.active_set_collpoints[0]],
          XColl[res.active_set_collpoints[res.active_set_collpoints.size() -
                                          1]],
          0, 0, il::io);
    }

    // Calculate the current shear crack velocity
    if (SL == 0) {
      crack_vel = 0;
    } else {
      crack_vel = (SL - sl) / res.dt;
    }

    // Assign the desired outputs to structure's members
    res.friction = fric;
    res.iter = iter;
    res.d_tot = flatten1(res_nonlinearsystem.new_tot_slip, il::io);
    res.tot_stress_state = sigma_tot_new;
    res.slippagezone = SL;
    res.P = press_prof_new;
    res.Pcm = Pcm_new;
    res.crack_velocity = crack_vel;

    // Update the iterations
    ++iter;
  }
}

///// OTHER UTILITIES /////

/// 1
// This routine turns a 2D array of double precision values into a vector of
// double precision values.
// Input:
//   - Arr -> 2D array that we want to "flat"
il::Array<double> flatten1(const il::Array2D<double> &Arr, il::io_t) {

  il::Array2D<il::int_t> A{Arr.size(0), 2, 0};

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

  il::Array2D<il::int_t> A{Arr.size(0), 2, 0};

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
//  - sigma_eff -> matrix of effective stress
//  - fric -> vector of friction coefficient at each collocation point
//  - Cohes -> cohesion (double precision value)
//  - NCollPoints -> number of collocation points
MCcheck boole_mc(il::Array2D<double> &sigma_eff, il::Array<double> fric,
                 il::Array<double> cohes, il::io_t) {

  MCcheck result_MCcheck;
  il::Array<int> T{sigma_eff.size(0), 0};

  for (il::int_t i = 0; i < T.size(); ++i) {
    if (sigma_eff(i, 0) <= (cohes[i] + (fric[i] * sigma_eff(i, 1)))) {

      T[i] = 1;
    }
  }

  for (il::int_t m = 0, k = 0; m < T.size(); ++m) {
    if (T[m] == 1) {

      result_MCcheck.Ncollpoint_satisfMC =
          result_MCcheck.Ncollpoint_satisfMC + T[m];

    } else {

      result_MCcheck.CollPoint_notsatisfMC.resize(k + 1);
      result_MCcheck.CollPoint_notsatisfMC[k] = m;
      k = k + 1;
    }
  }

  return result_MCcheck;
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
// Function that return the index ( IndRow,IndColumn ) of a given value ("seek")
// in an array2D.
// arr2D -> array2D in which we want to find the index of  given value
// seek -> value for which we want to find out the index
// It return an array that contain the {N.row, N.col} of the seek value

il::Array<int> find_2d(const il::Array2D<double> &arr2D, double_t seek,
                       il::io_t) {

  il::Array<int> outp{2};

  for (int i = 0; i < arr2D.size(0); ++i) {

    for (int j = 0; j < arr2D.size(1); ++j) {

      if (arr2D(i, j) == seek)
        outp[0] = i, outp[1] = j;
    }
  }

  return outp;
}

/// 10

il::Array<il::int_t> delete_duplicates(const il::Array<il::int_t> &arr,
                                       il::io_t) {

  il::Array<il::int_t> out{0};

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

/// 11

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

/// 12

void sort_ascending_order(il::Array<il::int_t> &arr, il::io_t) {

  il::int_t deposito;

  for (il::int_t m3 = 0; m3 < arr.size() - 1; ++m3) {
    for (il::int_t i = m3 + 1; i < arr.size(); ++i) {

      if (arr[m3] > arr[i]) {

        deposito = arr[m3];
        arr[m3] = arr[i];
        arr[i] = deposito;
      }
    }
  }
}
}
