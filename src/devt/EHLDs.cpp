//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 03.03.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from the project
#include "EHLDs.h"
#include "src/Elasticity/AssemblyDDM.h"
#include "src/Mesh/DOF_Handles.h"
#include "FromEdgeToCol.h"

namespace hfp2d {

Results_solution_nonlinearsystem
EHLDs(Mesh mesh, il::Array2D<double> &kmatd, il::Array2D<double> &Npc,
      LayerParameters1 &layer_parameters1, LayerParameters2 &layer_parameters2,
      LayerParameters3 &layer_parameters3, il::Array<il::int_t> id_layers,
      Parameters_dilatancy &dilat_parameters,
      Parameters_fluid &fluid_parameters,
      Results_one_timeincrement &SolutionAtTj, il::Array<double> press_prof,
      il::Array<double> tot_slip, int dof_dim, int p, il::Array<double> cohes,
      il::Status &status, il::Norm norm, int inj_point, il::Array<double> S,
      il::Array<il::int_t> Dof_slip_coll, il::Array2D<double> Sigma0,
      il::Array2D<double> sigma_tot,
      hfp2d::simulation_parameters simulation_parameters,
      Parameters_permeability &permeab_parameters, il::io_t) {

  Results_solution_nonlinearsystem Results_iterations;

  //// FULLY IMPLICIT SOLUTION OF THE COUPLED PROBLEM ////
  // Initialization of the system BigA*BigX = BigB -> Just for shear DDs for now
  il::Array2D<double> BigA{
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1,
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1, 0};
  il::Array<double> BigB{
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1, 0};

  // Assembling just elasticity part of BigA
  hfp2d::set_submatrix_non_linear_system(BigA, 0, 0, kmatd);

  /// Solution of the non-linear system of equations via fixed point iterations

  // Initialization of Finite Volume matrices
  il::Array2D<double> Vd;
  il::Array2D<double> Vp;
  il::Array2D<double> L{mesh.nelts() + 1, mesh.nelts() + 1, 0};

  // Initialization of the while loop
  double errDd = 2;
  double errDp = 2;
  int j = 0;
  il::Array<double> BigX{BigA.size(0), 0};
  il::Array<double> Dd{tot_slip.size(), 0};
  il::Array<double> DiffDd{Dd.size(), 0};
  il::Array<double> Ddold = Dd;
  il::Array<double> Dp{press_prof.size(), 0};
  il::Array<double> Dpold{press_prof.size(), 0};
  il::Array<double> DiffDp{Dp.size(), 0};
  il::Array<double> Ddd{Dd.size(), 0};
  il::Array2D<double> incrdk{mesh.nelts(), 2, 0};
  il::Array<double> tot_slipk{SolutionAtTj.active_set_collpoints.size(), 0};
  il::Array<double> tot_slip_prev{SolutionAtTj.active_set_collpoints.size(), 0};

  while (j < simulation_parameters.itermax_nonlinsystem &&
         (errDd > simulation_parameters.tolerance ||
          errDp > simulation_parameters.tolerance)) {
    ++j;

    // Current total slip
    for (il::int_t i = 0; i < incrdk.size(0); ++i) {
      for (il::int_t l = 0; l < incrdk.size(1); ++l) {
        incrdk(i, l) =
            SolutionAtTj.d_tot[hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io)(i, l)] +
            Dd[hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io)(i, l)];
      }
    }

    // Mass matrix "Vd"
    Vd = hfp2d::build_vd_matrix_p1(
        mesh, dilat_parameters,
        hfp2d::dofhandle_dp(dof_dim, mesh.nelts(), p, il::io), fluid_parameters,
        incrdk, il::io);

    // Pressure matrix "P"
    Vp = hfp2d::build_vp_matrix_p1(mesh, dilat_parameters, fluid_parameters,
                                   incrdk, il::io);

    // Finite Difference matrix "L"
    L = hfp2d::build_l_matrix(mesh, incrdk, fluid_parameters, dilat_parameters,
                              SolutionAtTj.dt, permeab_parameters, il::io);

    // Nf -> diagonal matrix that contains the current friction coefficient of
    //       the slipping collocation points
    il::Array2D<double> Nf{SolutionAtTj.active_set_collpoints.size(),
                           SolutionAtTj.active_set_collpoints.size(), 0};
    il::Array2D<double> fetc_dg = hfp2d::from_edge_to_col_dg(
        dof_dim, hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io), il::io);

    il::Array<double> incrdk_flat = flatten1(incrdk, il::io);
    il::Array<double> incrdk_coll = il::dot(fetc_dg, incrdk_flat);
    il::Array<double> frick_coll = hfp2d::lin_friction(
        layer_parameters1, layer_parameters2, layer_parameters3, id_layers,
        hfp2d::dofhandle_dp(1, mesh.nelts(), p, il::io), incrdk_coll,
        il::io);
    for (il::int_t m2 = 0; m2 < Nf.size(0); ++m2) {
      Nf(m2, m2) = frick_coll[SolutionAtTj.active_set_collpoints[m2]];
    }

    /// Assembling the system ///

    // Matrix of coefficient

    il::Array2D<double> Npk = il::dot(Nf, Npc);
    for (il::int_t k2 = 0; k2 < SolutionAtTj.active_set_collpoints.size();
         ++k2) {
      for (il::int_t i = 0; i < mesh.nelts() + 1; ++i) {
        BigA(k2, SolutionAtTj.active_set_collpoints.size() + i) = Npk(k2, i);
      }
    }

    for (il::int_t k = 0; k < Vd.size(0); ++k) {
      for (il::int_t i = 0, q2 = 0;
           i < SolutionAtTj.active_set_collpoints.size(); ++i, q2 = q2 + 2) {
        BigA(SolutionAtTj.active_set_collpoints.size() + k, i) =
            Vd(k, Dof_slip_coll[q2]);
      }
    }

    for (il::int_t m = 0; m < Vp.size(0); ++m) {
      for (il::int_t i = 0; i < Vp.size(1); ++i) {
        BigA(SolutionAtTj.active_set_collpoints.size() + m,
             SolutionAtTj.active_set_collpoints.size() + i) =
            Vp(m, i) - L(m, i);
      }
    }

    // Right hand side vector

    for (il::int_t l2 = 0; l2 < tot_slipk.size(); ++l2) {
      tot_slipk[l2] = incrdk_flat[SolutionAtTj.active_set_collpoints[l2]];
    }

    auto tau_prev = il::dot(kmatd, tot_slipk);
    for (il::int_t n = 0, n1 = 1; n < SolutionAtTj.active_set_collpoints.size();
         ++n, n1 = n1 + 2) {
      BigB[n] =
          (cohes[SolutionAtTj.active_set_collpoints[n]] +
           frick_coll[SolutionAtTj.active_set_collpoints[n]] *
               (sigma_tot(SolutionAtTj.active_set_collpoints[n], 1) -
                SolutionAtTj.Pcm(SolutionAtTj.active_set_collpoints[n], 1))) -
          tau_prev[n] - Sigma0(SolutionAtTj.active_set_collpoints[n], 0);
    }

    auto Lp = il::dot(L, press_prof);
    for (il::int_t i1 = 0; i1 < L.size(0); ++i1) {
      BigB[SolutionAtTj.active_set_collpoints.size() + i1] =
          Lp[i1] + (SolutionAtTj.dt * S[i1]);
    }

    /// Boundary conditions ///

    for (il::int_t k1 = 0; k1 < BigA.size(1); ++k1) {
      BigA(SolutionAtTj.active_set_collpoints.size() + inj_point, k1) = 0;
    }

    for (il::int_t l1 = 0; l1 < BigA.size(0); ++l1) {
      BigA(l1, SolutionAtTj.active_set_collpoints.size() + inj_point) = 0;
    }

    BigA(SolutionAtTj.active_set_collpoints.size() + inj_point,
         SolutionAtTj.active_set_collpoints.size() + inj_point) = 1;

    BigB[SolutionAtTj.active_set_collpoints.size() + inj_point] = 0;

    /// Solve the system ///

    BigX = il::linear_solve(BigA, BigB, il::io, status);
    status.abort_on_error();

    ///  Under relaxation technique & updating ///

    for (il::int_t m1 = 0; m1 < SolutionAtTj.active_set_collpoints.size();
         ++m1) {
      Ddd[SolutionAtTj.active_set_collpoints[m1]] = BigX[m1];
    }

    for (il::int_t m1 = 0; m1 < Dd.size(); ++m1) {
      Dd[m1] = (1 - simulation_parameters.under_relax_param) * Ddold[m1] +
               simulation_parameters.under_relax_param * Ddd[m1];
    }

    for (il::int_t n1 = 0; n1 < Dp.size(); ++n1) {
      Dp[n1] = (1 - simulation_parameters.under_relax_param) * Dpold[n1] +
               simulation_parameters.under_relax_param *
                   BigX[SolutionAtTj.active_set_collpoints.size() + n1];
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

  Results_iterations.Dd = Ddd;
  Results_iterations.Dp = Dp;
  Results_iterations.errDd = errDd;
  Results_iterations.errDp = errDp;
  Results_iterations.new_tot_slip = incrdk;
  Results_iterations.Niterations = j;

  return Results_iterations;
}
}