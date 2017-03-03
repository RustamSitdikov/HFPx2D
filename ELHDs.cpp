//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 03.03.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "ELHDs.h"

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "FromEdgeToCol.h"
#include "TimeIncr.h"

namespace hfp2d {

Results_solution_nonlinearsystem
ELHDs(Mesh mesh, il::Array2D<double> &kmatd, il::Array2D<double> &Npc,
      Parameters_friction fric_parameters,
      Parameters_dilatancy dilat_parameters, Parameters_fluid fluid_parameters,
      Results_one_timeincrement &SolutionAtTj, il::Array2D<double> sigma_tot,
      il::Array<double> press_prof, il::Array2D<double> Pcm,
      il::Array<double> tot_slip, int itermax, int dof_dim, int p,
      il::Array<double> cohes, il::Status status, il::Norm norm, int inj_point,
      il::Array<double> S, il::Array<int> Dof_slip_coll, il::io_t) {

  Results_solution_nonlinearsystem Results_iterations;

  //// FULLY IMPLICIT SOLUTION OF THE COUPLED PROBLEM ////
  // Initialization of the system BigA*BigX = BigB -> Just for shear DDs!!
  il::Array2D<double> BigA{
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1,
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1, 0};
  il::Array<double> BigB{
      SolutionAtTj.active_set_collpoints.size() + mesh.nelts() + 1, 0};

  // Assembling just elasticity part of BigA
  hfp2d::set_submatrix_non_linear_system(BigA, 0, 0, kmatd);

  /// Solution of the non-linear system of equations by using
  /// fixed point iterations

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
  il::Array<double> Dd{tot_slip.size(), 0};
  il::Array<double> DiffDd{Dd.size(), 0};
  il::Array<double> Ddold = Dd;
  il::Array<double> Dp{press_prof.size(), 0};
  il::Array<double> Dpold{press_prof.size(), 0};
  il::Array<double> DiffDp{Dp.size(), 0};
  il::Array<double> Ddd{Dd.size(), 0};

  while (j < itermax && (errDd > tolerance || errDp > tolerance)) {
    ++j;

    for (il::int_t i = 0; i < incrdk.size(0); ++i) {
      for (il::int_t l = 0; l < incrdk.size(1); ++l) {
        incrdk(i, l) =
            tot_slip[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io)(i, l)] +
            Dd[hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io)(i, l)];
      }
    }

    // Mass matrix "Vd"
    Vd = hfp2d::build_vd_matrix_p1(
        mesh, dilat_parameters,
        hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io),
        fluid_parameters, incrdk, il::io);

    // Pressure matrix "P"
    Vp = hfp2d::build_vp_matrix_p1(mesh, dilat_parameters, fluid_parameters,
                                   incrdk, il::io);
    // Finite Difference matrix "L"
    L = hfp2d::build_l_matrix(mesh, incrdk, fluid_parameters, dilat_parameters,
                              SolutionAtTj.dt, il::io);

    il::Array2D<double> Nf{SolutionAtTj.active_set_collpoints.size(),
                           SolutionAtTj.active_set_collpoints.size(), 0};
    for (il::int_t m2 = 0; m2 < Nf.size(0); ++m2) {
      Nf(m2, m2) = hfp2d::lin_friction(
          fric_parameters,
          il::dot(hfp2d::from_edge_to_col_dg(
                      dof_dim,
                      hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io),
                      il::io),
                  flatten1(incrdk, il::io)),
          il::io)[SolutionAtTj.active_set_collpoints[m2]];
    }

    /// Assembling the system ///

    // Matrix of coefficient
    for (il::int_t k2 = 0; k2 < SolutionAtTj.active_set_collpoints.size();
         ++k2) {
      for (il::int_t i = 0; i < mesh.nelts() + 1; ++i) {
        BigA(k2, SolutionAtTj.active_set_collpoints.size() + i) =
            il::dot(Nf, Npc)(k2, i);
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
    for (il::int_t n = 0, n1 = 1; n < SolutionAtTj.active_set_collpoints.size();
         ++n, n1 = n1 + 2) {
      BigB[n] =
          (cohes[SolutionAtTj.active_set_collpoints[n]] +
           hfp2d::lin_friction(
               fric_parameters,
               il::dot(hfp2d::from_edge_to_col_dg(
                           dof_dim,
                           hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io),
                           il::io),
                       flatten1(incrdk, il::io)),
               il::io)[SolutionAtTj.active_set_collpoints[n]] *
               (sigma_tot(SolutionAtTj.active_set_collpoints[n], 1) -
                Pcm(SolutionAtTj.active_set_collpoints[n], 1))) -
          (sigma_tot(SolutionAtTj.active_set_collpoints[n], 0));

      std::cout << BigB[n] << " ";
    }

    std::cout << "\n";

    for (il::int_t i1 = 0; i1 < L.size(0); ++i1) {
      BigB[SolutionAtTj.active_set_collpoints.size() + i1] =
          il::dot(L, press_prof)[i1] + (SolutionAtTj.dt * S[i1]);
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
      Dd[m1] = +(1 - betarela) * Ddold[m1] + betarela * Ddd[m1];
    }

    for (il::int_t n1 = 0; n1 < Dp.size(); ++n1) {
      Dp[n1] = (1 - betarela) * Dpold[n1] +
               betarela * BigX[SolutionAtTj.active_set_collpoints.size() + n1];
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

  Results_iterations.Dd = Dd;
  Results_iterations.Dp = Dp;
  Results_iterations.errDd = errDd;
  Results_iterations.errDp = errDp;
  Results_iterations.new_tot_slip = incrdk;
  Results_iterations.Niterations = j;

  return Results_iterations;
}
}