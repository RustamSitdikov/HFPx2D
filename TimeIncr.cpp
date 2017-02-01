//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 31.01.17.
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
#include "TimeIncr.h"

namespace hfp2d {

void time_incr(Mesh mesh, const int p, const double Cohes, const int itermax,
               il::Array2D<double> &kmat, const double tolerance,
               const double Incr_dil, const double d_wd, il::Array2D<int> &Dof,
               il::Array2D<double> rho, const double Init_dil,
               const double CompressFluid, const double Visc,
               il::Array<double> S, const int InjPoint, const int dof_dim,
               const double Peak_fric, const double Resid_fric,
               const double d_wf, il::Array2D<double> Sigma0,
               il::Array2D<double> Amb_press, il::Array2D<double> Pinit) {

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

  /// Initialization ///
  Result SolutionAtTj;
  SolutionAtTj.slippagezone = 0.;
  il::Array<double> In1{NCollPoints, 0.};
  il::Array2D<double> In2{mesh.nelts(), 2, 0.};
  SolutionAtTj.friction = hfp2d::exp_friction(Peak_fric, Resid_fric, d_wf, In1);
  SolutionAtTj.dilatancy = hfp2d::dilatancy(Init_dil, Incr_dil, d_wd, In1);
  SolutionAtTj.incr_d = In2;
  SolutionAtTj.tot_stress_state = Sigma0;

  il::Array<double> InitP1{2 * mesh.nelts(), 0.};
  il::Array<double> InitP2{2 * mesh.nelts(), 0.};

  InitP1 = flatten1(Amb_press);
  InitP2 = flatten1(Pinit);

  il::Array<double> InitPP{2 * mesh.nelts(), 0.};
  for (il::int_t i = 0; i < InitPP.size(); ++i) {

    InitPP[i] = InitP1[i] + InitP2[i];
  }

  SolutionAtTj.P = InitPP;

  double t = 0.5;
  double tmax = 1.;
  double TimeStep = 0.005;

  while (t <= tmax) {

    SolutionAtTj =
        elhds(SolutionAtTj, mesh, p, Cohes, itermax, kmat, tolerance, Incr_dil,
              d_wd, Dof, rho, Init_dil, CompressFluid, TimeStep, Visc, S,
              InjPoint, dof_dim, t, Peak_fric, Resid_fric, d_wf);

    t = t + SolutionAtTj.t;
  }
}
}
