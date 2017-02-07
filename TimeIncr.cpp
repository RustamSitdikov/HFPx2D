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

void time_incr(Mesh mesh, const int p, const double Cohes,
               il::Array2D<double> &kmat, const double Incr_dil,
               const double d_wd, il::Array2D<double> rho,
               const double Init_dil, const double CompressFluid,
               const double Visc, il::Array<double> S, const int dof_dim,
               const double Peak_fric, const double Resid_fric,
               const double d_wf, il::Array2D<double> Sigma0,
               il::Array<double> Amb_press, il::Array<double> Pinit) {

  // Total numbers of collocation points
  il::int_t NCollPoints = 2 * mesh.nelts();

  /// Initialization ///

  // Initialization of the structure of module "ELHDs"
  Result SolutionAtTj;

  // Initialization of pore pressure profile at nodal points
  il::Array<double> Pin{mesh.nelts() + 1, 0.};
  for (il::int_t j = 0; j < Pin.size(); ++j) {

    Pin[j] = Pinit[j] + Amb_press[j];
  }

  SolutionAtTj.P = Pin;

  // Injection point
  il::int_t InjPoint;
  InjPoint = hfp2d::find(SolutionAtTj.P, max_1d(SolutionAtTj.P));

  // Initialization of slippage length (No slip condition before fluid
  // injection)
  SolutionAtTj.slippagezone = 0.;

  // Initialization of friction vector
  il::Array<double> In1{NCollPoints, 0.};
  SolutionAtTj.friction = hfp2d::exp_friction(Peak_fric, Resid_fric, d_wf, In1);

  // Initialization of dilatancy vector
  SolutionAtTj.dilatancy = hfp2d::dilatancy(Init_dil, Incr_dil, d_wd, In1);

  // Initialization of vector total slip at nodal points
  // Remember: piecewise linear shear DDs
  il::Array<double> in_incr_d{2 * mesh.nelts(), 0.};
  SolutionAtTj.incr_d = in_incr_d;

  // Initialization of matrix of total stress
  // {{tau_1,sigma_n1},{tau_2, sigma_n2},{tau3, sigma_n3} ..}
  SolutionAtTj.tot_stress_state = Sigma0;

  double t = 0.005;
  double tmax = 0.01;
  double TimeStep = 0.005;

  while (t <= tmax) {

    hfp2d::elhds(SolutionAtTj, mesh, p, Cohes, kmat, Incr_dil, d_wd, rho,
                 Init_dil, CompressFluid, TimeStep, Visc, S, InjPoint, dof_dim,
                 Peak_fric, Resid_fric, d_wf);

    t = t + SolutionAtTj.dt;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the index of a given value ("seek") in a vector.
// arr -> 1D array in which we want to find the index of a given value
// seek ->  value for which we want to find out the index

il::int_t find(il::Array<double> arr, double_t seek) {

  for (il::int_t i = 0; i < arr.size(); ++i) {

    if (arr[i] == seek)
      return i;
  }

  return -1;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the max value in an vector.
// arr1D -> vector in which we want to find the max

double_t max_1d(il::Array<double> &arr1D) {

  double_t max;
  max = arr1D[0];

  for (il::int_t i = 0; i < arr1D.size(); ++i) {

    if (max < arr1D[i]) {

      max = arr1D[i];
    }
  }

  return max;
}
}
