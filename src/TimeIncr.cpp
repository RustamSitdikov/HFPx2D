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
#include <il/Timer.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/norm.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "Dilatancy.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
#include "MC_criterion.h"
#include "Output_results.h"

namespace hfp2d {

void time_incr(
    int inj_point, il::int_t NCollPoints, Mesh mesh, int p,
    il::Array<double> cohes, const il::Array2D<double> &kmat,
    LayerParameters1 layer_parameters1, LayerParameters2 layer_parameters2,
    LayerParameters3 layer_parameters3, il::Array<il::int_t> id_layers,
    Parameters_dilatancy dilat_parameters, Parameters_fluid &fluid_parameters,
    il::Array<double> S, int dof_dim, il::Array2D<double> Sigma0,
    il::Array<double> Amb_press, il::Array<double> Pinit,
    const std::string &Directory_results, il::Array<double> XColl,
    il::Array2D<double> &Fetc, double h,
    simulation_parameters simulation_parameters, double kf, il::io_t) {

  /// Initialization ///

  // Initialization of the structure of module "MC_criterion"
  Results_one_timeincrement SolutionAtTj;

  // Initialization of matrix of TOTAL stress at time t_0 (before the injection)
  // {{tau_1,sigma_n1},{tau_2, sigma_n2},{tau3, sigma_n3} ..}
  SolutionAtTj.tot_stress_state = Sigma0;

  // Initialization of vector total slip at nodal points at time t_0 (before the
  // injection)
  // Remember: piecewise linear DDs
  il::Array<double> no_slip{NCollPoints, 0};
  SolutionAtTj.d_tot = no_slip;

  // Initialization of friction vector at time t_0 (before the injection)
  SolutionAtTj.friction = hfp2d::lin_friction(
      layer_parameters1, layer_parameters2, layer_parameters3, id_layers,
      hfp2d::dofhandle_dg(dof_dim, mesh.nelts(), il::io), no_slip, il::io);

  // Initialization of slippage length (before the injection)
  SolutionAtTj.slippagezone = 0;

  // Initialization of pressure profile at nodal points & collocation points at
  // time t_0plus.
  //  -> sum of the ambient pressure profile + initial pore press perturbation
  // (the latter is needed to activate the shear crack)
  il::Array<double> Pin{mesh.nelts() + 1, 0};
  for (il::int_t j = 0; j < Pin.size(); ++j) {
    Pin[j] = Pinit[j] + Amb_press[j];
  }
  SolutionAtTj.P = Pin;

  il::Array2D<double> Pcm{2 * mesh.nelts(), 2, 0};
  il::Array2D<double> fetc_cg = hfp2d::from_edge_to_col_cg(
      dof_dim, hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io),
      hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io), il::io);
  il::Array<double> pcoll = il::dot(fetc_cg, Pin);

  for (il::int_t l = 0, k = 1; l < Pcm.size(0); ++l) {
    Pcm(l, 1) = pcoll[k];
    k = k + 2;
  }
  SolutionAtTj.Pcm = Pcm;

  // Initialization of active set of collocation points at t_0plus
  for (il::int_t i = 0, k = 0; i < NCollPoints; ++i) {
    if (SolutionAtTj.tot_stress_state(i, 0) >=
        cohes[i] +
            SolutionAtTj.friction[i] * (SolutionAtTj.tot_stress_state(i, 1) -
                                        SolutionAtTj.Pcm(i, 1))) {
      SolutionAtTj.active_set_collpoints.resize(k + 1);
      SolutionAtTj.active_set_collpoints[k] = i;
      k = k + 1;
    }
  }

  // Initial time step
  double TimeStep = simulation_parameters.TimeStep_min;
  SolutionAtTj.dt = TimeStep;

  // Initialization of time
  double t;
  t = simulation_parameters.t_0plus;

  // Parameter for adaptive time stepping
  double psi = 1;

  while (
      t <= simulation_parameters.t_max &&
      SolutionAtTj.slippagezone <=
          euclidean_distance(XColl[0], XColl[XColl.size() - 1], 0, 0, il::io)) {

    std::cout << "******** Time increment *******"
              << " "
              << "t = " << t << "\n";

    hfp2d::MC_criterion(mesh, p, cohes, kmat, layer_parameters1,
                        layer_parameters2, layer_parameters3, id_layers,
                        dilat_parameters, fluid_parameters, S, inj_point,
                        dof_dim, XColl, Fetc, Sigma0, simulation_parameters, kf,
                        il::io, SolutionAtTj);

    // Adaptive time stepping
    if (psi * (h / SolutionAtTj.crack_velocity) >=
        simulation_parameters.TimeStep_max) {
      SolutionAtTj.dt = simulation_parameters.TimeStep_max;
    } else if (psi * (h / SolutionAtTj.crack_velocity) <=
               simulation_parameters.TimeStep_min) {
      SolutionAtTj.dt = simulation_parameters.TimeStep_min;
    } else {
      SolutionAtTj.dt = psi * (h / SolutionAtTj.crack_velocity);
    }

    // Export to get a different output file per each iteration
    hfp2d::export_results(SolutionAtTj, t, Directory_results,
                          std::string{"Test"} + std::to_string(t) +
                              std::string{".txt"});

    t = t + SolutionAtTj.dt;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the index of a given value ("seek") in a vector.
// arr -> 1D array in which we want to find the index of a given value
// seek ->  value for which we want to find out the index
int find(const il::Array<double> &arr, double_t seek, il::io_t) {

  for (int i = 0; i < arr.size(); ++i) {
    if (arr[i] == seek)
      return i;
  }

  return -1;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the max value in an vector.
// arr1D -> vector in which we want to find the max
double_t max_1d(const il::Array<double> &arr1D, il::io_t) {

  double_t max;
  max = arr1D[0];

  for (il::int_t i = 0; i < arr1D.size(); ++i) {

    if (max < arr1D[i]) {

      max = arr1D[i];
    }
  }

  return max;
}

////////////////////////////////////////////////////////////////////////////////
// This function calculates the 2D euclidean distance between two points
// {x1,y1}, {x2,y2}
// Input -> x- y- coordinates of the two points
// Output -> double precision values that represents the euclidean distance
//           between the two aforementioned points
double euclidean_distance(double x1, double x2, double y1, double y2,
                          il::io_t) {

  double dist;

  double x = x1 - x2;
  double y = y1 - y2;

  dist = (x * x) + (y * y);
  dist = sqrt(dist);

  return dist;
}
}
