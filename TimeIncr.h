//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 31.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_TIMEINCR_H
#define HFPX2D_TIMEINCR_H

namespace hfp2d {

struct simulation_parameters {

  // initial time (starting time)
  double t_0plus;

  // Maximum time
  double t_max;

  // The time stepping is adapted depending on the current shear crack velocity;
  // for this reason we need to constrain the time stepping in order to avoid a
  // very small or very large time step

  // Minimum time step
  double TimeStep_min;

  // Maximum time step
  double TimeStep_max;

  // Maximum number of iterations to satisfy MC criterion
  il::int_t itermax_MCcriterion;

  // Maximum number of iterations to solve non linear system of equations via
  // fixed point iterations combined with under relaxation technique
  il::int_t itermax_nonlinsystem;

  // Tolerance for errors in the solution of the non linear system of equations
  double tolerance;

  // Coefficient for under relaxation technique
  double under_relax_param;
};

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
    simulation_parameters simulation_parameters, double kf, il::io_t);

int find(const il::Array<double> &arr, double_t seek, il::io_t);

double_t max_1d(const il::Array<double> &arr1D, il::io_t);

double euclidean_distance(double x1, double x2, double y1, double y2, il::io_t);
}

#endif // HFPX2D_TIMEINCR_H
