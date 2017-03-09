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

void time_incr(double t_0plus, int inj_point, il::int_t NCollPoints, Mesh mesh,
               int p, il::Array<double> cohes, const il::Array2D<double> &kmat,
               Parameters_friction fric_parameters,
               Parameters_dilatancy dilat_parameters,
               Parameters_fluid &fluid_parameters, il::Array<double> S,
               int dof_dim, il::Array2D<double> Sigma0,
               il::Array<double> Amb_press, il::Array<double> Pinit,
               const std::string &Directory_results, il::Array<double> XColl,
               il::Array2D<double> &Fetc, double h, double TimeStep_max,
               double TimeStep_min, il::io_t);

int find(const il::Array<double> &arr, double_t seek, il::io_t);

double_t max_1d(const il::Array<double> &arr1D, il::io_t);

double euclidean_distance(double x1, double x2, double y1, double y2, il::io_t);
}

#endif // HFPX2D_TIMEINCR_H
