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

void time_incr(Mesh mesh, int p, double Cohes, const il::Array2D<double> &kmat,
               double Incr_dil, double d_wd, il::Array2D<double> rho,
               double Init_dil, double CompressFluid, double Visc,
               il::Array<double> S, int dof_dim, double Peak_fric,
               double Resid_fric, double d_wf, il::Array2D<double> Sigma0,
               il::Array<double> Amb_press, il::Array<double> Pinit,
               const std::string &Directory_results, il::Array<double> XColl,
               il::io_t);

int find(const il::Array<double> &arr, double_t seek, il::io_t);

double_t max_1d(const il::Array<double> &arr1D, il::io_t);
}

#endif // HFPX2D_TIMEINCR_H
