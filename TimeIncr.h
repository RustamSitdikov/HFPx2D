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

void time_incr(Mesh mesh, const int p, const double Cohes, const int itermax,
               il::Array2D<double> &kmat, const double tolerance,
               const double Incr_dil, const double d_wd, il::Array2D<int> &Dof,
               il::Array2D<double> rho, const double Init_dil,
               const double CompressFluid, const double Visc,
               il::Array<double> S, const int InjPoint, const int dof_dim,
               const double Peak_fric, const double Resid_fric,
               const double d_wf, il::Array2D<double> Sigma0, il::Array2D<double> Amb_press, il::Array2D<double> Pinit);
}

#endif // HFPX2D_TIMEINCR_H
