//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FVM_H
#define HFPX2D_FVM_H

// Inclusion from standard library
#include <cmath>

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

// Inclusion from the project
#include "Mesh.h"

namespace hfp2d {

il::Array<double> average(const il::Array2D<double> &d, il::io_t);

il::Array<double> quarter(const il::Array2D<double> &d, il::io_t);

il::Array2D<int> position_2d_array(const il::Array2D<int> &arr2D, int seek,
                                   il::io_t);

il::Array2D<int> search(const il::Array2D<int> &matrix, int x, il::io_t);

il::Array<int> row_selection(const il::Array2D<int> &arr, il::int_t idx,
                             il::io_t);

il::Array<double> shear_conductivities_newtonian(double Visc, Mesh mesh,
                                                 il::Array2D<double> rho,
                                                 const il::Array2D<double> &d,
                                                 double Incr_dil, double d_wd,
                                                 double Init_dil, il::io_t);

il::Array2D<double> build_l_matrix(Mesh mesh, const il::Array2D<double> &d,
                                   const il::Array2D<double> &rho, double Visc,
                                   double Incr_dil, double d_wd,
                                   double Init_dil, const double &TimeStep,
                                   il::io_t);

il::Array2D<double> build_vp_matrix_p1(Mesh mesh, double Incr_dil,
                                       double Init_dil, double CompressFluid,
                                       const il::Array2D<double> &d,
                                       double d_wd, il::io_t);

il::Array2D<double> build_vd_matrix_p1(Mesh mesh, double Incr_dil, double d_wd,
                                       il::Array2D<int> Dof,
                                       il::Array2D<double> rho,
                                       const il::Array2D<double> &d, il::io_t);

il::Array2D<double> build_vp_matrix_p0(Mesh mesh, double Incr_dil,
                                       double Init_dil, double CompressFluid,
                                       const il::Array<double> &d, double d_wd,
                                       il::io_t);

il::Array2D<double> build_vd_matrix_p0(Mesh mesh, double Incr_dil, double d_wd,
                                       il::Array2D<int> &Dof,
                                       il::Array2D<double> rho,
                                       il::Array<double> &d, il::io_t);
}

#endif // HFPX2D_FVM_H
