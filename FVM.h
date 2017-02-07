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

il::Array<double> average(const il::Array2D<double> &d);

il::Array<double> quarter(const il::Array2D<double> &d);

il::Array2D<int> position_2d_array(const il::Array2D<int> &arr2D, int seek);

il::Array2D<int> search(const il::Array2D<int> &matrix, int x);

il::Array<int> row_selection(il::Array2D<int> &arr, il::int_t idx);

il::Array<double> shear_conductivities_newtonian(
    const double Visc, Mesh mesh, il::Array2D<double> rho, il::Array2D<double> &d,
    const double Incr_dil, const double d_wd, const double Init_dil);

il::Array2D<double> build_l_matrix(Mesh mesh, il::Array2D<double> &d,
                                   il::Array2D<double> &rho, const double Visc,
                                   const double Incr_dil, const double d_wd,
                                   const double Init_dil,
                                   const double &TimeStep);

il::Array2D<double> build_vp_matrix_p1(Mesh mesh, const double Incr_dil,
                                       const double Init_dil,
                                       const double CompressFluid,
                                       il::Array2D<double> &d,
                                       const double d_wd);

il::Array2D<double> build_vd_matrix_p1(Mesh mesh, const double Incr_dil,
                                       const double d_wd, il::Array2D<int> Dof,
                                       il::Array2D<double> rho,
                                       il::Array2D<double> &d);

il::Array2D<double> build_vp_matrix_p0(Mesh mesh, const double Incr_dil,
                                       const double Init_dil,
                                       const double CompressFluid,
                                       il::Array<double> &d,
                                       const double d_wd);

il::Array2D<double> build_vd_matrix_p0(Mesh mesh, const double Incr_dil,
                                       const double d_wd, il::Array2D<int> &Dof,
                                       il::Array2D<double> rho,
                                       il::Array<double> &d);
}

#endif // HFPX2D_FVM_H
