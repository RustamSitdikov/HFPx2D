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

il::Array<double> Average(const il::Array2D<double> &d);

il::Array<double> Quarter(const il::Array2D<double> &d);

il::Array2D<int> Position_2DArray(const il::Array2D<int> &arr2D, int seek);

il::Array2D<int> Search(const il::Array2D<int> &matrix, int x);

il::Array<int> RowSelection(il::Array2D<int> &arr, il::int_t idx);

il::Array<double>
ShearConductivitiesNewtonian(const int Visc, Mesh mesh, il::Array2D<double> rho,
                             il::Array2D<double> &d, const double Incr_dil,
                             const double d_wd, const double Init_dil);

il::Array2D<double> BuildVpMatrix(Mesh mesh, const double Incr_dil,
                                  const double Init_dil,
                                  const double CompressFluid,
                                  il::Array2D<double> &d, const double d_wd);

il::Array2D<double> BuildVdMatrix(Mesh mesh, const double Incr_dil,
                                  const double d_wd, il::Array2D<int> &Dof,
                                  il::Array2D<double> rho,
                                  il::Array2D<double> &d);
}

#endif // HFPX2D_FVM_H
