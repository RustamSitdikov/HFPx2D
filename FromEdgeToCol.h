//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FROMEDGETOCOL_H
#define HFPX2D_FROMEDGETOCOL_H

// Inclusion from standard library
#include <cmath>

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

namespace hfp2d {

il::Array2D<double> from_edge_to_col_dg_full2d(int dof_dim,
                                               il::Array2D<int> Dof, il::io_t);

il::Array2D<double> from_edge_to_col_dg(int dof_dim,
                                        il::Array2D<int> Dofw, il::io_t);

il::Array2D<double> from_edge_to_col_cg(int dof_dim,
                                        il::Array2D<int> Dof,
                                        il::Array2D<int> Dofp, il::io_t);
}

#endif // HFPX2D_FROMEDGETOCOL_H
