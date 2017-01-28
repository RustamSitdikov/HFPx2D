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

il::Array<double> from_edge_to_col(il::Array<double> &d_edge, const int Nelts,
                                const int dof_dim);
}

#endif // HFPX2D_FROMEDGETOCOL_H
