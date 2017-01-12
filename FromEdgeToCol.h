//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_FROMEDGETOCOL_H
#define HFPX2D_FROMEDGETOCOL_H

#include <il/linear_algebra.h>
#include <cmath>

il::Array<double > FromEdgeToCol(il::Array<double> &d_edge, il::Array2D<int> id, int Nelts, int dof_dim);

#endif //HFPX2D_FROMEDGETOCOL_H
