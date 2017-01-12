//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_FVM_H
#define HFPX2D_FVM_H

#include <il/linear_algebra.h>
#include <cmath>

il::Array<double> Average(il::Array2D<double> &d);

il::Array<double> Quarter(il::Array2D<double> &d);

#endif //HFPX2D_FVM_H
