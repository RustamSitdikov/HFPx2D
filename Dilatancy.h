//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_DILATANCY_H
#define HFPX2D_DILATANCY_H

#include <il/linear_algebra.h>

il::Array<double> Dilatancy(const double Init_dil, const double Incr_dil, const double d_wd, il::Array<double> &d);

il::Array<double> DDilatancy(const double Incr_dil, const double d_wd, il::Array<double> &d);

#endif //HFPX2D_DILATANCY_H


