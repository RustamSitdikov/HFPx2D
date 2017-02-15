//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FRICTION_H
#define HFPX2D_FRICTION_H

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

namespace hfp2d {

il::Array<double> exp_friction(double Peak_fric, double Resid_fric, double d_wf,
                               const il::Array<double> &d, il::io_t);

il::Array<double> lin_friction(double Peak_fric, double Resid_fric, double d_wf,
                               const il::Array<double> &d, il::io_t);
}

#endif // HFPX2D_FRICTION_H
