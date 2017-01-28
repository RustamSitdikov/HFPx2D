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

il::Array<double> exp_friction(const double Peak_fric, const double Resid_fric,
                               const double d_wf, il::Array<double> &d);

il::Array<double> lin_friction(const double Peak_fric, const double Resid_fric,
                               const double d_wf, il::Array<double> &d);
}

#endif // HFPX2D_FRICTION_H
