//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 28.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_CONDUCTIVITIESNEWTONIAN_H
#define HFPX2D_CONDUCTIVITIESNEWTONIAN_H

// Inclusion from standard library
#include <cmath>

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

namespace hfp2d {

    il::Array<double> ConductivitiesNewtonian(const il::Array<double> &rho,
                                              const il::Array<double> &vector,
                                              const il::Array<double> EltSizes,
                                              const double Visc);
}

#endif // HFPX2D_CONDUCTIVITIESNEWTONIAN_H
