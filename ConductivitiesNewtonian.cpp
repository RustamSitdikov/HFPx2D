//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 28.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>
#include <il/math.h>

// Inclusion from the project
#include "ConductivitiesNewtonian.h"

namespace hfp2d {

// Output: array (vector) that contains volume of fluid for each element
il::Array<double> conductivities_newtonian(const il::Array<double> &rho,
                                           const il::Array<double> &vector,
                                           il::Array<double> EltSizes,
                                           Parameters_fluid &fluid_parameters,
                                           il::io_t) {

  // Inputs:
  //  - rho -> vector of fluid density values at the middle of each element
  //  (size -> Nelts)
  //  - vector -> vector of shear DD or opening DD at the middle of the elements
  //  (size -> Nelts)
  //  - EltSizes -> vector that contains the element sizes
  //  - Visc -> fluid viscosity (floating point value)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> Res{EltSizes.size(), 0};

  for (il::int_t i = 0; i < Res.size(); ++i) {

//    Res[i] = ((rho[i] * (pow(vector[i], 3))) / EltSizes[i]) *
//             (1 / (12 * fluid_parameters.viscosity));
    Res[i] = ((rho[i] * (vector[i]*vector[i]*vector[i])) / EltSizes[i]) *
             (1 / (12 * fluid_parameters.viscosity));
  }

  return Res;
}
}