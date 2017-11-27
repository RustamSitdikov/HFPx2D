//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 28.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <il/math.h>
#include <iostream>

// Inclusion from the project
#include "ConductivitiesNewtonian.h"

namespace hfp2d {

// Output: array (vector) that contains volume of fluid for each element
il::Array<double> conductivities_newtonian(const il::Array<double> &rho,
                                           const il::Array<double> &vector,
                                           il::Array<double> EltSizes,
                                           Parameters_fluid &fluid_parameters,
                                           const il::Array<double> &permeab,
                                           il::io_t) {

  // Inputs:
  //  - rho -> vector of fluid density values at the middle of each element
  //  (size -> Nelts)
  //  - vector -> vector of dilatancy value for slip at the middle of the
  //  elements (size -> Nelts)
  //  - EltSizes -> vector that contains the element sizes
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - permeab -> vector that contains the permeability for each element (size
  //  -> Nelts)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> Res{EltSizes.size(), 0};

  for (il::int_t i = 0; i < Res.size(); ++i) {

    Res[i] = ((rho[i] * (vector[i] * permeab[i])) / EltSizes[i]) *
             (1 / (12 * fluid_parameters.viscosity));
  }

  return Res;
}
}