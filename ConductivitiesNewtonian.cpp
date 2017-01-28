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

// Inclusion from the project
#include "ConductivitiesNewtonian.h"

// Output: array (vector) that contains volume of fluid for each element
il::Array<double> ConductivitiesNewtonian(const il::Array<double> &rho,
                                          const il::Array<double> &vector,
                                          const il::Array<double> EltSizes,
                                          const double Visc) {

  // Inputs:
  //  - rho -> vector of fluid density values (size -> Nelts)
  //  - vector -> vector of shear DD or opening DD at the middle of the elements
  //  (size -> Nelts)
  //  - EltSizes -> vector that contains the element sizes
  //  - Visc -> fluid viscosity (floating point value)

  il::Array<double> Res{EltSizes.size(), 0.};

  for (il::int_t i = 0; i < Res.size(); ++i) {

    Res[i] = ((rho[i] * (pow(vector[i], 3))) / EltSizes[i]) * (1 / (12 * Visc));
  }

  return Res;
}