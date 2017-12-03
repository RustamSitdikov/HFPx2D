//
// HFPx2D project.
//
// Created by Federico Ciardo on 08.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include "src/solvers/FluidInjFrictWeakDilatFault.h"

//////////////////////// Numerical code /////////////////////////////////

int main(int argc, char const *argv[]) {

  std::cout << "\n";
  std::cout << " ->    Begin of simulation    <-";
    std::cout << "\n\n";

  // Call the solver of the problem of fluid injection into a frictional
  // weakening dilatant fault
  hfp2d::fluidInjFrictWeakDilatFault(argc - 1, argv + 1);

  std::cout << " ->    End of simulation    <-";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
