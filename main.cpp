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

int main(int argc, char const *argv[]){

  hfp2d::fluidInjFrictWeakDilatFault(argc - 1, argv + 1);

  std::cout << " End of code \n\n\n";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

