//
// HFPx2D_Shear project.
//
// Created by Federico Ciardo on 03.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include "src/Problems/FluidInj_FricDilatFault.h"

////////////////////////////////////////////////////////////////////////////////

int main() {
  std::cout << "\n\n ----- Problem: fluid injection into a frictional "
               "weakening dilatant fault ----- \n\n"
            << std::endl;

  int Nelts = 3;

  int ret = hfp2d::FluidInjectionFrictionalWeakeningFault(Nelts);
  std::cout << "return " << ret;

  std::cout << " ----- End of simulation ----- \n\n\n";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////
