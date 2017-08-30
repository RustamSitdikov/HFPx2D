//
// HFPx2D project.
//
// Created by Brice Lecampion on 06.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <iostream>

#include "src/Solvers/SimpleElastic.h"


////////////////////////////////////////////////////////////////////////////////
int main() {

  int nelts = 10;

  double ret = hfp2d::SimpleGriffithExample(nelts);

  std::cout << "\n rel error L2 norm: " << ret << "\n";

  std::cout << " end of code";

  return ret;
}
////////////////////////////////////////////////////////////////////////////////

