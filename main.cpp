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
#include "src/core/SolutionClass.h"


////////////////////////////////////////////////////////////////////////////////
int main() {

  int nelts = 10;

  double ret = hfp2d::SimpleGriffithExample(nelts);

  std::cout << "\n rel error L2 norm: " << ret << "\n";

  std::cout << " end of code \n\n\n";




  hfp2d::Solution a(5,1);
  hfp2d::Solution a1(7,1,3.0);
  hfp2d::Solution a2(a);

    il::int_t numElem=10;

  il::Array<double> b(numElem,0.0);




  //std::cout << a.displacement().size() << std::endl;
  for(il::int_t i=0; i<10; i++)
  {
    std::cout << a1.displacement(i) << std::endl;
  }

//  for (int i=0; i<5; i++) {
//    a.displacement(i,2.0*i);
//  }

  return ret;
}
////////////////////////////////////////////////////////////////////////////////

