//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 24.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>

#include <gtest/gtest.h>
#include <il/Array2D.h>
#include <il/Array.h>

#include <src/core/Mesh.h>


TEST(Mesh, test_ribbon){
  //  a very simple wellMesh with 4 elements  (0,1,2,3)
  // ensure the ribbon elements are 1 and 2

  // create the wellMesh.
   il::int_t nelts = 4;

   il::Array2D<double> xy{nelts+1,2,0.};

//  // create a basic 1D wellMesh .... first fracture
   for (il::int_t i = 0; i < nelts + 1; ++i) {
      xy(i, 0) = -1.0 + 1.*i  ;
      xy(i, 1) = 0.;
   };

  il::Array2D<il::int_t> myconn{nelts,2,0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, 0);
//
  auto ribbon = mesh.getRibbonElements();
//
//  std::cout << "ribbon size: " << ribbon.size() << " "<<ribbon[0] << " " << ribbon[1];
//

  ASSERT_TRUE((ribbon[0]==1)&& (ribbon[1]==2)  );//&& (ribbon[1]==2)

}