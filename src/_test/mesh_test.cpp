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
#include <il/base.h>

#include <src/core/Mesh.h>


hfp2d::Mesh simplemesh(il::int_t nelts){

  il::Array2D<double> xy(nelts, 0);
  il::Array2D<il::int_t> myconn(nelts,0);

  // create a basic 1D mesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + i  ;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, 0);

  return mesh;
}

TEST(Mesh, test_ribbon){
  //  a very simple mesh with 4 elements  (0,1,2,3)
  // ensure the ribbon elements are 1 and 2

  // create the mesh.

  il::int_t nelts = 4;

  hfp2d::Mesh mesh =simplemesh(nelts);

  il::Array<il::int_t> ribbon = mesh.getRibbonElements();

  std::cout << "ribbon size: " << ribbon.size() << " "<<ribbon[0] << " " << ribbon[1];

 ASSERT_TRUE((ribbon[0]==1)  );//&& (ribbon[1]==2)

}