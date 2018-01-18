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
#include <il/Array.h>
#include <il/Array2D.h>

#include <src/core/Mesh.h>

/// TEST 1
TEST(Mesh, test_ribbon) {
  //  a very simple wellMesh with 4 elements  (0,1,2,3)
  // ensure the ribbon elements are 1 and 2

  // create the wellMesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  //  // create a basic 1D wellMesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, 0);
  //
  auto ribbon = mesh.getRibbonElements();
  //
  //  std::cout << "ribbon size: " << ribbon.size() << " "<<ribbon[0] << " " <<
  //  ribbon[1];
  //

  ASSERT_TRUE((ribbon[0] == 1) && (ribbon[1] == 2));  //&& (ribbon[1]==2)
}

/// TEST 2
TEST(Mesh, test_mesh_object) {
  //  a very simple Mesh with 4 elements  (0,1,2,3)
  // ensure the mesh class works correctly

  // create the Mesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  // create a basic 1D Mesh
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  il::int_t e = 0;

  hfp2d::Mesh mesh(xy, myconn, 0);

  ASSERT_TRUE((mesh.numberOfElts() == nelts) &&
              (mesh.numberOfElts() == myconn.size(0)) &&
              (mesh.numberOfNodes() == xy.size(0)) &&
              (mesh.interpolationOrder() == 0) && (mesh.tipElts(0) == 0) &&
              (mesh.tipElts(1) == 3) && (mesh.tipNodes(0) == 0) &&
              (mesh.tipNodes(1) == 4) && (mesh.eltSize(e) == 1.) &&
              (mesh.allEltSize()[0] == 1.));
}