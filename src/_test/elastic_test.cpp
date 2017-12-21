//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <gtest/gtest.h>

#include <src/solvers/SimpleElasticBenchmarks.h>

////////////////////////////////////////////////////////////////////////////////

//   TEST griffith crack with 10 P1 elements
TEST(UnitP1, griffith_10) {
  double ret = hfp2d::SimpleGriffithExampleLinearElement(10);

  ASSERT_NEAR(0.0333916, ret, 0.001); // norm on all DD excluding tip ends

}

TEST(UnitP1,griffith_10_add) {
  double ret = hfp2d::SimpleGriffithExampleLinearElement_AddMesh(10);

  ASSERT_NEAR(0.00625229, ret, 0.001);// norm on just the DD before the tip

}

//   TEST griffith crack with 10 P1 elements
TEST(UnitP1, griffith_10_by_nodes) {
  double ret = hfp2d::SimpleGriffithExampleLinearElement_byNodes(10);

  ASSERT_NEAR(0.0333916, ret, 0.001); // norm on all DD excluding tip ends

}


TEST( UnitP0,griffith_10) {
  double ret = hfp2d::SimpleGriffithExampleS3D_P0(10);

  ASSERT_NEAR(0.0045515, ret, 0.001);// norm on all DD excluding tip ends
}


TEST(UnitP0, griffith_10_add) {
  double ret = hfp2d::SimpleGriffithExampleS3D_P0_AddMesh(10);

  ASSERT_NEAR(0.0809615, ret, 0.001);   // norm on all elements including tip
}


TEST( UnitP0,griffith_10_by_nodes) {
  double ret = hfp2d::SimpleGriffithExampleS3D_P0_byNodes(10);

  ASSERT_NEAR(0.0045515, ret, 0.001);// norm on all DD excluding tip ends
}
