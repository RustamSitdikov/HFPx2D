//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <gtest/gtest.h>

#include "src/Solvers/SimpleElastic.h"

////////////////////////////////////////////////////////////////////////////////

//   TEST griffith crack with 10 P1 elements
TEST(griffith_crack_P1, UnitP1) {
  double ret = hfp2d::SimpleGriffithExampleLinearElement(10);

  ASSERT_NEAR(0.0333916, ret, 0.001);

}


TEST(griffith_crack_P0, UnitP0) {
  double ret = hfp2d::SimpleGriffithExampleS3D_P0(10);

  ASSERT_NEAR(0.0045515, ret, 0.001);
}


