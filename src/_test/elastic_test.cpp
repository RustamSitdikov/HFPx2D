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

TEST(griffith_crack_1, UnitP) {
  double ret = hfp2d::SimpleGriffithExample(10);

  ASSERT_NEAR(0.0333916, ret, 0.001);
}
