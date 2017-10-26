//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#include <gtest/gtest.h>

#include <src/tip/tipAsymptote.h>

//  UNIT TESTS FOR TIP Asymptote Libary
TEST(tip_inversion_1,t1) {
  //
   double ret = 0.0333916;

  ASSERT_NEAR(0.0333916, ret, 0.001);

}
