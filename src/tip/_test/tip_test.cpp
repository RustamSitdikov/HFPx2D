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
TEST(tip_inversion_1, t1) {
//
bool mute = false;
tip::TipParameters tipPar;

tipPar.k1c = 1.4E6;
tipPar.e_p = 3.2E10;
tipPar.cl = 3.0E-3;
tipPar.mu = 1.0e-3;

tipPar.s0 = 0.844598;
tipPar.vt = 0.0;

  tipPar.wa = 0.00026694;
  tipPar.dt = 0.312799;

    tip::tipInversion(tip::res_u_0_m, tipPar, 1000, 1E-6, 50, mute);
// the new tip distance is now tipPar.st and the new velocity is tipPar.vt
double rm = tip::res_u_0_m(tipPar.st, tipPar);

ASSERT_NEAR(0.84508554222845511, tipPar.st, 0.001);

ASSERT_NEAR(0.0, rm, 0.001);

}
