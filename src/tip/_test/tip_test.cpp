//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#include <gtest/gtest.h>

#include <src/tip/tipAsymptote.h>
#include "../tipAsymptote.h"

//  UNIT TESTS FOR TIP Asymptote Libary
TEST(tip_inversion_1, t1) {
//
bool mute = false;
tip::TipParameters prevTipPar;
prevTipPar.k1c = 1.4E6;
prevTipPar.e_p = 3.2E10;
prevTipPar.cl = 3.0E-3;
prevTipPar.mu = 1.0e-3;
prevTipPar.wa = 0.0002;
prevTipPar.st = 0.844598;
prevTipPar.dt = 0.3127;
prevTipPar.vt = 0.0;

double dt = 0.312799;
double wa = 0.00026694;

tip::TipParameters newTipPar = tip::propagateTip
        (tip::res_u_0_m, prevTipPar,
         dt, wa, 1000,
         1E-6, 50, mute);

double rm = tip::res_u_0_m(newTipPar.st, newTipPar);

ASSERT_NEAR(0.84508554222845511, newTipPar.st, 0.001);

ASSERT_NEAR(0.0, rm, 0.001);

}
