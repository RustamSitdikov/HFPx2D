//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "Dilatancy.h"
#include <cmath>

// Function that return an array that contain dilatancy values (according to exponential dilatant hardening law)
il::Array<double> Dilatancy(const double Init_dil, const double Incr_dil, const double d_wd, il::Array<double> &d){

    il::Array<double> D{d.size(),0.};

    for (il::int_t i = 0; i < D.size(); ++i) {

       D[i] = Init_dil + (Incr_dil*(1-exp(-d[i]/d_wd)));
    }

    return D;

};

// Function that return an array that contain the derivative w.r.t slip of dilatancy values (according to exponential dilatant hardening law)
il::Array<double> DDilatancy(const double Incr_dil, const double d_wd, il::Array<double> &d){

    il::Array<double> DD{d.size(),0.};

    for (il::int_t i = 0; i < DD.size(); ++i) {

      DD[i] = (Incr_dil/d_wd)*(exp(-d[i]/d_wd));

    }

    return DD;

};
