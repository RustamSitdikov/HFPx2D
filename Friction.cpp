//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "Friction.h"
#include <cmath>


// Function that return an array that contain friction coefficient (according to EXPONENTIAL friction weakening law)
il::Array<double> Exp_Friction(const double Peak_fric, const double Resid_fric, const double d_wf, il::Array<double> &d){

    il::Array<double> f{d.size(),0.};

    for (il::int_t i = 0; i < f.size(); ++i) {

        f[i] = Peak_fric - ((Peak_fric - Resid_fric)*(1 - exp(-(d[i]/d_wf))));
    }

    return f;

};

// Function the return an array that contain friction coefficient (according to LINEAR friction weakening law)
il::Array<double> Lin_Friction(const double Peak_fric, const double Resid_fric, const double d_wf, il::Array<double> &d){

    il::Array<double> f{d.size(),0.};
    double_t sl;
    sl = Peak_fric/d_wf;

    for (il::int_t i = 0; i < f.size(); ++i) {

        if (d[i] < d_wf){

            f[i] = Peak_fric - (sl*d[i]);
        }

        else {

        f[i] = Resid_fric;

        }
//
    }

    return f;

};
