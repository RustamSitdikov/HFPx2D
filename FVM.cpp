//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "FVM.h"
#include <iostream>

// Some utilities
il::Array<double> Average(il::Array2D<double> &d){

    il::Array<double> Average{d.size(0),0.};

    for (il::int_t i = 0; i < d.size(0); ++i) {

            Average[i] = (d(i,0) + d(i,1))/2;

    }

    return Average;

};

//
il::Array<double> Quarter(il::Array2D<double> &d){


  il::Array<double> Quarter(2*d.size(0),0.);

    for (il::int_t i = 0, j=0; i < (d.size(0)); ++i, j=j+2) {


        Quarter[j] = ((3*d(i,0))+d(i,1))/4;
        Quarter[j+1] = (d(i,0)+(3*d(i,1)))/4;

    }

    return Quarter;

};
