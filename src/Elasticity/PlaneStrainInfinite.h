//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_PLANESTRAININFINITE_H
#define HFPX2D_PLANESTRAININFINITE_H

#include <il/StaticArray2D.h>
//#include <il/Array.h>
#include <il/StaticArray.h>

namespace hfp2d {


il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(  double h,
                                                         double Ep,
                                                         double x,
                                                         double y);


il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    const il::StaticArray<double, 2> &xe, double h,
    const il::StaticArray<double, 2> &s, const il::StaticArray<double, 2> &n,
    double Ep);

}

#endif //HFPX2D_PLANESTRAININFINITE_H
