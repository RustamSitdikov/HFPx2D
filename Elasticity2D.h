//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_ELASTICITY2D_H
#define HFPX2D_ELASTICITY2D_H

// Inclusion from Inside Loop library
//#include <il/Array2D.h>
#include <il/StaticArray2D.h>
//#include <il/Array.h>
#include <il/StaticArray.h>

namespace hfp2d {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(const double h,
                                                       const double Ep,
                                                       const double x,
                                                       const double y);

// void normal_shear_stress_kernel_dp1_dd(il::StaticArray2D<double, 2, 4> &St,
//                                       const il::StaticArray<double, 2> xe,
//                                       const double &h,
//                                       const il::StaticArray<double, 2> s,
//                                       const il::StaticArray<double, 2> n,
//                                       const double Ep);
il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    const il::StaticArray<double, 2> xe, const double &h,
    const il::StaticArray<double, 2> s, const il::StaticArray<double, 2> n,
    const double Ep);
}

#endif // HFPX2D_ELASTICITY2D_H