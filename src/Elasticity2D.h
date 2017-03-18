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

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd( double h,
                                                        double Ep,
                                                        double x,
                                                        double y);

il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    const il::StaticArray<double, 2>& xe, double h,
    const il::StaticArray<double, 2>& s, const il::StaticArray<double, 2>& n,
    double Ep);

// Simplified 3D kernel - piece wise constant
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(double a, double b,
                                                          double G, double nu,
                                                          double xx, double yy);

il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_s3d_dp0_dd(
    const il::StaticArray<double, 2>& xe, double hx, double height,
    const il::StaticArray<double, 2>& s, const il::StaticArray<double, 2>& n,
    double G, double nu);
}

#endif  // HFPX2D_ELASTICITY2D_H