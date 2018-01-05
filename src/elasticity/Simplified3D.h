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

#include <src/core/Mesh.h>
#include <src/core/ElasticProperties.h>
#include <src/core/SegmentData.h>

namespace hfp2d {


// Simplified 3D kernel - piece wise constant
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(double a, double b,
                                                          double G, double nu,
                                                          double xx, double yy);

il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_s3d_dp0_dd(
    SegmentData source_elt, SegmentData receiver_elt,  il::int_t i_col,
    ElasticProperties Elas, double ker_options);


il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_s3d_dp0_dd_nodal(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t s_col,
    il::int_t i_col, ElasticProperties Elas, double ker_options);

}

#endif // HFPX2D_ELASTICITY2D_H