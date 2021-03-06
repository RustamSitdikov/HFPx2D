//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_PLANESTRAININFINITE_H
#define HFPX2D_PLANESTRAININFINITE_H

#include <il/StaticArray2D.h>
//#include <il/Array.h>
#include <il/StaticArray.h>

#include <src/core/ElasticProperties.h>
#include <src/core/Mesh.h>
#include <src/core/SegmentData.h>

namespace hfp2d {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(double h, double Ep,
                                                       double x, double y);


il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    SegmentData source_elt, SegmentData receiver_elt,  il::int_t i_col,
    ElasticProperties Elas, double ker_options);



// by nodal effect
il::StaticArray<double, 4> stresses_kernel_dp1_dd_nodal(il::int_t local_node_i,
                                                        double h, double Ep,
                                                        double x, double y);

il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_dp1_dd_nodal(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t s_col,
    il::int_t i_col, ElasticProperties Elas, double ker_options);

}

#endif  // HFPX2D_PLANESTRAININFINITE_H
