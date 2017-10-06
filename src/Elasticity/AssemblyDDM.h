//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_ASSEMBLYDDM_H
#define HFPX2D_ASSEMBLYDDM_H

// Inclusion from Inside Loop library
#include <il/Array2D.h>

// Inclusion from the project
#include <src/core/Mesh.h>
#include <src/core/ElasticProperties.h>
#include <src/core/segmentData.h>

namespace hfp2d {

typedef  il::StaticArray2D<double, 2, 4> (*vKernelCall)(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t i_col,
    ElasticProperties Elas, double ker_options);


il::Array2D<double> basic_assembly(Mesh &mesh, ElasticProperties &elas,
                                   vKernelCall KernelCall, double ker_options);


void take_submatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
                    const il::Array2D<double> &A);

void set_submatrix(il::Array2D<double> &A, int i0, int i1,
                   const il::StaticArray2D<double, 2, 4> &B);

void AddTipCorrectionP0(const Mesh &mesh, const ElasticProperties &elas,
                        il::int_t tipElt, il::Array2D<double> &Kmat );

void RemoveTipCorrectionP0(const Mesh &mesh, const ElasticProperties &elas,
                           il::int_t tipElt, il::Array2D<double> &Kmat );

il::Array2D<double> ReArrangeKP0(const Mesh &mesh,il::Array2D<double> &Kmat);

}

#endif  // HFPX2D_ASSEMBLYDDM_H
