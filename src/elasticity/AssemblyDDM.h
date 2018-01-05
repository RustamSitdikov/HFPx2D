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
#include <src/core/SegmentData.h>

namespace hfp2d {

typedef  il::StaticArray2D<double, 2, 4> (*vKernelCall)(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t i_col,
    ElasticProperties Elas, double ker_options);


il::Array2D<double> basic_assembly(Mesh &mesh, hfp2d::ElasticProperties &elas,
                                   vKernelCall KernelCall, double ker_options);

void basic_assembly_add_elts(Mesh &new_mesh, il::int_t n_add,
                             hfp2d::ElasticProperties &elas, vKernelCall KernelCall,
                             double ker_options, il::io_t,
                             il::Array2D<double> &K);


void AddTipCorrectionP0(hfp2d::Mesh &mesh, const hfp2d::ElasticProperties &elas,
                        il::int_t tipElt, il::Array2D<double> &Kmat );

void RemoveTipCorrectionP0(hfp2d::Mesh &mesh, const hfp2d::ElasticProperties &elas,
                           il::int_t tipElt, il::Array2D<double> &Kmat );

il::Array2D<double> ReArrangeKP0(const Mesh &mesh,il::Array2D<double> &Kmat);



typedef  il::StaticArray2D<double, 2,2> (*vKernelCallNode)(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t i_col,
    il::int_t j_col,  ElasticProperties Elas, double ker_options);


il::Array2D<double> basic_assembly_nodal(Mesh &mesh,
                                         hfp2d::ElasticProperties &elas,
                                         vKernelCallNode KernelCall,
                                         double ker_options);


}

#endif  // HFPX2D_ASSEMBLYDDM_H
