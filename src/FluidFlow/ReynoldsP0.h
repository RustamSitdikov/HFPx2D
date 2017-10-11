//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#ifndef HFPX2DUNITTEST_REYNOLDSP0_H

#define HFPX2DUNITTEST_REYNOLDSP0_H

#include "il/Array.h"


#include <src/core/Fluid.h>
#include <src/core/Mesh.h>
#include <src/core/SolutionAtT.h>
#include <src/core/RockProperties.h>
#include <src/core_dev/Sources.h>

namespace hfp2d{

il::Array<double>  EdgeConductivitiesP0Newtonian(il::Array2D<il::int_t> &edgeAdj,il::Array<double> &width,hfp2d::Fluid &fluid);

il::Array2D<il::int_t> GetEdgesSharing2(hfp2d::Mesh &mesh);


il::Array2D<double> BuildFD_P0(hfp2d::Mesh &mesh, hfp2d::Fluid &fluid,
                               il::Array<double> &hydraulicwidth, double coef);


hfp2d::SolutionAtT ReynoldsSolverP0(hfp2d::SolutionAtT &soln,
                                    il::Array2D<double> &ElasMat,
                                    hfp2d::Fluid &fluid,
                                    hfp2d::RockProperties &rock,
                                    hfp2d::Sources &source, double timestep) ;

}

#endif //HFPX2DUNITTEST_REYNOLDSP0_H
