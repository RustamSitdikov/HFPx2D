//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2DUNITTEST_REYNOLDSP0_H

#define HFPX2DUNITTEST_REYNOLDSP0_H

#include "il/Array.h"

#include <src/core/Fluid.h>
#include <src/core/Mesh.h>
#include <src/core/SimulationParameters.h>
#include <src/core/SolidProperties.h>
#include <src/core/Solution.h>
#include <src/core/Sources.h>

namespace hfp2d {

il::Array<double> EdgeConductivitiesP0Newtonian(il::Array2D<il::int_t> &edgeAdj,
                                                il::Array<double> &width,
                                                hfp2d::Fluid &fluid);

il::Array2D<double> BuildFD_P0(hfp2d::Mesh &mesh, hfp2d::Fluid &fluid,
                               il::Array<double> &hydraulicwidth, double coef);

Solution ReynoldsSolverP0(hfp2d::Solution &soln,
                          il::Array2D<double> &ElasMat,
                           hfp2d::Fluid &fluid,
                          const hfp2d::SolidProperties &rock,
                          const hfp2d::Sources &source, double timestep,
                          bool imp_tip_width,
                          il::Array<il::int_t> &tip_region_elt,
                          il::Array<double> &tip_width,
                          hfp2d::SimulationParameters &simulParams, bool mute);
}

#endif  // HFPX2DUNITTEST_REYNOLDSP0_H
