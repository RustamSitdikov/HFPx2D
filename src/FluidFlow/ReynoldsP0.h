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

namespace hfp2d{

il::Array<double>  EdgeConductivitiesP0Newtonian(il::Array2D<il::int_t> &edgeAdj,il::Array<double> &width,hfp2d::Fluid &fluid);




}

#endif //HFPX2DUNITTEST_REYNOLDSP0_H
