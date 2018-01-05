//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//


#ifndef HFPX2DUNITTEST_DEVLHFMSOLVER_H
#define HFPX2DUNITTEST_DEVLHFMSOLVER_H

#include <src/core/SimulationParameters.h>
#include <src/core/Solution.h>
#include <src/core/Fluid.h>
#include <src/core/SolidProperties.h>
#include <src/core/Sources.h>

namespace hfp2d {

int TwoParallelHFs(int nelts, double dist);


hfp2d::Solution FractureFrontLoop(
    const  hfp2d::Solution &Sol_n, il::Array2D<double> &ElasMat, hfp2d::Fluid &fluid,
    const hfp2d::SolidProperties &rock,const hfp2d::Sources &source, double frac_height,
    double timestep, hfp2d::SimulationParameters &simulParams, bool mute);

}


#endif //HFPX2DUNITTEST_DEVLHFMSOLVER_H
