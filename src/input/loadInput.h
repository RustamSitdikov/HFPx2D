//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 04.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_LOADINPUT_H
#define HFPX2D_LOADINPUT_H

// Std libs
#include <iostream>

// IL libs
#include <il/toml.h>

// HFP2D libs
#include "src/core/Mesh.h"
#include "src/core_dev/Properties.h"
#include "src/core_dev/Simulation.h"
//#include "src/core/InSituStress.h"
#include "src/core_dev/Sources.h"
#include "src/input/geometry/loadGeometry.h"
#include "src/input/properties/loadProperties.h"
//#include "src/input/conditions/loadConditions.h"
#include "src/input/Sources/loadSources.h"

namespace hfp2d{

void loadInput(const il::String &inputFileName,
               il::io_t,
               Mesh &theMesh,
               Properties &theProperties,
               //InSituStress &theConditions,
               Sources &theSources,
               simulationParams &simParameters);

}

#endif //HFPX2D_LOADINPUT_H
