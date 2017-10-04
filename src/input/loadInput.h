//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 04.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_LOADINPUT_H
#define HFPX2DUNITTEST_LOADINPUT_H

#include <iostream>
#include "il/toml.h"
#include "src/core/Mesh.h"
#include "src/core/Properties.h"
#include "src/core/Simulation.h"
#include "src/core/Conditions.h"
#include "src/core/Sources.h"
#include "src/input/geometry/loadGeometry.h"
#include "src/input/properties/loadProperties.h"
#include "src/input/conditions/loadConditions.h"
#include "src/input/Sources/loadSources.h"

namespace hfp2d{

void loadInput(const il::String &inputFileName,
               il::io_t,
               Mesh &theMesh,
               Properties &theProperties,
               Conditions &theConditions,
               Sources &theSources,
               Simulation &simParameters);

}

#endif //HFPX2DUNITTEST_LOADINPUT_H
