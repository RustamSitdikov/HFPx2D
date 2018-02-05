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

// Inclusion from standard library
#include <iostream>

// Inclusion from Inside Loop library
#include "il/toml.h"

// Inclusion from the project
#include <src/core/ElasticProperties.h>
#include <src/core/Mesh.h>
#include <src/core/SimulationParameters.h>
#include <src/core/Sources.h>
#include <src/core_dev/FractureEvolution.h>
#include <src/input/conditions/LoadConditions.h>
#include <src/input/geometry/LoadGeometry.h>
#include <src/input/properties/LoadProperties.h>

namespace hfp2d {

void loadInput(const il::String &config_filename, il::io_t, Mesh &MyMesh,
               ElasticProperties &ElasticProperties,
               FluidProperties &FluidProperties, SolidEvolution &SolidEvolution,
               FractureEvolution &FractureEvolution,
               InSituStress &BackgroundLoadingConditions,
               double &const_overpress, double &t_0plus1, double &time_step,
               double &time_step_max, double &final_time, bool &expl_impl,
               bool &quasi_dynamic, double &dilat_plast);

}

#endif  // HFPX2DUNITTEST_LOADINPUT_H
