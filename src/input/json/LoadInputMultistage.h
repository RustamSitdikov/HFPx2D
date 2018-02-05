//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 07.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_LOADINPUTMULTISTAGE_H
#define HFPX2D_LOADINPUTMULTISTAGE_H


#include <src/core/SolidProperties.h>
#include <src/core/Fluid.h>
#include <src/core/InSituConditions.h>
#include <src/wellbore/WellMesh.h>
#include <src/wellbore/WellInjection.h>

#include <src/util/json.hpp>

namespace hfp2d {
using json = nlohmann::json;


hfp2d::WellMesh  loadWellMesh(json &j_wmesh);

hfp2d::WellInjection loadWellParameters(json &j_params,
                                        hfp2d::WellMesh &the_well);

hfp2d::Fluid loadFluidProperties(json &j_fluid);

hfp2d::Sources loadWellSource(json &j_params,hfp2d::WellMesh &the_well);

hfp2d::SolidProperties loadSolidProperties(json &j_rock);

hfp2d::InSituConditions loadInSitu(json &j_insitu);

};

#endif //HFPX2D_LOADINPUTMULTISTAGE_H
