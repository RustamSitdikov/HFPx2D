//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 07.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_LOADINPUTMULTISTAGE_H
#define HFPX2D_LOADINPUTMULTISTAGE_H

#include <src/wellbore/WellMesh.h>
#include <src/wellbore/WellInjection.h>
#include <src/util/json.hpp>

namespace hfp2d {
using json = nlohmann::json;


hfp2d::WellMesh  LoadWellMesh(json &j_wmesh);

hfp2d::WellInjection LoadWellParameters(json &j_params,hfp2d::WellMesh &the_well);

};

#endif //HFPX2D_LOADINPUTMULTISTAGE_H
