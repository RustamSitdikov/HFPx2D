//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_GEOMETRY_H
#define HFPX2D_GEOMETRY_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include "src/core/Mesh.h"
#include "autoLineMeshCreation.h"
#include "loadMeshFile.h"


namespace hfp2d {

Mesh loadGeometry(const il::String &inputFileName, const il::MapArray<il::String, il::Dynamic> &meshCreationMap);

il::int_t findInterpOrder(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);
il::int_t findNumFractures(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);

}

#endif //HFPX2D_GEOMETRY_H
