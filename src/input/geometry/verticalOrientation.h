//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_VERTICALORIENTATIONMESH_H
#define HFPX2DUNITTEST_VERTICALORIENTATIONMESH_H

#include <iostream>
#include "il/toml.h"
#include "src/core/Mesh.h"
#include "src/input/geometry/createMesh.h"
#include "src/input/geometry/autoMeshUtilities.h"

namespace hfp2d {

////////////// VERTICAL MESH //////////////

Mesh verticalOrientationMesh(const il::String &inputFileName,
                                   il::int_t fractureID,
                             const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

}

#endif //HFPX2DUNITTEST_VERTICALORIENTATIONMESH_H
