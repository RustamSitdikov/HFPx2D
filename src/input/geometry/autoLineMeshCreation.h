//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_CUSTOMORIENTATIONMESH_H
#define HFPX2DUNITTEST_CUSTOMORIENTATIONMESH_H

#include <iostream>
#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/toml.h>
#include <il/String.h>
#include <il/linear_algebra.h>
#include "src/core/Mesh.h"
#include "src/input/geometry/autoLineMeshUtilities.h"
#include "src/input/geometry/autoLineMeshInfo.h"
#include "src/input/findUtilities.h"


namespace hfp2d {

////////////// CUSTOM MESH //////////////

Mesh autoLineMesh(const il::String &inputFileName,
                  const il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                  const il::int_t interpOrder);

}

#endif //HFPX2DUNITTEST_CUSTOMORIENTATIONMESH_H
