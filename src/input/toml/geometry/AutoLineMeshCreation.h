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

// Inclusion from standard library
#include <iostream>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>
#include <il/linear_algebra.h>

// Inclusion from the projec
#include "src/core/Mesh.h"
#include "AutoLineMeshUtilities.h"
#include "AutoLineMeshInfo.h"
#include "src/input/toml/findUtilities.h"

namespace hfp2d {

////////////// CUSTOM MESH //////////////

Mesh autoLineMesh(const il::String &inputFileName,
                  const il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                  const il::int_t interpOrder);

}

#endif //HFPX2DUNITTEST_CUSTOMORIENTATIONMESH_H
