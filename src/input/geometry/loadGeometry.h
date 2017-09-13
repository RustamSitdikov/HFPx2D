//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_GEOMETRY_H
#define HFPX2DUNITTEST_GEOMETRY_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include "src/core/Mesh.h"
#include "verticalOrientation.h"
#include "horizontalOrientation.h"
#include "diagonalOrientation.h"
#include "customOrientation.h"
#include "loadMeshFile.h"


namespace hfp2d {

void loadGeometry(const il::String &inputFileName,
              const il::MapArray<il::String, il::Dynamic> &meshCreationMap,
              il::io_t,
              Mesh &theMesh);

}

#endif //HFPX2DUNITTEST_GEOMETRY_H
