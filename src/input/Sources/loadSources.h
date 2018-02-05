//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 05.02.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_LOADSOURCES_H
#define HFPX2D_LOADSOURCES_H

#include <src/core/Mesh.h>
#include <src/core/Sources.h>
#include <src/input/findUtilities.h>
#include <iostream>
#include "il/base.h"
#include "il/toml.h"

namespace hfp2d {

Sources loadSources(const Mesh &theLoadedMesh, const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &sourcesMap);

il::int_t findSourceLocation(const double locationX, const double locationY,
                             const Mesh &theLoadedMesh);
}

#endif  // HFPX2D_LOADSOURCES_H
