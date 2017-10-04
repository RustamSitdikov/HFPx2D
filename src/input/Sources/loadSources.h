//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2DUNITTEST_LOADSOURCES_H
#define HFPX2DUNITTEST_LOADSOURCES_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include <src/input/findUtilities.h>
#include <src/core/Mesh.h>
#include <src/core_dev/Sources.h>

namespace hfp2d{

Sources loadSources(const Mesh &theLoadedMesh,
                    const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &sourcesMap);

il::int_t findSourceLocation(const double locationX,
                             const double locationY,
                             const Mesh &theLoadedMesh);

}

#endif //HFPX2DUNITTEST_LOADSOURCES_H
