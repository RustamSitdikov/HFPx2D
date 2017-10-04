//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2DUNITTEST_LOADCONDITIONS_H
#define HFPX2DUNITTEST_LOADCONDITIONS_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include "src/core/Properties.h"
#include "src/core/Conditions.h"
#include <src/input/findUtilities.h>
#include <src/core/Mesh.h>

namespace hfp2d {

Conditions loadConditions(const Mesh &theLoadedMesh,
                          const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &conditionsMap);

}
#endif //HFPX2DUNITTEST_LOADCONDITIONS_H
