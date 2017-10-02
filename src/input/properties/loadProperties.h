//
// Created by lorenzo on 9/9/17.
//

#ifndef HFPX2DUNITTEST_LOADSOLID_H
#define HFPX2DUNITTEST_LOADSOLID_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include "src/core/Properties.h"
#include "loadSingleMaterial.h"
#include <src/input/findUtilities.h>
#include <src/core/Mesh.h>

namespace hfp2d {

Properties loadProperties(const Mesh &theLoadedMesh,
                          const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &solidMaterialMap);

}
#endif //HFPX2DUNITTEST_LOADSOLID_H
