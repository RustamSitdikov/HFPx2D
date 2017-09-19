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

namespace hfp2d {

Properties loadProperties(const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &solidMaterialMap,
                          const il::int_t numOfMats);
/*
double findYoungModulus(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);
double findPoissonRatio(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);

double findFluidDensity(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);
double findFluidViscosity(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);
double findFluidCompressibility(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);

il::int_t findNumMaterials(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName);
*/
}
#endif //HFPX2DUNITTEST_LOADSOLID_H
