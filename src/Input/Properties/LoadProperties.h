//
// Created by lorenzo on 9/9/17.
//

#ifndef HFPX2DUNITTEST_LOADSOLID_H
#define HFPX2DUNITTEST_LOADSOLID_H

// Inclusion from standard library
#include <iostream>

// Inclusion from Inside Loop library
#include "il/base.h"
#include "il/toml.h"

// Inclusion from the project
#include <src/Core/FluidProperties.h>
#include <src/Core/Mesh.h>
#include <src/Core_dev/FractureEvolution.h>
#include <src/Core_dev/SolidEvolution.h>
#include <src/Input/FindUtilities.h>
#include "LoadSingleMaterial.h"
#include "src/Core/ElasticProperties.h"

namespace hfp2d {

void loadProperties(const Mesh &theLoadedMesh, const il::String &input_filename,
                    const il::MapArray<il::String, il::Dynamic> &propertiesMap,
                    ElasticProperties &ElasticProperties,
                    FluidProperties &FluidProperties,
                    SolidEvolution &SolidEvolution,
                    FractureEvolution &FractureEvolution);
}
#endif  // HFPX2DUNITTEST_LOADSOLID_H
