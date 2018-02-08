//
// Created by Federico Ciardo on 13/11/17.
//

#ifndef HFPX2DUNITTEST_LOADCONDITIONS_H
#define HFPX2DUNITTEST_LOADCONDITIONS_H

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include <src/core/Mesh.h>
#include <src/input/toml/findUtilities.h>
#include <src/core/InSituStress.h>

namespace hfp2d {

InSituStress loadConditions(Mesh &theLoadedMesh, const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &ConditionsMap);

}
#endif  // HFPX2DUNITTEST_LOADCONDITIONS_H
