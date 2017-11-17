//
// Created by Federico Ciardo on 13/11/17.
//

#ifndef HFPX2DUNITTEST_LOADCONDITIONS_H
#define HFPX2DUNITTEST_LOADCONDITIONS_H

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include <src/Core/Mesh.h>
#include <src/Input/FindUtilities.h>
#include <src/Core/InSituStress.h>

namespace hfp2d {

InSituStress loadConditions(Mesh &theLoadedMesh, const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &ConditionsMap);
}
#endif  // HFPX2DUNITTEST_LOADCONDITIONS_H
