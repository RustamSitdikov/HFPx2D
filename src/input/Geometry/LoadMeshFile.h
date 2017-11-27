//
// Created by lorenzo on 9/8/17.
//

#ifndef HFPX2DUNITTEST_LOADMESHFILE_H
#define HFPX2DUNITTEST_LOADMESHFILE_H

// Inclusion from standard library
#include <iostream>

// Inclusion from Inside Loop library
#include "il/Array.h"
#include "il/String.h"

// Inclusion from the project
#include "src/core/Mesh.h"

namespace hfp2d {

void loadMeshFile(const il::String &meshFileName, il::io_t, Mesh &theMesh);

}

#endif //HFPX2DUNITTEST_LOADMESHFILE_H
