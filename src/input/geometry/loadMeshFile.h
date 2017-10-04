//
// Created by lorenzo on 9/8/17.
//

#ifndef HFPX2D_LOADMESHFILE_H
#define HFPX2D_LOADMESHFILE_H

#include <iostream>
#include "il/base.h"
#include "il/toml.h"
#include "il/Array.h"
#include "il/String.h"
#include "src/core/Mesh.h"

namespace hfp2d {

void loadMeshFile(const il::String &meshFileName, il::io_t, Mesh &theMesh);

}

#endif //HFPX2D_LOADMESHFILE_H
