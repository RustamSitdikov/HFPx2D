//
// Created by lorenzo on 9/4/17.
//

#ifndef HFPX2DUNITTEST_LOADINPUT_H
#define HFPX2DUNITTEST_LOADINPUT_H

#include <iostream>
#include "il/toml.h
#include "Mesh.h"
#include "SimulationParameters.h"

namespace hfp2d{

void loadInput(std::string inputFileName,il::io_t,
               Mesh mesh,
               Properties matProperties,
               Simulation simParameters);

}

#endif //HFPX2DUNITTEST_LOADINPUT_H
