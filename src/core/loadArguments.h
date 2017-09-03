//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 02.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_LOADARGUMENTS_H
#define HFPX2DUNITTEST_LOADARGUMENTS_H

#include <iostream>
#include <fstream>
#include <sys/stat.h>

namespace hfp2d {

void loadArguments(int argc, char *argv[],
                   bool &checkInput, std::string &inputFileName,
                   bool &checkRestart, std::string &restartFileName,
                   bool &checkOutput, std::string &outputDirectory);

void cleanOutputDir(const char *path);
}
#endif //HFPX2DUNITTEST_LOADARGUMENTS_H
