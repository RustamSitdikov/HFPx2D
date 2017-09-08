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
#include <cstring>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <il/String.h>

namespace hfp2d {

void loadArguments(const int argc, const char* const* argv, il::io_t,
                   bool &checkInput, il::String &inputFileName,
                   bool &checkRestart, il::String &restartFileName,
                   bool &checkOutput, il::String &outputDirectory);

void cleanOutputDir(const char *path);

}
#endif //HFPX2DUNITTEST_LOADARGUMENTS_H
