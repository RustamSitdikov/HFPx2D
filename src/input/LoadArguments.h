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

#include <dirent.h>
#include <il/String.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iostream>

namespace hfp2d {

void loadArguments(int argc, char const *argv[], il::io_t,
                   il::String &input_filename,
                   il::String &path_output_directory);

void cleanOutputDir(const char *path);
}
#endif  // HFPX2DUNITTEST_LOADARGUMENTS_H
