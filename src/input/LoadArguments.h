//
// HFPx2D project.
//
// Created by Federico Ciardo on 02.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_LOADARGUMENTS_H
#define HFPX2DUNITTEST_LOADARGUMENTS_H

// Inclusion from standard library
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iostream>

// Inclusion from Inside Loop library
#include <il/String.h>

namespace hfp2d {

void loadArguments(int argc, char const *argv[], il::io_t,
                   il::String &config_filename,
                   il::String &path_output_directory);

void cleanOutputDir(const char *path);

}
#endif  // HFPX2DUNITTEST_LOADARGUMENTS_H
