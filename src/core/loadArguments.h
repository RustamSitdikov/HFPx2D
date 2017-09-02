//
// Created by lorenzo on 9/2/17.
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
}

#endif //HFPX2DUNITTEST_LOADARGUMENTS_H
