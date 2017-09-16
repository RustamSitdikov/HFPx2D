//
// HFPx2D project.
//
// Created by Brice Lecampion on 06.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <iostream>
#include <fstream>

#include "src/input/loadInput.h"
#include "src/Solvers/SimpleElastic.h"
#include "src/core/SolutionClass.h"
#include "src/core/loadArguments.h"


////////////////////////////////////////////////////////////////////////////////

/// Proposed Structure of the code
//  Introduce input data in the code
//  .. Load arguments passed to the main and deal with the different options
//  .... -i --> input file
//  .... -o --> output directory
//  .... -r --> restart file
//  .. Load variables in the input file (mesh, connectivity, material) or
//     provide default values/errors if something is missing
//  .. Create DOF handles
//  Setup problem (if needed prepare all constant matrices)
//  Start timestep loop
//  .. Compute system matrix
//  .. Compute force vector
//  .. Solve for the new displacements and pressures
//  .. Update variables (dependent variable, historical variables, etc)
//  Reject/Accept solution
//  .. Compute the required criteria
//  .. Adjust time step
//  Output and save the results (on not!)


////////////////////////////////////////////////////////////////////////////////
//int const argc, char const* const* argv
int main() {
/*
  // Creating variables to deal with input arguments
  il::String inputFileName, outputDirectory, restartFileName;
  bool checkInput = false;
  bool checkOutput = false;
  bool checkRestart = false;

  hfp2d::loadArguments(argc, argv, il::io,
                       checkInput, inputFileName,
                       checkRestart, restartFileName,
                       checkOutput, outputDirectory);

  hfp2d::Mesh theMesh;
  hfp2d::Properties matProperties;
  hfp2d::Simulation simParameters;

  /// TAKING CARE OF THE INPUT IN CASE OF NEW ANALYSIS OR RESTART
  if (checkInput) {

    hfp2d::loadInput(inputFileName, il::io, theMesh, matProperties, simParameters);

  } else if (checkRestart) {
    // hfp2d::loadRestart(resetFileName,il::io_t,Mesh,Properties,SimulationParam,Solution);
  }

  if (checkOutput) { // Eliminate this once a script for the output is done

    // Example output to DUMMY file
    il::String outputFile = il::join(outputDirectory, "/", "cracklength.txt");
    std::cout << outputDirectory << std::endl;
    std::cout << outputFile << std::endl;

    std::ofstream foutlc;
    foutlc.open(outputFile.asCString());
    for (int i = 0; i < 10; i++) {
      foutlc << i << " again good" << "\n";
    }
    foutlc << "Good bye, once again" << "\n";
    foutlc.close();

  }*/

//////////////////////// Prepare data for computation /////////////////////////////////

//////////////////////// Initiate the computational loop /////////////////////////////////

//////////////////////// Previous code snippet /////////////////////////////////

  int nelts = 10;

  double ret = hfp2d::SimpleGriffithExampleLinearElement(nelts);

  std::cout << "\n rel error L2 norm: " << ret << "\n";

  std::cout << " end of code \n\n\n";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

