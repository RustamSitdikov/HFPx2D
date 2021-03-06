//
// HFPx2D project.
//
// Created by Brice Lecampion on 06.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cmath>
#include <fstream>
#include <iostream>

//#include "src/input/loadInput.h"
#include "src/solvers/SimpleElasticBenchmarks.h"
//#include "src/core_dev/OldSolutionClass.h"
//#include "src/input/loadArguments.h"


#include <src/solvers/HFPropagationP0.h>
#include <src/wellbore/SimpleWellFlowSolverBenchmark.h>
#include <src/solvers/MultiFracsSolver.h>

////////////////////////////////////////////////////////////////////////////////

/// Proposed Structure of the code
//  Introduce input data in the code
//  .. Load arguments passed to the main and deal with the different options
//  .... -i --> input file
//  .... -o --> output directory
//  .... -r --> restart file
//  .. Load variables in the input file (wellMesh, connectivity, material) or
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
//
//
// int main(int const argc, char const* const* argv) {
//
//  // Creating variables to deal with input arguments
//  il::String inputFileName, outputDirectory, restartFileName;
//  bool checkInput = false;
//  bool checkOutput = false;
//  bool checkRestart = false;
//
//  hfp2d::loadArguments(argc, argv, il::io,
//                       checkInput, inputFileName,
//                       checkRestart, restartFileName,
//                       checkOutput, outputDirectory);
//
//
//
//  hfp2d::Mesh theMesh;
//  hfp2d::Properties matProperties;
//  hfp2d::InSituStress simConditions;
//  hfp2d::Sources sources;
//  hfp2d::SimulationParameters simParameters;

//  /// TAKING CARE OF THE INPUT IN CASE OF NEW ANALYSIS OR RESTART
//  if (checkInput) {
//
//    hfp2d::loadInput(inputFileName,
//                     il::io,
//                     theMesh,
//                     matProperties,
//                     simConditions,
//                     sources,
//                     simParameters);
//
//  } else if (checkRestart) {
//    //
//    hfp2d::loadRestart(resetFileName,il::io_t,Mesh,Properties,SimulationParam,Solution);
//  }
//
//  if (checkOutput) { // Eliminate this once a script for the output is done
//
////    // Example output to DUMMY file
////    il::String outputFile = il::join(outputDirectory, "/",
///"cracklength.txt");
////    std::cout << outputDirectory << std::endl;
////    std::cout << outputFile << std::endl;
////
////    std::ofstream foutlc;
////    foutlc.open(outputFile.asCString());
////    for (int i = 0; i < 10; i++) {
////      foutlc << i << " again good" << "\n";
////    }
////    foutlc << "Good bye, once again" << "\n";
////    foutlc.close();
//
//
//    // prepare data for output
//
//  }

//////////////////////// Previous code snippet /////////////////////////////////

//#include <il/Array.h>
//#include <src/core/Mesh.h>

//#include <src/input/json/loadJsonMesh.h>
//
//using json = nlohmann::json;

int main() {
  //  std::cout << "\n\n ----- Simple Griffith crack examples ----- \n\n" <<
  //  std::endl;

  int nelts = 5;
  double dist = 1e8;
  //

  std::cout << "-----------------" << std::endl;
  std::cout << "Parallel HFs test" << std::endl;
  std::string filename = "../Debug/ParallelHFTestsMvertex.json";

    int ret = hfp2d::ParallelHFs(filename);

  ////  il::Array<double> w{10};
  //
  //  std::cout << "\n rel error L2 norm in Linear Elements: " << ret1 << "\n";
  //  std::cout << "\n rel error L2 norm in Constant Elements (with tip
  //  correction): " << ret2 << "\n";

//  std::cout << "------------------" << std::endl;
//  std::cout << "Wellbore flow test" << std::endl;
//  int test= hfp2d::WellboreFlowBenchmark();

//  std::cout << "---------------------" << std::endl;
//  std::cout << "Well-HF coupling test" << std::endl;
//  int test = hfp2d::MultipleFracsPropagation();

  //std::cout << " end of code \n\n\n";

//   std::string dir = "../Debug/";
//   std::string meshfilename = dir + "TestMesh.json";
//
//   hfp2d::Mesh mymesh=hfp2d::loadJsonMesh(meshfilename);

  std::cout << "end of code" << std::endl;


  return 0;
}
////////////////////////////////////////////////////////////////////////////////
