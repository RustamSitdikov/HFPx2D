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
#include <cmath>
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
//
int main(int const argc, char const* const* argv) {

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
  hfp2d::Conditions simConditions;
  hfp2d::Sources sources;
  hfp2d::Simulation simParameters;

  /// TAKING CARE OF THE INPUT IN CASE OF NEW ANALYSIS OR RESTART
  if (checkInput) {

    hfp2d::loadInput(inputFileName,
                     il::io,
                     theMesh,
                     matProperties,
                     simConditions,
                     sources,
                     simParameters);

  } else if (checkRestart) {
    // hfp2d::loadRestart(resetFileName,il::io_t,Mesh,Properties,SimulationParam,Solution);
  }

  if (checkOutput) { // Eliminate this once a script for the output is done

//    // Example output to DUMMY file
//    il::String outputFile = il::join(outputDirectory, "/", "cracklength.txt");
//    std::cout << outputDirectory << std::endl;
//    std::cout << outputFile << std::endl;
//
//    std::ofstream foutlc;
//    foutlc.open(outputFile.asCString());
//    for (int i = 0; i < 10; i++) {
//      foutlc << i << " again good" << "\n";
//    }
//    foutlc << "Good bye, once again" << "\n";
//    foutlc.close();


    // prepare data for output

  }

//////////////////////// Prepare data for computation /////////////////////////////////

  // create the source vector for displacement+pressure dofs
  //il::int_t totalNumDofs= theMesh.numberOfDisplDofsPerElement()+ theMesh.numberOfPressDofsPerElement();
  // the source vector (or forcing vector) will be created after we checked for the position of the injection
  // there will be a method in the class to give that vector

  // create the stiffness matrix for the computation
  //il::Array2D<double> kmat = basic_assembly_new(mesh_total, id, p, material.Ep);
  //il::Array2D<double> kmat = basic_assembly_new( mesh, id, p, Ep);  // from Dong
  // how I would like it
  // kmat = basic_assembly (mesh, properties, kernel);

//////////////////////// Initiate the computational loop /////////////////////////////////

//////////////////////// Previous code snippet /////////////////////////////////

  std::cout << "\n\n ----- Simple Griffith crack examples ----- \n\n" << std::endl;

  int nelts = 10;

  double ret1 = hfp2d::SimpleGriffithExampleLinearElement(nelts);
  double ret2 = hfp2d::SimpleGriffithExampleS3D_P0(nelts);

  std::cout << "\n rel error L2 norm in Linear Elements: " << ret1 << "\n";
  std::cout << "\n rel error L2 norm in Constant Elements: " << ret2 << "\n";

  std::cout << " end of code \n\n\n";

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

