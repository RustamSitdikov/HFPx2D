//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 04.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "loadInput.h"

namespace hfp2d {

void loadInput(const il::String &inputFileName,
               il::io_t,
               Mesh &theMesh,
               Properties &matProperties,
               Simulation &simParameters)
{

  ///  **** Read the input data from TOML input file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
//  auto config =
//      il::load<il::MapArray<il::String, il::Dynamic>>(inputFileName, il::io, status);
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(inputFileName, il::io, status);
  status.abortOnError();

  il::int_t keyFound;
  // il::int_t theInteger;
  // double theDouble;
  // il::String theString;

  ////////// GEOMETRY KEYWORD //////////
  keyFound = config.search("geometry");

  // If "geometry" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    loadGeometry(inputFileName,
            meshCreationMap,
            il::io,
            theMesh);

  } else {
    std::cerr << "ERROR: Geometry not found in input file " << inputFileName << std::endl;
    exit(2);
  }

  ////////// Materials: SOLID KEYWORD //////////
  // placeholder

  ////////// Materials: FLUID KEYWORD //////////
  // placeholder

  ////////// Conditions: IN-SITU KEYWORD //////////
  // placeholder

  ////////// Conditions: FAULT KEYWORD //////////
  // placeholder

  ////////// Conditions: INJECTION KEYWORD //////////
  // placeholder

  ////////// Solution: STRATEGY KEYWORD //////////
  // placeholder


  // next cards
}
}