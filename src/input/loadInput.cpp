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
               Properties &theProperties,
               Simulation &simParameters) {

  ///  **** Read the input data from TOML input file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(inputFileName, il::io, status);
  status.abortOnError();

  il::int_t keyFound;
  // il::int_t theInteger;
  // double theDouble;
  // il::String theString;

/*
  std::cout << "passato input" << std::endl;
  std::cout << inputFileName << std::endl;
*/
  ////////// GEOMETRY KEYWORD //////////
  keyFound = config.search("geometry");

  // If "geometry" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theMesh = loadGeometry(inputFileName, meshCreationMap);

    // check of the mesh
//    for(il::int_t i=0; i < theMesh.numberOfNodes(); i++){
//      std :: cout << i << " " << theMesh.node(i,0) << " " << theMesh.node(i,1) <<std::endl;
//    }
//    for(il::int_t i=0; i < theMesh.numberOfElements(); i++){
//      std :: cout << i << " " << theMesh.matID(i) << " " << theMesh.connectivity(i,0) << " " << theMesh.connectivity(i,1) <<std::endl;
//    }
//    for(il::int_t i=0; i < theMesh.numberOfElements(); i++){
//      std :: cout << i << " " << theMesh.fracID(i) << " " << theMesh.dofDispl(i,0) << " "
//                              << theMesh.dofDispl(i,1) << " "
//                              << theMesh.dofDispl(i,2) << " "
//                              << theMesh.dofDispl(i,3) << " " <<std::endl;
//    }
//    for(il::int_t i=0; i < theMesh.numberOfElements(); i++){
//      std :: cout << i << " " << theMesh.fracID(i) << " " << theMesh.dofPress(i,0) << " "
//                              << theMesh.dofPress(i,1) << " " <<std::endl;
//    }


  } else {
    std::cerr << "ERROR: Geometry not found in input file " << inputFileName << std::endl;
    exit(2);
  }


  ////////// Materials: PROPERTIES KEYWORD //////////
  il::int_t numOfMats=theMesh.numberOfMaterials();

  keyFound = config.search("properties");

  // If "solid" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &PropertiesMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theProperties=loadProperties(inputFileName, PropertiesMap, numOfMats);

    // Once the properties are ready, create a vectors of properties
    // at each collocation point (for stress related properties)
    // and each node (for pore pressure/flow related properties)
    // included the in-situ conditions
  } else {
    std::cerr << "ERROR: Solid properties not found in input file " << inputFileName << std::endl;
    exit(3);
  }

  ////////// Conditions: INJECTION KEYWORD //////////
  // placeholder

  ////////// Solution: STRATEGY KEYWORD //////////
  // placeholder


  // next cards

}

}