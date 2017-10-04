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
               InSituStress &theConditions,
               Sources &theSources,
               Simulation &simParameters) {

  ///  **** Read the input data from TOML input file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(inputFileName, il::io, status);
  status.abortOnError();

  il::int_t keyFound;

  ////////// GEOMETRY KEYWORD //////////
  keyFound = config.search("geometry");

  // If "geometry" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theMesh = loadGeometry(inputFileName, meshCreationMap);

  } else {
    std::cerr << "ERROR: Geometry not found in input file " << inputFileName << std::endl;
    exit(2);
  }


  ////////// Materials: PROPERTIES KEYWORD //////////
  keyFound = config.search("properties");

  // If "solid" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &PropertiesMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theProperties=loadProperties(theMesh, inputFileName, PropertiesMap);

    // Once the properties are ready, create a vectors of properties
    // at each collocation point (for stress related properties)
    // and each node (for pore pressure/flow related properties)
    // included the in-situ conditions
  } else {
    std::cerr << "ERROR: Solid properties not found in input file " << inputFileName << std::endl;
    exit(3);
  }

  ////////// Conditions: IN-SITU KEYWORD //////////

  /// Create a vector as in the case of the materials for the materialID
  /// then create a collocation point vector for the stress distribution
  /// and a nodal vector for the pore pressure distribution
  /// TO BE CHECKED if the pore pressure is needed only at nodes

  keyFound = config.search("in-situ_conditions");

  // If "solid" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &conditionsMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theConditions=loadConditions(theMesh, inputFileName, conditionsMap);

    // Once the properties are ready, create a vectors of properties
    // at each collocation point (for stress related properties)
    // and each node (for pore pressure/flow related properties)
    // included the in-situ conditions
  } else {
    std::cerr << "ERROR: Solid properties not found in input file " << inputFileName << std::endl;
    exit(3);
  }

  ////////// Conditions: INJECTION KEYWORD //////////
  keyFound = config.search("injection");

  // If "solid" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // Save the geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &injectionMap = config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    theSources=loadSources(theMesh, inputFileName, injectionMap);

    // Once the properties are ready, create a vectors of properties
    // at each collocation point (for stress related properties)
    // and each node (for pore pressure/flow related properties)
    // included the in-situ conditions
  } else {
    std::cerr << "ERROR: Solid properties not found in input file " << inputFileName << std::endl;
    exit(3);
  }

  /// The injection take place at a number of locations. To make it as
  /// easy as possible, givem the location, we consider the closest nodes
  /// to the injection point.
  /// A cool check would be if the distance is less than the segment size.
  /// However in that case the segment size should be available for every
  /// element.

  ////////// Solution: STRATEGY KEYWORD //////////
  // placeholder

  /// Only one solver for the moment, so we just load the tolerances and
  /// the maximum number of iterations for each loop.

  ////////// Solution: OUTPUT KEYWORD //////////
  // placeholder

  /// it should load which variable should be printed



}

}