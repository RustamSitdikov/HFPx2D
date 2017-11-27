//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 04.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

// Inclusion from the project
#include "LoadInput.h"

namespace hfp2d {

void loadInput(const il::String &input_filename, il::io_t, Mesh &MyMesh,
               ElasticProperties &ElasticProperties,
               FluidProperties &FluidProperties, SolidEvolution &SolidEvolution,
               FractureEvolution &FractureEvolution,
               InSituStress &BackgroundLoadingConditions) {
  ///  **** Read the input data from TOML input file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(input_filename,
                                                                il::io, status);
  status.abortOnError();

  il::int_t keyFound;

  ////////// GEOMETRY KEYWORD //////////
  keyFound = config.search("Geometry");

  // If "Geometry" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Geometry map, i.e. the data to create the mesh
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap =
        config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    MyMesh = loadGeometry(input_filename, meshCreationMap);

  } else {
    std::cerr << "ERROR: 'Geometry' not found in input file " << input_filename
              << std::endl;
    exit(EXIT_FAILURE);
  }

  ////////// Materials: PROPERTIES KEYWORD //////////
  keyFound = config.search("Properties");

  // If "Properties" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Properties map
    const il::MapArray<il::String, il::Dynamic> &PropertiesMap =
        config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    loadProperties(MyMesh, input_filename, PropertiesMap, ElasticProperties,
                   FluidProperties, SolidEvolution, FractureEvolution);

    // Once the Properties are ready, create a vectors of Properties
    // at each collocation point (for stress related Properties)
    // and each node (for pore pressure/flow related Properties)
    // included the in-situ Conditions
  } else {
    std::cerr << "ERROR: 'Properties' not found in input file "
              << input_filename << std::endl;
    exit(EXIT_FAILURE);
  }

  ////////// Materials: PROPERTIES KEYWORD //////////
  keyFound = config.search("Loading_Conditions");

  // If "Properties" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Properties map
    const il::MapArray<il::String, il::Dynamic> &ConditionsMap =
        config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadGeometry script
    BackgroundLoadingConditions =
        loadConditions(MyMesh, input_filename, ConditionsMap);

    // Once the Properties are ready, create a vectors of Properties
    // at each collocation point (for stress related Properties)
    // and each node (for pore pressure/flow related Properties)
    // included the in-situ Conditions
  } else {
    std::cerr << "ERROR: 'Loading Conditions' not found in input file "
              << input_filename << std::endl;
    exit(EXIT_FAILURE);
  }
}
}