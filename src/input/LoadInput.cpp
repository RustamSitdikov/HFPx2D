//
// HFPx2D project.
//
// Created by Federico Ciardo on 04.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

// Inclusion from the project
#include "LoadInput.h"

namespace hfp2d {

void loadInput(const il::String &config_filename, il::io_t, Mesh &MyMesh,
               ElasticProperties &ElasticProperties,
               FluidProperties &FluidProperties, SolidEvolution &SolidEvolution,
               FractureEvolution &FractureEvolution,
               InSituStress &BackgroundLoadingConditions,
               double &const_overpress, double &t_0plus1, double &time_step,
               double &final_time, bool &expl_impl, bool &quasi_dynamic) {
  ///  **** Read the input data from TOML configuration file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(config_filename,
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
    MyMesh = loadGeometry(config_filename, meshCreationMap);

  } else {
    std::cerr << "ERROR: 'Geometry' not found in input file " << config_filename
              << std::endl;
    exit(EXIT_FAILURE);
  }

  ////////// PROPERTIES KEYWORD //////////
  keyFound = config.search("Properties");

  // If "Properties" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Properties map
    const il::MapArray<il::String, il::Dynamic> &PropertiesMap =
        config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadProperties script
    loadProperties(MyMesh, config_filename, PropertiesMap, ElasticProperties,
                   FluidProperties, SolidEvolution, FractureEvolution);

  } else {
    std::cerr << "ERROR: 'Properties' not found in input file "
              << config_filename << std::endl;
    exit(EXIT_FAILURE);
  }

  ////////// LOADING_CONDITIONS KEYWORD //////////
  keyFound = config.search("Loading_Conditions");

  // If "Loading_Conditions" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Loading_Conditions map
    const il::MapArray<il::String, il::Dynamic> &ConditionsMap =
        config.value(keyFound).asMapArray();

    // Send the data in meshCreationMap to loadConditions script
    BackgroundLoadingConditions =
        loadConditions(MyMesh, config_filename, ConditionsMap);

  } else {
    std::cerr << "ERROR: 'Loading Conditions' not found in input file "
              << config_filename << std::endl;
    exit(EXIT_FAILURE);
  }

  ////////// SIMULATION_PARAMETERS KEYWORD //////////
  keyFound = config.search("Simulation_Parameters");

  // If "Simulation_Parameters" is found and it is a map
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    // Save the Simulation_Paramaters map
    const il::MapArray<il::String, il::Dynamic> &SimulationParametersMap =
        config.value(keyFound).asMapArray();

    const_overpress = findDouble("constant_overpressure",
                                 SimulationParametersMap, config_filename);

    t_0plus1 =
        findDouble("initial_time", SimulationParametersMap, config_filename);

    time_step =
        findDouble("time_step", SimulationParametersMap, config_filename);

    final_time =
        findDouble("final_time", SimulationParametersMap, config_filename);

    expl_impl = findBool("expl_impl", SimulationParametersMap, config_filename);

    quasi_dynamic =
        findBool("QD", SimulationParametersMap, config_filename);

  } else {
    std::cerr << "ERROR: 'Simulation_Parameters' not found in input file "
              << config_filename << std::endl;
    exit(EXIT_FAILURE);
  }
}
}