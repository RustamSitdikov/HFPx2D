//
// Created by lorenzo on 9/4/17.
//

#include "loadInput.h"

#include <iostream>
#include <il/Array.h>
#include <il/toml.h>
#include "Mesh.h"
#include "SimulationParameters.h"



namespace hfp2d {

void loadInput(std::string inputFileName){
//                , il::io_t,
//                Mesh mesh,
//                Properties matProperties,
//                Simulation simParameters)


  ///  **** Read the input data from TOML input file **** ///

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto config =
      il::load<il::MapArray<il::String, il::Dynamic>>(inputFileName, il::io, status);
  status.abortOnError();

  il::int_t keyFound;
  il::int_t theInteger;
  double theDouble;
  std::string theString;

  ////////// Start searching for the geometry keyword //////////

  il::int_t numLayers; // number of layers in mesh

  keyFound = config.search("geometry");
  // If the "geometry" is found
  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

    // load number of layers
    keyFound = config.search("number_of_layers");
    if (config.found(keyFound) && config.value(keyFound).isInteger()) {

      numLayers = config.value(keyFound).toInteger();
    } else {
      std::cerr << "ERROR: missing the number of layers in geometry." << std::endl;
      std::cerr << "file: "<< inputFileName << std::endl;
      exit(2);
    }

    // check type mesh generation
    if (config.found(config.search("automatic"))) {
      // if mesh generation is "automatic", start loading the layers

      ////////// Automatic creation of the mesh //////////
      for (il::int_t idLayer = 0; idLayer < numLayers; idLayer++) {
        // TODO:::RESTART FROM HERE
        // now we create a string with the layer name
        std::string layerName = "layer" + std::to_string(idLayer);

        // search for the layer name
        keyFound = config.search(layerName);

        if (config.found(keyFound) && config.value(keyFound).isMapArray()) {

          const il::MapArray<il::String, il::Dynamic> &automaticMesh = config.value(keyFound).asMapArray();

          keyFound = config.search("name");
          keyFound = config.search("orientation"); // horizontal/vertical/diagonal/custom
          // for the corresponding orientation
          keyFound = config.search("x_c");
          keyFound = config.search("y_c");
          keyFound = config.search("angle");
          keyFound = config.search("length");
          keyFound = config.search("x_1");
          keyFound = config.search("y_1");
          keyFound = config.search("y_2");
          keyFound = config.search("x_2");
          keyFound = config.search("number_of_elements");
          keyFound = config.search("material_ID");
          keyFound = config.search("joined_with_layer")

        }
        else
        {
          std::cerr << "ERROR: missing layer/mismatch in number of layers." << std::endl;
          std::cerr << "file: "<< inputFileName << std::endl;
          exit(2);
        }
      }
    }
    else if(config.found(config.search("manual")))
    {
      keyFound = config.search("mesh_input_file");
      if (config.found(keyFound) && config.value(keyFound).isString()){

        // hfp2d::loadMeshFile(config.value(keyFound)); // every key in TOML is first a string
        // TODO: complete with total mesh load

      } else {
        std::cerr << "ERROR: manual mode in mesh but not file specified." << std::endl;
        std::cerr << "file: "<< inputFileName << std::endl;
        exit(2);
      }
      } else {
      std::cerr << "ERROR: mesh generation not supported." << std::endl;
      std::cerr << "file: "<< inputFileName << std::endl;
      exit(2);
    }


  } else {
    std::cerr << "ERROR: Geometry not found in input file " << inputFileName << std::endl;
    exit(2);
  }

  // next cards
}