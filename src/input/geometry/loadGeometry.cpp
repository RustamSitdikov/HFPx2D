//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "loadGeometry.h"

namespace hfp2d {

Mesh loadGeometry(const il::String &inputFileName, const il::MapArray<il::String, il::Dynamic> &meshCreationMap) {

  Mesh theMesh;
  il::int_t keyFound;
  il::int_t interpOrder;
  il::int_t numFractures;

  // check type mesh generation
  if (meshCreationMap.found(meshCreationMap.search("automatic"))) {

    ////////// Automatic creation of the mesh //////////
    // if mesh generation is "automatic", start loading the interpolation order and the number of fractures

    // load interpolation order of the mesh
    keyFound = meshCreationMap.search("interpolation_order");

    if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isInteger()) {

      interpOrder = meshCreationMap.value(keyFound).toInteger();
      IL_EXPECT_FAST(interpOrder >= 0); // we do not want a negative interpolation order

    } else {

      std::cerr << "ERROR: missing the interpolation order in geometry." << std::endl;
      std::cerr << "file: " << inputFileName << std::endl;
      exit(2);

    }

    // load number of fractures
    keyFound = meshCreationMap.search("number_of_fractures");

    if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isInteger()) {

      numFractures = meshCreationMap.value(keyFound).toInteger();
      IL_EXPECT_FAST(numFractures > 0); // we want at least 1 fracture

    } else {

      std::cerr << "ERROR: missing the number of fractures in geometry." << std::endl;
      std::cerr << "file: " << inputFileName << std::endl;
      exit(2);

    }

    // For every fracture, search for the name "fracture" + number of fracture
    for (il::int_t fractureID = 0; fractureID < numFractures; fractureID++) {

      // now we create a string with the layer name
      const il::String fractureName = il::join("fracture", il::toString(fractureID));

      // search for the layer name
      keyFound = meshCreationMap.search(fractureName);

      // here we check that the "fracture#" is found and that it is a map array
      if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isMapArray()) {

        // we save the map array into autoCreationMap
        const il::MapArray<il::String, il::Dynamic> &autoCreationMap = meshCreationMap.value(keyFound).asMapArray();

        // and we pass it to the automatic mesher
        theMesh = autoLineMesh(inputFileName, fractureID, autoCreationMap, interpOrder);

      } else {
        std::cerr << "ERROR: missing fracture/mismatch in number of fractures." << std::endl;
        std::cerr << "fracture: " << fractureID << "file: " << inputFileName << std::endl;
        exit(2);
      }
    }
  } else if (meshCreationMap.found(meshCreationMap.search(il::toString("manual")))) {

    keyFound = meshCreationMap.search(il::toString("mesh_input_file"));

    if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isString()) {

      ////////////// MANUAL MESH PROVIDED BY THE USER //////////////

      // every key in TOML can be mapped to il::String
      const il::String meshFileName = meshCreationMap.value(keyFound).asString();

      // load the custom mesh
      loadMeshFile(meshFileName, il::io, theMesh);

    } else {
      std::cerr << "ERROR: manual mode in mesh but not file specified." << std::endl;
      std::cerr << "file: " << inputFileName << std::endl;
      exit(2);
    }

  } else {
    std::cerr << "ERROR: mesh generation not supported." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(2);
  }

  return theMesh;
}
}