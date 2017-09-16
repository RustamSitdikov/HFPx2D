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
  il::int_t numFractures;

  // load number of layers
  keyFound = meshCreationMap.search(il::toString("number_of_fractures"));
  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isInteger()) {

    numFractures = meshCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the number of layers in geometry." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(2);
  }

  // check type mesh generation
  if (meshCreationMap.found(meshCreationMap.search(il::toString("automatic")))) {
    // if mesh generation is "automatic", start loading the layers

    ////////// Automatic creation of the mesh //////////
    for (il::int_t fractureID = 0; fractureID < numFractures; fractureID++) {

      // now we create a string with the layer name
      const il::String fractureName = il::join(il::toString("fracture"), il::toString(fractureID));

      // search for the layer name
      keyFound = meshCreationMap.search(fractureName);

      if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isMapArray()) {

        const il::MapArray<il::String, il::Dynamic> &autoCreationMap = meshCreationMap.value(keyFound).asMapArray();

        //keyFound = meshCreationMap.search(il::toString("name"));

        keyFound = meshCreationMap.search(il::toString("orientation")); // horizontal/vertical/diagonal/custom

        if (autoCreationMap.found(keyFound) && autoCreationMap.value(keyFound).isString()) {
          if ((autoCreationMap.value(keyFound).asString()) == "vertical") {

            theMesh = verticalOrientationMesh(inputFileName,
                                              fractureID,
                                              autoCreationMap);

          } else if (autoCreationMap.value(keyFound).asString() == "horizontal") {

            theMesh = horizontalOrientationMesh(inputFileName,
                                                fractureID,
                                                autoCreationMap);

          } else if (autoCreationMap.value(keyFound).asString() == "diagonal") {

            theMesh = diagonalOrientationMesh(inputFileName,
                                              fractureID,
                                              autoCreationMap);

          } else if (autoCreationMap.value(keyFound).asString() == "custom") {

            theMesh = customOrientationMesh(inputFileName,
                                            fractureID,
                                            autoCreationMap);

          } else {
            std::cerr << "ERROR: orientation type not recognized." << std::endl;
            std::cerr << "file: " << inputFileName << std::endl;
            exit(2);
          }
        }
      } else {
        std::cerr << "ERROR: missing layer/mismatch in number of layers." << std::endl;
        std::cerr << "file: " << inputFileName << std::endl;
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