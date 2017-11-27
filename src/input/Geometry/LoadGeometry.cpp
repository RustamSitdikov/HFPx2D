//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

// Inclusion from the project
#include "LoadGeometry.h"

namespace hfp2d {

Mesh loadGeometry(
    const il::String &inputFileName,
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap) {
  Mesh theGlobalMesh;
  Mesh theSingleMesh;
  il::int_t keyFound;
  il::int_t interpOrder;
  il::int_t numFractures;

  keyFound = meshCreationMap.search("mesh_generation");

  // Check existence of mesh generation keyword
  if (meshCreationMap.found(keyFound)) {
    ////////// Automatic creation of the mesh //////////
    // if mesh generation is "automatic", start loading the interpolation order
    // and the number of fractures

    if (meshCreationMap.value(keyFound).asString() == "automatic") {
      interpOrder = findInterpOrder(meshCreationMap, inputFileName);

      numFractures = findNumFractures(meshCreationMap, inputFileName);

      // For every fracture, search for the name "fracture" + number of fracture
      for (il::int_t fractureID = 1, number_fractures = 1;
           fractureID <= numFractures; fractureID++) {
        // now we create a string with the layer name
        const il::String fractureName =
            il::join("fracture", il::toString(fractureID));

        // search for the layer name
        keyFound = meshCreationMap.search(fractureName);

        // here we check that the "fracture#" is found and that it is a map
        // array
        if (meshCreationMap.found(keyFound) &&
            meshCreationMap.value(keyFound).isMapArray()) {
          // we save the map array into autoCreationMap
          const il::MapArray<il::String, il::Dynamic> &autoCreationMap =
              meshCreationMap.value(keyFound).asMapArray();

          // and we pass it to the automatic mesher
          theSingleMesh = autoLineMesh(inputFileName, fractureID,
                                       autoCreationMap, interpOrder);

//          std::cout << " Mesh number " << fractureID << " has "
//                    << number_fractures << " fractures" << std::endl;
          if (theGlobalMesh.numberOfElts() == 0) {
            theGlobalMesh = theSingleMesh;

          } else {
            number_fractures = number_fractures + fractureID;
            // theGlobalMesh.appendMesh(theSingleMesh,false); TODO: appendMesh
            // method
          }

//          std::cout << " after fracture " << fractureID
//                    << " the global mesh has " << number_fractures
//                    << " fractures" << std::endl;
        } else {
          std::cerr
              << "ERROR: missing fracture/mismatch in number of fractures."
              << std::endl;
          std::cerr << "fracture: " << fractureID << "file: " << inputFileName
                    << std::endl;
          exit(2);
        }
      }
    } else if (meshCreationMap.value(keyFound).asString() == "manual") {
      keyFound = meshCreationMap.search(il::toString("mesh_input_file"));

      if (meshCreationMap.found(keyFound) &&
          meshCreationMap.value(keyFound).isString()) {
        ////////////// MANUAL MESH PROVIDED BY THE USER //////////////

        // every key in TOML can be mapped to il::String
        const il::String meshFileName =
            meshCreationMap.value(keyFound).asString();

        // load the custom mesh
        loadMeshFile(meshFileName, il::io, theGlobalMesh);

      } else {
        std::cerr << "ERROR: manual mode in mesh but not file specified."
                  << std::endl;
        std::cerr << "file: " << inputFileName << std::endl;
        exit(2);
      }

    } else {
      std::cerr << "ERROR: mesh generation not supported." << std::endl;
      std::cerr << "file: " << inputFileName << std::endl;

      exit(2);
    }

  } else {
    std::cerr << "ERROR: no mesh generation specified." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(2);
  }
  return theGlobalMesh;
}

/////////////////////////

il::int_t findInterpOrder(
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap,
    const il::String &inputFileName) {
  il::int_t keyFound;
  il::int_t interpOrder;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("interpolation_order");

  if (meshCreationMap.found(keyFound) &&
      meshCreationMap.value(keyFound).isInteger()) {
    interpOrder = meshCreationMap.value(keyFound).toInteger();
    IL_EXPECT_FAST(interpOrder >=
                   0);  // we do not want a negative interpolation order

  } else {
    std::cerr << "ERROR: missing the interpolation order in Geometry."
              << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(2);
  }

  return interpOrder;
}

////////////////////////////////

il::int_t findNumFractures(
    const il::MapArray<il::String, il::Dynamic> &meshCreationMap,
    const il::String &inputFileName) {
  il::int_t keyFound;
  il::int_t numFractures;
  keyFound = meshCreationMap.search("number_of_fractures");

  if (meshCreationMap.found(keyFound) &&
      meshCreationMap.value(keyFound).isInteger()) {
    numFractures = meshCreationMap.value(keyFound).toInteger();

    IL_EXPECT_FAST(numFractures > 0);  // we want at least 1 fracture

  } else {
    std::cerr << "ERROR: missing the number of fractures in Geometry."
              << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(2);
  }

  return numFractures;
}
}