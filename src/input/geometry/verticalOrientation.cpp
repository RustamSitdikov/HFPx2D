//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "verticalOrientation.h"

namespace hfp2d {

////////////// VERTICAL MESH //////////////

void verticalOrientationMesh(const il::String &inputFileName,
                             const il::int_t &fractureID,
                             const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                             il::io_t,
                             Mesh &theMesh) {

  il::int_t keyFound;

  double x_c;
  keyFound = autoCreationMap.search(il::toString("x_c"));
  if (autoCreationMap.found(keyFound)) {
    x_c = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_c in vertical automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  double y_c;
  keyFound = autoCreationMap.search(il::toString("y_c"));
  if (autoCreationMap.found(keyFound)) {
    y_c = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_c in vertical automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  double length;
  keyFound = autoCreationMap.search(il::toString("length"));
  if (autoCreationMap.found(keyFound)) {
    length = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing length in vertical automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  il::int_t numElements;
  keyFound = autoCreationMap.search(il::toString("number_of_elements"));
  if (autoCreationMap.found(keyFound)) {
    numElements = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the number of elements in vertical automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  //// Here we check for optional arguments to the automatic generation of mesh
  // in particular we are dealing with:
  // - materialID
  // - farFieldStressID
  // - porePressCondID

  il::int_t materialID;
  keyFound = autoCreationMap.search(il::toString("material_ID"));
  if (autoCreationMap.found(keyFound)) {
    materialID = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the material ID in vertical automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

// TODO: connect with the intersection check or remove
//  il::int_t joinedWith;
//  keyFound = autoCreationMap.search(il::toString("joined_with_layer"));
//  if (autoCreationMap.found(keyFound)) {
//    joinedWith = autoCreationMap.value(keyFound).toInteger();
//  } // optional argument

///// Create Mesh
  createVerticalMesh(x_c, y_c, length, numElements, materialID, il::io, theMesh);

//  check intersections if needed, TBD
//  groupLayers.append({layerID,joinedWith});

}
}