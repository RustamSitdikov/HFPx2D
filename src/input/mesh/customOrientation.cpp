//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "customOrientation.h"

namespace hfp2d {

////////////// CUSTOM MESH //////////////

void customOrientationMesh(const il::String &inputFileName,
                           const il::int_t &idLayer,
                           const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                           il::io_t,
                           Mesh &theMesh) {

  il::int_t keyFound;

  double x_1;
  keyFound = autoCreationMap.search(il::toString("x_c"));
  if (autoCreationMap.found(keyFound)) {
    x_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  double y_1;
  keyFound = autoCreationMap.search(il::toString("y_c"));
  if (autoCreationMap.found(keyFound)) {
    y_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_2 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  double x_2;
  keyFound = autoCreationMap.search(il::toString("x_c"));
  if (autoCreationMap.found(keyFound)) {
    x_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  double y_2;
  keyFound = autoCreationMap.search(il::toString("y_c"));
  if (autoCreationMap.found(keyFound)) {
    y_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_2 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  il::int_t numElements;
  keyFound = autoCreationMap.search(il::toString("number_of_elements"));
  if (autoCreationMap.found(keyFound)) {
    numElements = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the number of elements in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  il::int_t materialID;
  keyFound = autoCreationMap.search(il::toString("material_ID"));
  if (autoCreationMap.found(keyFound)) {
    materialID = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the material ID in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << idLayer << ", file: " << inputFileName << std::endl;
    exit(2);
  }

// TODO: connect with the intersection check or remove
//  il::int_t joinedWith;
//  keyFound = autoCreationMap.search(il::toString("joined_with_layer"));
//  if (autoCreationMap.found(keyFound)) {
//    joinedWith = autoCreationMap.value(keyFound).toInteger();
//  } // optional argument

///// Create Mesh
  createCustomMesh(x_1, y_1, x_2, y_2, numElements, materialID, il::io, theMesh);

//  check intersections if needed, TBD
//  groupLayers.append({layerID,joinedWith});

}
}