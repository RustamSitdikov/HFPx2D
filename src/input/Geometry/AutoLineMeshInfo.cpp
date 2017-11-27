//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 13.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "AutoLineMeshInfo.h"


// (x_1, y_1) and (x_2, y_2) are the start and the end point of the line to be meshed
double findX1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  double x_1;
  keyFound = autoCreationMap.search(il::toString("x_1"));
  if (autoCreationMap.found(keyFound)) {
    x_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return x_1;
}

double findY1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double y_1;
  keyFound = autoCreationMap.search(il::toString("y_1"));
  if (autoCreationMap.found(keyFound)) {
    y_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return y_1;
}


double findX2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double x_2;
  keyFound = autoCreationMap.search(il::toString("x_2"));
  if (autoCreationMap.found(keyFound)) {
    x_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return x_2;
}

double findY2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double y_2;
  keyFound = autoCreationMap.search(il::toString("y_2"));
  if (autoCreationMap.found(keyFound)) {
    y_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_2 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return y_2;
}


// Find number of elements for the mesh
il::int_t findNumElem(const il::String &inputFileName,
                      const il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  il::int_t numElements;
  keyFound = autoCreationMap.search("number_of_elements");
  if (autoCreationMap.found(keyFound)) {
    numElements = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the number of elements in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return numElements;
}


// Find the material ID of the mesh
il::int_t findMaterialID(const il::String &inputFileName,
                         const il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  // Material ID is not necessary. If it is not set, then it is zero.
  il::int_t materialID;
  // if the material_ID keyword is found
  keyFound = autoCreationMap.search("material_ID");
  if (autoCreationMap.found(keyFound)) {
    // save its value
    materialID = autoCreationMap.value(keyFound).toInteger();
  } else {
    // save zero as default value
    materialID = 0;
  }

  return materialID;
}

