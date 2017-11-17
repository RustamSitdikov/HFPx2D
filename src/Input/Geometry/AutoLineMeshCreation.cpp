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
#include "AutoLineMeshCreation.h"

namespace hfp2d {

////////////// AUTOMATIC LINE MESH CREATION //////////////

Mesh autoLineMesh(const il::String &inputFileName, const il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                  const il::int_t interpOrder) {
  // Find start and end coordinates of line mesh
  double x_1 = findDouble("x_1", autoCreationMap, inputFileName);
  double y_1 = findDouble("y_1", autoCreationMap, inputFileName);
  double x_2 = findDouble("x_2", autoCreationMap, inputFileName);
  double y_2 = findDouble("y_2", autoCreationMap, inputFileName);

  // Recover number of elements
  il::int_t numElements =
      findInteger("number_of_elements", autoCreationMap, inputFileName);

  // Recover material ID
  il::int_t materialID =
      findInteger("material_ID", autoCreationMap, inputFileName);

  // create coordinates and connectivity matrices for the mesh
  il::Array2D<double> nodesCoordinates =
      createCustomMesh(x_1, y_1, x_2, y_2, numElements, interpOrder);
  il::Array2D<il::int_t> elementsConnectivity =
      createAutoConnectivity(interpOrder, numElements);
  il::Array2D<il::int_t> displ_dof_handle =
      createAutoDisplacementDofHandle(interpOrder, numElements);
  il::Array2D<il::int_t> press_dof_handle =
      createAutoPressureDofHandle(interpOrder, numElements);

  il::Array<int> vectorMaterialID(elementsConnectivity.size(0));
  il::Array<int> vectorFractureID(elementsConnectivity.size(0));
  il::Array<int> vectorConditionID(elementsConnectivity.size(0));

// Creating the vectors of fracture ID and material ID
#pragma omp parallel for
  for (il::int_t i = 0; i < elementsConnectivity.size(0); i++) {
    vectorMaterialID[i] = materialID;
    vectorFractureID[i] = fractureID;
  }

  ///// Create Mesh /////

  /// Note that vectorMateriaID is a simply array of int materialID
  /// -> homogeneous material (materialID should be a vector!)
  Mesh theMesh(nodesCoordinates, elementsConnectivity, vectorMaterialID,
               interpOrder);

  return theMesh;
}
}