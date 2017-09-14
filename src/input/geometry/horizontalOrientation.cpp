//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "horizontalOrientation.h"

namespace hfp2d {

////////////// HORIZONTAL MESH //////////////

Mesh horizontalOrientationMesh(const il::String &inputFileName,
                               const il::int_t fractureID,
                               const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {


  double x_c = findXC(inputFileName, fractureID, autoCreationMap);
  double y_c = findYC(inputFileName, fractureID, autoCreationMap);
  double length = findLength(inputFileName, fractureID, autoCreationMap);
  il::int_t numElements = findNumElem(inputFileName, fractureID, autoCreationMap);
  il::int_t interpOrder = findInterpOrder(inputFileName, fractureID, autoCreationMap);
  il::String sourceLocation = findSource(inputFileName, fractureID, autoCreationMap);


  //// Here we check for optional arguments to the automatic generation of mesh
  // in particular we are dealing with:
  // - materialID
  // - farFieldStressID
  // - porePressCondID

  il::int_t materialID = findMaterialID(inputFileName, fractureID, autoCreationMap);
  il::int_t farFieldID = findFarFieldID(inputFileName, fractureID, autoCreationMap);
  il::int_t porePresID = findPorePresID(inputFileName, fractureID, autoCreationMap);

  // create coordinates and connectivity matrices for the mesh
  il::Array2D<double> nodesCoordinates = createHorizontalMesh(x_c, y_c, length, numElements, interpOrder);
  il::Array2D<il::int_t> elementsConnectivity = createAutoConnectivity(interpOrder, numElements);
  il::Array2D<il::int_t> displ_dof_handle = createAutoDisplacementDofHandle(interpOrder, numElements);
  il::Array2D<il::int_t> press_dof_handle = createAutoPressureDofHandle(interpOrder, numElements);


///// Create Mesh

  Mesh theMesh(interpOrder,
               nodesCoordinates,
               elementsConnectivity,
               displ_dof_handle,
               press_dof_handle,
               fractureID,
               materialID,
               farFieldID,
               porePresID,
               sourceLocation);

  return theMesh;
}
}