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

Mesh verticalOrientationMesh(const il::String &inputFileName,
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
  il::Array2D<double> nodesCoordinates = createVerticalMesh(x_c, y_c, length, numElements, interpOrder);
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

//// TODO: connect with the intersection - check or remove
//
//////// HERE WE SAVE EVERYTHING TO THE MESH CLASS VARIABLE
//// NOTE:  we should take into account previous generated meshes, in order to avoid
////        overwriting them. Additionally, we will have to check that we do not have
////        repeated nodes, in particular when the joined_with flag is ON!!
//
//if(theMesh.numberOfElements()>0)
//{
//theMesh=Mesh(1,coordinates, connectivity, sourceID);
//}
//else
//{
//// TODO: check intersections here?! before appending rather than later
//// theMesh.appendMesh(coordinates, connectivity, matID);
//}
//
//// TODO: ADD THE DOF_HANDLES IN THE MESH
//// TODO: TWO DIFFERENT FRACS SHOULD HAVE TWO DIFFERENT FRAC_ID, always
//// so when we add a new one we do not have troubles with the tips
//// TODO: ADD THE IS_TIP BOOLEAN VECTOR