//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 08.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "autoLineMeshCreation.h"

namespace hfp2d {

////////////// AUTOMATIC LINE MESH CREATION //////////////

Mesh autoLineMesh(const il::String &inputFileName,
                  const il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap,
                  const il::int_t interpOrder) {

  // Find start and end coordinates of line mesh
  double x_1 = findDouble("x_1", autoCreationMap, inputFileName);
  double y_1 = findDouble("y_1", autoCreationMap, inputFileName);
  double x_2 = findDouble("x_2", autoCreationMap, inputFileName);
  double y_2 = findDouble("y_2", autoCreationMap, inputFileName);

  // Recover number of elements
  il::int_t numElements = findInteger("number_of_elements", autoCreationMap, inputFileName);

  // Recover material ID
  il::int_t materialID = findInteger("material_ID", autoCreationMap, inputFileName);

  // create coordinates and connectivity matrices for the mesh
  il::Array2D<double> nodesCoordinates = createCustomMesh(x_1, y_1, x_2, y_2, numElements, interpOrder);
  il::Array2D<il::int_t> elementsConnectivity = createAutoConnectivity(interpOrder, numElements);
  il::Array2D<il::int_t> displ_dof_handle = createAutoDisplacementDofHandle(interpOrder, numElements);
  il::Array2D<il::int_t> press_dof_handle = createAutoPressureDofHandle(interpOrder, numElements);


  il::Array<il::int_t> vectorMaterialID(elementsConnectivity.size(0));
  il::Array<il::int_t> vectorFractureID(elementsConnectivity.size(0));

  // Creating the vectors of fracture ID and material ID
#pragma omp parallel for
  for(il::int_t i=0; i < elementsConnectivity.size(0); i++)
  {
    vectorMaterialID[i]=materialID;
    vectorFractureID[i]=fractureID;
  }

///// Create Mesh

  Mesh theMesh(interpOrder, nodesCoordinates, elementsConnectivity,
               displ_dof_handle, press_dof_handle, vectorFractureID, vectorMaterialID);

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