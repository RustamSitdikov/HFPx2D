//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 06.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "createMesh.h"

namespace hfp2d {

////////////// VERTICAL MESH //////////////
void createVerticalMesh(double& x_c,            // Center of line, x coordinate
                        double& y_c,            // Center of line, y coordinate
                        double& length,         // Length of line
                        il::int_t& numElements, // Number of elements to be generated
                        il::int_t& materialID,  // material ID of elements
                        il::io_t,
                        Mesh &theMesh) {        // Final mesh to be returned

  // first compute beginning and end of the mesh
  double x_1 = x_c;
  double y_1 = y_c - length / 2.0;
  double x_2 = x_c;
  double y_2 = y_c + length / 2.0;

  createCustomMesh(x_1, y_1, x_2, y_2, numElements, materialID, il::io, theMesh);

}

////////////// HORIZONTAL MESH //////////////
void createHorizontalMesh(double& x_c,            // Center of line, x coordinate
                          double& y_c,            // Center of line, y coordinate
                          double& length,         // Length of line
                          il::int_t& numElements, // Number of elements to be generated
                          il::int_t& materialID,  // material ID of elements
                          il::io_t,
                          Mesh &theMesh) {        // Final mesh to be returned

  // first compute beginning and end of the mesh
  double x_1 = x_c - length / 2.0;
  double y_1 = y_c;
  double x_2 = x_c + length / 2.0;
  double y_2 = y_c;

  createCustomMesh(x_1, y_1, x_2, y_2, numElements, materialID, il::io, theMesh);

}

////////////// DIAGONAL MESH //////////////
//theMesh = hfp2d::createDiagonalMesh(x_c,y_c,angle,length,numElements,materialID);


void createDiagonalMesh(double& x_c,            // Center of line, x coordinate
                        double& y_c,            // Center of line, y coordinate
                        double& angle,          // Angle of line
                        double& length,         // Length of line
                        il::int_t& numElements, // Number of elements to be generated
                        il::int_t& materialID,  // material ID of elements
                        il::io_t,
                        Mesh &theMesh) {        // Final mesh to be returned

  // first compute beginning and end of the mesh
  double delta_x= length * cos(angle);
  double delta_y= length * sin(angle);

  double x_1 = x_c - delta_x;
  double y_1 = y_c - delta_y;
  double x_2 = x_c + delta_x;
  double y_2 = y_c + delta_y;

  createCustomMesh(x_1, y_1, x_2, y_2, numElements, materialID, il::io, theMesh);

}


////////////// CUSTOM MESH //////////////
//theMesh = hfp2d::createCustomMesh(x_1,y_1,x_2,y_2,numElements,materialID);

void createCustomMesh(double& x_1,             // First point of the line, x coordinate
                      double& y_1,             // First point of the line, y coordinate
                      double& x_2,             // Second point of the line, x coordinate
                      double& y_2,             // Second point of the line, y coordinate
                      il::int_t& numElements,  // Number of elements to be generated
                      il::int_t& materialID,   // material ID of elements
                      il::io_t,
                      Mesh& theMesh) {         // Final mesh to be returned

  il::Array2D<double> coordinates(numElements + 1,2);
  il::Array2D<il::int_t> connectivity(numElements, 2);
  il::Array<il::int_t> matID(numElements);

  double delta_x = (x_2 - x_1) / numElements;
  double delta_y = (y_2 - y_1) / numElements;

  for (il::int_t i = 0; i < numElements; i++) {

    coordinates(i,0) = x_1 + delta_x * i;
    coordinates(i,1) = y_1 + delta_y * i;

    connectivity(i, 0) = i;
    connectivity(i, 0) = i + 1;

    matID[i] = materialID;

  }

  ////// HERE WE SAVE EVERYTHING TO THE MESH CLASS VARIABLE
  // NOTE:  we should take into account previous generated meshes, in order to avoid
  //        overwriting them. Additionally, we will have to check that we do not have
  //        repeated nodes, in particular when the joined_with flag is ON!!

  if(theMesh.numberOfElements()>0)
  {
    theMesh.init1DMesh(coordinates, connectivity, matID);
  }
  else
  {
    // TODO: check intersections here?! before appending rather than later
    theMesh.appendMesh(coordinates, connectivity, matID);
  }



}

}