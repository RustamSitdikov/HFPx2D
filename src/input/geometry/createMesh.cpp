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
il::Array2D<double> createVerticalMesh(const double x_c,            // Center of line, x coordinate
                                       const double y_c,            // Center of line, y coordinate
                                       const double length,         // Length of line
                                       const il::int_t numElements, // Number of elements to be generated
                                       const il::int_t interpOrder) // Element interpolation order
{

  il::Array2D<double> nodesCoordinates;

  // first compute beginning and end of the mesh
  double x_1 = x_c;
  double y_1 = y_c - length / 2.0;
  double x_2 = x_c;
  double y_2 = y_c + length / 2.0;

  nodesCoordinates = createCustomMesh(x_1, y_1, x_2, y_2, numElements, interpOrder);

  return nodesCoordinates;
};

////////////// HORIZONTAL MESH //////////////
il::Array2D<double> createHorizontalMesh(const double x_c,            // Center of line, x coordinate
                                         const double y_c,            // Center of line, y coordinate
                                         const double length,         // Length of line
                                         const il::int_t numElements, // Number of elements to be generated
                                         const il::int_t interpOrder) // Element interpolation order
{

  il::Array2D<double> nodesCoordinates;

  // first compute beginning and end of the mesh
  double x_1 = x_c - length / 2.0;
  double y_1 = y_c;
  double x_2 = x_c + length / 2.0;
  double y_2 = y_c;

  nodesCoordinates = createCustomMesh(x_1, y_1, x_2, y_2, numElements, interpOrder);

  return nodesCoordinates;

};

////////////// DIAGONAL MESH //////////////
il::Array2D<double> createDiagonalMesh(const double x_c,            // Center of line, x coordinate
                                       const double y_c,            // Center of line, y coordinate
                                       const double angle,          // Angle of line
                                       const double length,         // Length of line
                                       const il::int_t numElements, // Number of elements to be generated
                                       const il::int_t interpOrder) // Element interpolation order
{

  il::Array2D<double> nodesCoordinates;

  // first compute beginning and end of the mesh
  double delta_x = length * cos(angle);
  double delta_y = length * sin(angle);

  double x_1 = x_c - delta_x;
  double y_1 = y_c - delta_y;
  double x_2 = x_c + delta_x;
  double y_2 = y_c + delta_y;

  nodesCoordinates = createCustomMesh(x_1, y_1, x_2, y_2, numElements, interpOrder);

  return nodesCoordinates;

};

////////////// CUSTOM MESH //////////////
il::Array2D<double> createCustomMesh(const double x_1,             // First point of the line, x coordinate
                                     const double y_1,             // First point of the line, y coordinate
                                     const double x_2,             // Second point of the line, x coordinate
                                     const double y_2,             // Second point of the line, y coordinate
                                     const il::int_t numElements,  // Number of elements to be generated
                                     const il::int_t interpOrder)  // Element interpolation order
{
  il::int_t number_of_nodes = numElements * interpOrder + 1;
  il::Array2D<double> nodesCoordinates(number_of_nodes, 2);

  double delta_x = (x_2 - x_1) / numElements / interpOrder;
  double delta_y = (y_2 - y_1) / numElements / interpOrder;

  for (il::int_t i = 0; i < number_of_nodes; i++) {

    nodesCoordinates(i, 0) = x_1 + delta_x * i;
    nodesCoordinates(i, 1) = y_1 + delta_y * i;

  }

  return nodesCoordinates;

};


////////////// CREATE MESH CONNECTIVITY //////////////

il::Array2D<il::int_t> createAutoConnectivity(const il::int_t interpolationOrder,
                                              const il::int_t numberOfElements) {
  il::Array2D<il::int_t> connectivityMatrix(numberOfElements, interpolationOrder + 1);

  for (il::int_t i = 0; i < numberOfElements; i++) {
    for (il::int_t j = 0; j < interpolationOrder + 1; j++) {
      connectivityMatrix(i, j) = i * interpolationOrder + j;
    }
  }

  return connectivityMatrix;
}

////////////// CREATE DISPLACEMENTS DOF HANDLES //////////////

il::Array2D<il::int_t> createAutoDisplacementDofHandle(const il::int_t interpolationOrder,
                                                       const il::int_t numberOfElements) {

  il::Array2D<il::int_t> dof_handle_displacement(numberOfElements, 2 * (interpolationOrder + 1))

// filling the dof_handle_displacements (2D discontinuous galerkin)
  for (il::int_t i = 0, j; i < numberOfElements; i++) {

    j = i * 2 * (interpolationOrder + 1);

    for (int k = 0; k < 2 * (interpolationOrder + 1); k++) {
      dof_handle_displacement(i, k) = j + k;
    }

  }
  return dof_handle_displacement;
}

////////////// CREATE PRESSURE DOF HANDLES //////////////

il::Array2D<il::int_t> createAutoPressureDofHandle(const il::int_t interpolationOrder,
                                                   const il::int_t numberOfElements) {

  il::Array2D<il::int_t> dof_handle_pressure(numberOfElements, interpolationOrder + 1);

// fillind the dof_handle_pressure (continuous galerkin)
  for (il::int_t i = 0; i < numberOfElements; i++) {
    for (il::int_t j = 0; j < numberOfElements; j++) {
      dof_handle_pressure(i, j) = i * interpolationOrder + j;
    }
  }

  return dof_handle_pressure;
}

}