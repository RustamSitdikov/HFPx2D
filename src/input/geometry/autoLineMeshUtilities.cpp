//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 06.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "autoLineMeshUtilities.h"

namespace hfp2d {

////////////// CUSTOM MESH NODES CREATION //////////////
il::Array2D<double> createCustomMesh(const double x_1,             // First point of the line, x coordinate
                                     const double y_1,             // First point of the line, y coordinate
                                     const double x_2,             // Second point of the line, x coordinate
                                     const double y_2,             // Second point of the line, y coordinate
                                     const il::int_t numElements,  // Number of elements to be generated
                                     const il::int_t interpOrder)  // Element interpolation order
{

  il::int_t number_of_nodes;
  double delta_x;
  double delta_y;

  if(interpOrder == 0) {

    number_of_nodes = numElements * 2 + 1; // particular case for the constant piecewise element
    delta_x = (x_2 - x_1) / numElements;
    delta_y = (y_2 - y_1) / numElements;

  } else {

    number_of_nodes = numElements * interpOrder + 1;
    delta_x = (x_2 - x_1) / numElements / interpOrder;
    delta_y = (y_2 - y_1) / numElements / interpOrder;

  }

  il::Array2D<double> nodesCoordinates(number_of_nodes, 2);

#pragma omp parallel for
  for (il::int_t i = 0; i < number_of_nodes; i++) {

    nodesCoordinates(i, 0) = x_1 + delta_x * i;
    nodesCoordinates(i, 1) = y_1 + delta_y * i;

  }

  return nodesCoordinates;

};


////////////// CREATE MESH CONNECTIVITY //////////////

il::Array2D<il::int_t> createAutoConnectivity(const il::int_t interpolationOrder,
                                              const il::int_t numberOfElements) {

  il::int_t columns_conn_mtx;

  if(interpolationOrder == 0){
    columns_conn_mtx = 2;
  } else {
    columns_conn_mtx = interpolationOrder + 1;
  }

  // Declaration of the connectivity matrix
  il::Array2D<il::int_t> connectivityMatrix(numberOfElements, columns_conn_mtx);

#pragma omp parallel for collapse(2)
  for (il::int_t i = 0; i < numberOfElements; i++) {

    for (il::int_t j = 0; j < columns_conn_mtx; j++) {

      connectivityMatrix(i, j) = i * (columns_conn_mtx-1) + j;

    }

  }

  return connectivityMatrix;
}

////////////// CREATE DISPLACEMENTS DOF HANDLES //////////////

il::Array2D<il::int_t> createAutoDisplacementDofHandle(const il::int_t interpolationOrder,
                                                       const il::int_t numberOfElements) {

  il::int_t columns_displ_dofh_mtx;

  if(interpolationOrder == 0){

    // 2 is the problem dimension (2D) which means 2 displacements at each node, even if constant element
    // Also, for the constant element, we consider 2 nodes per element.
    columns_displ_dofh_mtx = 2*2;

  } else {

    // In general, at each node we have 2 displacements (x,y). Every element possesses interpolationOrder + 1 nodes.
    columns_displ_dofh_mtx = 2*(interpolationOrder + 1);

  }

  il::Array2D<il::int_t> dof_handle_displacement(numberOfElements, columns_displ_dofh_mtx);

// filling the dof_handle_displacements (2D discontinuous Galerkin)
#pragma omp parallel for collapse(2)
  for (il::int_t i = 0, j; i < numberOfElements; i++) {

    for (int k = 0; k < columns_displ_dofh_mtx; k++) {

      dof_handle_displacement(i, k) =  i * columns_displ_dofh_mtx + k;

    }

  }
  return dof_handle_displacement;
}

////////////// CREATE PRESSURE DOF HANDLES //////////////

il::Array2D<il::int_t> createAutoPressureDofHandle(const il::int_t interpolationOrder,
                                                   const il::int_t numberOfElements) {

  il::int_t columns_press_dofh_mtx;

  if(interpolationOrder == 0){

    // 1 pressure dof at each node. In the constant element, the nodes are 2
    columns_press_dofh_mtx = 2;

  } else {

    // In general, every element possesses interpolationOrder+1 nodes with 1 pressure value per node
    columns_press_dofh_mtx = (interpolationOrder + 1); // 2D problem with generic interpolation order
  }

  il::Array2D<il::int_t> dof_handle_pressure(numberOfElements, columns_press_dofh_mtx);

// fillind the dof_handle_pressure (continuous Galerkin)
#pragma omp parallel for collapse(2)
  for (il::int_t i = 0; i < numberOfElements; i++) {

    for (il::int_t j = 0; j < columns_press_dofh_mtx; j++) {

      dof_handle_pressure(i, j) = i * (columns_press_dofh_mtx-1) + j;

    }

  }

  return dof_handle_pressure;
}

}