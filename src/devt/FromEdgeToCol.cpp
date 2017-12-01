//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from the standard library
#include <cmath>
#include <iostream>

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

// Inclusion from the project
#include "FromEdgeToCol.h"

namespace hfp2d {

// Function which allow us to switch from end points (two values per node ->
// discontinuous polinomial) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns a matrix (size 4Nelts x 4Nelts) that, if multiplied by nodal
// slip + opening vector (size 4Nelts), it returns the slip AND opening vector
// at collocation points

// Function which allow us to switch from end points (two values per node ->
// discontinuous polinomial) to collocation points
// Dof.size(0) = number of elements
// 4 * Dof.size(0) = number of collocation points for DD

// It returns a matrix (size 4Nelts x 4Nelts) that, if multiplied by nodal
// slip + opening vector (size 4Nelts), it returns the slip AND opening vector
// at collocation points

il::Array2D<double> from_edge_to_col_dg_full2d(Mesh &theMesh) {
  // input:
  //  - theMesh -> object of type Mesh

  il::Array2D<double> Fetc{theMesh.numberDDDofs(), theMesh.numberDDDofs(), 0};
  il::StaticArray2D<double, 4, 4> ShapeFunction{0};

  // contribution to slip at collocation node 1
  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 2) = (1 - (1 / sqrt(2))) / 2;

  // contribution to opening at collocation node 1
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(1, 3) = (1 - (1 / sqrt(2))) / 2;

  // contribution to slip at collocation node 2
  ShapeFunction(2, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(2, 2) = (1 + (1 / sqrt(2))) / 2;

  // contribution to opening at collocation node 2
  ShapeFunction(3, 1) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(3, 3) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < theMesh.numberOfElts(); ++i) {
    for (il::int_t j = 0; j < theMesh.numberDDDofsPerElt(); j++) {
      for (il::int_t k = 0; k < theMesh.numberDDDofsPerElt(); k++) {
        Fetc(theMesh.dofDD(i, j), theMesh.dofDD(i, k)) = ShapeFunction(j, k);
      }
    }
  }

  return Fetc;
};

// TODO store the matrix as sparse
// Function which allow us to switch from end points (two values per node ->
// discontinuous polinomial) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns a matrix (size 2Nelts x 2Nelts) that, if multiplied by EITHER
// slip OR opening vector (size 2*Nelts), it returns the slip OR opening vector
// at collocation points
il::Array2D<double> from_edge_to_col_dg(Mesh &theMesh) {
  // input:
  //  - theMesh -> object of type Mesh

  il::Array2D<il::int_t> dof_single_dd{theMesh.numberOfElts(),
                                       (theMesh.interpolationOrder() + 1), 0};
  for (il::int_t i = 0; i < theMesh.numberOfElts(); i++) {
    for (il::int_t j = 0; j < 1 * (theMesh.interpolationOrder() + 1); j++) {
      dof_single_dd(i, j) = i * 1 * (theMesh.interpolationOrder() + 1) + j;
    }
  }

  // Note matrix on all the DDs dofs
  il::Array2D<double> Fetc{2 * dof_single_dd.size(0), 2 * dof_single_dd.size(0),
                           0};
  il::StaticArray2D<double, 2, 2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < dof_single_dd.size(0); ++i) {
    for (il::int_t j = 0; j < theMesh.numberDDDofsPerElt() / 2; ++j) {
      Fetc(dof_single_dd(i, 0), dof_single_dd(i, j)) = ShapeFunction(0, j);
      Fetc(dof_single_dd(i, 1), dof_single_dd(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
};

// Function which allow us to switch from end points (one values per node
// -> continuous polinomial) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns a matrix (size 4Nelts x Nelts+1) that, if multiplied by nodal
// pressure values (size Nnodes), it returns the pressure values at
// collocation points
/// Note that all the dofs are considered -> not optimal
il::Array2D<double> from_edge_to_col_cg(Mesh &theMesh) {
  il::Array2D<int> dofhandle_press{theMesh.numberOfElts(),
                                   theMesh.interpolationOrder() + 1, 0};
  for (int i = 0; i < dofhandle_press.size(0); ++i) {
    dofhandle_press(i, 0) = i;
    dofhandle_press(i, 1) = i + 1;
  }

  il::Array2D<il::int_t> dof_single_dd{theMesh.numberOfElts(),
                                       (theMesh.interpolationOrder() + 1), 0};
  for (il::int_t i = 0; i < theMesh.numberOfElts(); i++) {
    for (il::int_t j = 0; j < 1 * (theMesh.interpolationOrder() + 1); j++) {
      dof_single_dd(i, j) = i * 1 * (theMesh.interpolationOrder() + 1) + j;
    }
  }

  il::Array2D<double> Fetc{theMesh.numberDDDofs()/2, theMesh.numberOfNodes(),
                           0};
  il::StaticArray2D<double, 2, 2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < theMesh.numberOfElts(); ++i) {
    for (il::int_t j = 0; j < (theMesh.interpolationOrder() + 1); ++j) {
      Fetc(dof_single_dd(i, 0), dofhandle_press(i, j)) = ShapeFunction(0, j);
      Fetc(dof_single_dd(i, 1), dofhandle_press(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
}
}