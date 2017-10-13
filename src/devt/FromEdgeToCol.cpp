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

il::Array2D<double> from_edge_to_col_dg_full2d(int dof_dim,
                                               il::Array2D<int> Dof, il::io_t) {

  // Inputs:
  //  - dof_dim -> degrees of freedom per nodes
  //  - Dofw -> DOFs handle for BOTH slip AND opening (size -> Neltsx4)

  // Note matrix on all the DDs dofs
  il::Array2D<double> Fetc{4 * Dof.size(0), 4 * Dof.size(0), 0};
  il::StaticArray2D<double,2,4> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 2) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(0, 3) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 2) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(1, 3) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < Dof.size(0); ++i) {

    for (il::int_t j = 0; j < 2 * dof_dim; ++j) {

      Fetc(Dof(i, 0), Dof(i, j)) = ShapeFunction(0, j);
      Fetc(Dof(i, 1), Dof(i, j)) = ShapeFunction(0, j);
      Fetc(Dof(i, 2), Dof(i, j)) = ShapeFunction(1, j);
      Fetc(Dof(i, 3), Dof(i, j)) = ShapeFunction(1, j);
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
il::Array2D<double> from_edge_to_col_dg(int dof_dim,
                                        il::Array2D<int> Dofw, il::io_t) {

  // Inputs:
  //  - dof_dim -> degrees of freedom per nodes
  //  - Dofw -> DOFs handle for EITHER slip OR opening (size -> Neltsx2)

  // Note matrix on all the DDs dofs
  il::Array2D<double> Fetc{2 * Dofw.size(0), 2 * Dofw.size(0), 0};
  il::StaticArray2D<double,2,2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < Dofw.size(0); ++i) {

    for (il::int_t j = 0; j < dof_dim; ++j) {

      Fetc(Dofw(i, 0), Dofw(i, j)) = ShapeFunction(0, j);
      Fetc(Dofw(i, 1), Dofw(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
};

// Function which allow us to switch from end points (one values per node
// -> continuous polinomial)
// to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns a matrix (size 4Nelts x Nelts+1) that, if multiplied by nodal
// pressure values (size Nelts + 1), it returns the pressure values at
// collocation points
il::Array2D<double> from_edge_to_col_cg(int dof_dim,
                                        il::Array2D<int> Dof,
                                        il::Array2D<int> Dofp, il::io_t) {

  // Note matrix on all the DDs dofs
  il::Array2D<double> Fetc{4 * Dof.size(0), Dof.size(0) + 1, 0};
  il::StaticArray2D<double,2,2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < Dof.size(0); ++i) {

    for (il::int_t j = 0; j < dof_dim; ++j) {

      Fetc(Dof(i, 1), Dofp(i, j)) = ShapeFunction(0, j);
      Fetc(Dof(i, 3), Dofp(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
}


il::Array2D<double> from_edge_to_col_cg_new(Mesh &theMesh) {

  // From mesh
  // Dof(i,j) = mesh.dofDispl(i,j)
  // Dofp(i,j) = mesh.dofPress(i,j)
  // dof_dim = 2 (number of dimensions of the problem)
  // Dof.size(0) = mesh.numNodes()
  // 4 = mesh.numDisplDofsPerElem()

  // Internal pressure contributes only to the opening degrees of freedom
  // so components 1 and 3

  // Note matrix on all the DDs dofs
  const il::int_t dim = 2; // 2 because the problem is 2D
  il::Array2D<double> Fetc{theMesh.numDisplDofs(), theMesh.numPressDofs(), 0};
  il::StaticArray2D<double,2,2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < theMesh.numElems(); ++i) {

    for (il::int_t j = 0; j < dim; ++j) {

      Fetc(theMesh.dofDispl(i, 1), theMesh.dofPress(i, j)) = ShapeFunction(0, j);
      Fetc(theMesh.dofDispl(i, 3), theMesh.dofPress(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
}
}