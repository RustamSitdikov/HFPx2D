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

// Function which allow us to switch from end points (two values per node ->
// discontinuous polinomial) to collocation points
// Dof.size(0) = number of elements
// 4 * Dof.size(0) = number of collocation points for DD

// It returns a matrix (size 4Nelts x 4Nelts) that, if multiplied by nodal
// slip + opening vector (size 4Nelts), it returns the slip AND opening vector
// at collocation points

il::Array2D<double> from_edge_to_col_dg_full2d_new(Mesh &mesh) {

  // Inputs:
  //  - dof_dim -> degrees of freedom per nodes
  //  - Dofw -> DOFs handle for BOTH slip AND opening (size -> Neltsx4)

  // Note matrix on all the DDs dofs
  const il::int_t DDxElem = mesh.numDDDofsPerElem();
  il::Array2D<double> Fetc{mesh.numDDDofs(), mesh.numDDDofs(), 0};
  il::StaticArray2D<double,4,4> ShapeFunction{0};

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

  for (il::int_t i = 0; i < mesh.numElems(); ++i)
  {

    for (il::int_t j = 0; j < DDxElem; j++)
    {
      for (il::int_t k = 0; k < DDxElem; k++)
      {
        Fetc(mesh.dofDD(i, j), mesh.dofDD(i, k)) = ShapeFunction(j, k);
      }
    }
  }
    // TODO: avoid copying zeros and use only the two main values of ShapeFunc.

//  il::StaticArray2D<double,2,4> ShapeFunction{0};
//
//  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
//  ShapeFunction(0, 1) = (1 + (1 / sqrt(2))) / 2;
//  ShapeFunction(0, 2) = (1 - (1 / sqrt(2))) / 2;
//  ShapeFunction(0, 3) = (1 - (1 / sqrt(2))) / 2;
//
//  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
//  ShapeFunction(1, 1) = (1 - (1 / sqrt(2))) / 2;
//  ShapeFunction(1, 2) = (1 + (1 / sqrt(2))) / 2;
//  ShapeFunction(1, 3) = (1 + (1 / sqrt(2))) / 2;
//
//  for (il::int_t i = 0; i < mesh.numElems(); ++i) {
//
//    for (il::int_t j = 0; j < mesh.numDDDofsPerElem(); ++j) {
//
//      Fetc(mesh.dofDispl(i, 0), mesh.dofDD(i, j)) = ShapeFunction(0, j);
//      Fetc(mesh.dofDispl(i, 1), mesh.dofDD(i, j)) = ShapeFunction(0, j);
//      Fetc(mesh.dofDispl(i, 2), mesh.dofDD(i, j)) = ShapeFunction(1, j);
//      Fetc(mesh.dofDispl(i, 3), mesh.dofDD(i, j)) = ShapeFunction(1, j);
//    }
//  }


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
  // Dof(i,j) = mesh.dofDD(i,j)
  // Dofp(i,j) = mesh.dofPress(i,j)
  // dof_dim = 2 (number of dimensions of the problem)
  // Dof.size(0) = mesh.numNodes()
  // 4 = mesh.numDDDofsPerElem()

  // Internal pressure contributes only to the opening degrees of freedom
  // so components 1 and 3

  // Note matrix on all the DDs dofs
  const il::int_t dim = 2; // 2 because the problem is 2D
  il::Array2D<double> Fetc{theMesh.numDDDofs(), theMesh.numPressDofs(), 0};
  il::StaticArray2D<double,2,2> ShapeFunction{0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  for (il::int_t i = 0; i < theMesh.numElems(); ++i) {

    for (il::int_t j = 0; j < dim; ++j) {

      Fetc(theMesh.dofDD(i, 1), theMesh.dofPress(i, j)) = ShapeFunction(0, j);
      Fetc(theMesh.dofDD(i, 3), theMesh.dofPress(i, j)) = ShapeFunction(1, j);
    }
  }

  return Fetc;
}
}