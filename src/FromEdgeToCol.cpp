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
#include "src/DOF_Handles.h"
#include "FromEdgeToCol.h"

namespace hfp2d {

// Function which allow us to switch from end points (two values per node ->
// dof_dim = 2) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns a matrix (size 2Nelts x 2Nelts) that, if multiplied by nodal
// slip OR opening vector (size 2Nelts), it returns the slip OR opening vector
// at collocation points

// this only works for p=1
il::Array2D<double> from_edge_to_col(const int nelts, const int dof_dim) {

  // Inputs:
  //  - Nelts -> number of elements
  //  - dof_dim -> degrees of freedom per nodes
  //
  // it should be a sparse matrix....
  il::Array2D<double> Fetc{2 * nelts, 2 * nelts, 0.};
  il::StaticArray2D<double> ShapeFunction{2, 2, .0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  il::Array2D<int> A{nelts, 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * dof_dim;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  for (il::int_t i = 0, k = 0, q = 1; i < nelts; ++i) {

    for (il::int_t j = 0; j < dof_dim; ++j) {

        Fetc(k, A(i, j)) = ShapeFunction(0, j);
        Fetc(q, A(i, j)) = ShapeFunction(1, j);

    }

    k = k + 2;
    q = q + 2;

  }

  return Fetc;
};
}
