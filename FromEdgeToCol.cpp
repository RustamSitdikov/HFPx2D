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
#include "DOF_Handles.h"
#include "FromEdgeToCol.h"

namespace hfp2d {

// Function which allow us to switch from end points (two values per node ->
// dof_dim = 2) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference unit element at location
// {-1/sqrt(2) , +1/sqrt(2)}
// It returns an array (vector) that includes EITHER shear (slip) OR opening at
// collocation points

il::Array<double> from_edge_to_col(il::Array<double> &d_edge, const int Nelts,
                                const int dof_dim) {

  // Inputs:
  //  - d_edge -> vector that contains EITHER slip OR opening at nodal points
  //  - Nelts -> number of elements
  //  - dof_dim -> degrees of freedom per nodes

  il::Array<double> Fetc{d_edge.size(), 0.};
  il::Array2D<double> ShapeFunction{2, 2, .0};

  ShapeFunction(0, 0) = (1 + (1 / sqrt(2))) / 2;
  ShapeFunction(0, 1) = (1 - (1 / sqrt(2))) / 2;

  ShapeFunction(1, 0) = (1 - (1 / sqrt(2))) / 2;
  ShapeFunction(1, 1) = (1 + (1 / sqrt(2))) / 2;

  il::Array2D<int> A{Nelts, 2, 0};

  for (il::int_t k = 0, j; k < A.size(0); ++k) {

    j = k * dof_dim;

    for (il::int_t i = 0; i < A.size(1); ++i) {

      A(k, i) = i + j;
    }
  }

  for (il::int_t i = 0; i < Nelts; ++i) {

    for (int j = 0; j < dof_dim; j = j + 2) {

      Fetc[A(i, 0)] = Fetc[A(i, 0)] + ShapeFunction(0, j) * d_edge[A(i, j)];
      Fetc[A(i, 1)] = Fetc[A(i, 1)] + ShapeFunction(1, j) * d_edge[A(i, j)];
    }
  }

  return Fetc;
};
}