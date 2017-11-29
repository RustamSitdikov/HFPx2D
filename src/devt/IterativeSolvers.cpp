//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 29.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "src/devt/IterativeSolvers.h"

namespace hfp2d {

il::Array<double> jacobiIterativeSolver(il::Array2D<double> &matrix_of_coeff,
                                        il::Array<double> &right_hand_side,
                                        il::Array<double> &initial_solution,
                                        Mesh &theMesh, double relax_parameter,
                                        il::int_t Niterations) {
  il::Array<double> BigX{initial_solution.size(), 0.};

  IL_EXPECT_FAST(matrix_of_coeff.size(0) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(1) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(0) == right_hand_side.size());
  IL_EXPECT_FAST(right_hand_side.size() == initial_solution.size());

  double a_i, b_i;

  while (Niterations > 0) {
    for (il::int_t i = 0; i < matrix_of_coeff.size(0); ++i) {
      a_i = matrix_of_coeff(i, i);
      b_i = right_hand_side[i];
      BigX[i] = relax_parameter * (b_i / a_i) +
                (1 - relax_parameter) * initial_solution[i];

      for (il::int_t j = 0; j < matrix_of_coeff.size(1); ++j) {
        if (j == i) continue;

        BigX[i] = BigX[i] -
                  relax_parameter *
                      (matrix_of_coeff(i, j) / matrix_of_coeff(i, i)) *
                      initial_solution[j];
        initial_solution[i] = BigX[i];
      }
    }

    Niterations--;
  }

  return BigX;
};

il::Array<double> gaussSeidelIterativeSolver(
    il::Array2D<double> &matrix_of_coeff, il::Array<double> &right_hand_side,
    il::Array<double> &initial_solution, Mesh &theMesh, double relax_parameter,
    il::int_t Niterations) {
  il::Array<double> BigX{initial_solution.size(), 0.};

  IL_EXPECT_FAST(matrix_of_coeff.size(0) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(1) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(0) == right_hand_side.size());
  IL_EXPECT_FAST(right_hand_side.size() == initial_solution.size());

  double a_i, b_i, sum;

  while (Niterations > 0) {
    for (il::int_t i = 0; i < matrix_of_coeff.size(0); ++i) {
      a_i = matrix_of_coeff(i, i);
      b_i = right_hand_side[i];
      sum = 0;

      for (il::int_t j = 0; j < i; ++j) {
        sum += (matrix_of_coeff(i, j) * BigX[j]);
      }

      for (il::int_t k = i + 1; k < matrix_of_coeff.size(1); ++k) {
        sum += (matrix_of_coeff(i, k)) * initial_solution[k];
      }

      BigX[i] = (relax_parameter / a_i) * (b_i - sum) +
                (1 - relax_parameter) * initial_solution[i];
      initial_solution[i] = BigX[i];
    }

    Niterations--;
  }

  return BigX;
};
}