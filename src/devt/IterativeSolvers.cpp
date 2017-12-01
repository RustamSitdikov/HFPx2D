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
                                        double relax_parameter,
                                        il::int_t numb_iterations) {

  // This function solve a linear system of equations Ax=b through Jacobi
  // iterative
  // method, combined with an Over Relaxation Method (JOR).
  // Reference: "Numerical Mathematics" by Quarteroni A., Sacco R. and Saleri F.
  // Inputs:
  //        - matrix_of_coeff -> A matrix
  //        - right_hand_side -> b vector
  //        - initial solution -> x_o
  //        - relax_parameter -> relaxation parameter
  //        - num_iterations -> total number of iterations
  // Output:
  //        - solution -> x vector

  il::Array<double> solution{initial_solution.size(), 0.};

  IL_EXPECT_FAST(matrix_of_coeff.size(0) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(1) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(0) == right_hand_side.size());
  IL_EXPECT_FAST(right_hand_side.size() == initial_solution.size());

  while (numb_iterations > 0) {
    for (il::int_t i = 0; i < matrix_of_coeff.size(0); ++i) {
      solution[i] =
          relax_parameter * (right_hand_side[i] / matrix_of_coeff(i, i)) +
          (1 - relax_parameter) * initial_solution[i];

      for (il::int_t j = 0; j < matrix_of_coeff.size(1); ++j) {
        if (j == i) continue;

        solution[i] = solution[i] -
                      relax_parameter *
                          (matrix_of_coeff(i, j) / matrix_of_coeff(i, i)) *
                          initial_solution[j];
        initial_solution[i] = solution[i];
      }
    }

    numb_iterations--;
  }
  return solution;
};

il::Array<double> gaussSeidelIterativeSolver(
    il::Array2D<double> &matrix_of_coeff, il::Array<double> &right_hand_side,
    il::Array<double> &initial_solution, double relax_parameter,
    il::int_t numb_iterations) {

  // This function solve a linear system of equations Ax=b through Gauss-Seidel
  // iterative method, combined with an Over Relaxation Method (JOR).
  // Reference: "Numerical Mathematics" by Quarteroni A., Sacco R. and Saleri F.
  // Inputs:
  //        - matrix_of_coeff -> A matrix
  //        - right_hand_side -> b vector
  //        - initial solution -> x_o
  //        - relax_parameter -> relaxation parameter
  //        - num_iterations -> total number of iterations
  // Output:
  //        - solution -> x vector

  il::Array<double> solution{initial_solution.size(), 0.};

  IL_EXPECT_FAST(matrix_of_coeff.size(0) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(1) == initial_solution.size());
  IL_EXPECT_FAST(matrix_of_coeff.size(0) == right_hand_side.size());
  IL_EXPECT_FAST(right_hand_side.size() == initial_solution.size());

  double sum;

  while (numb_iterations > 0) {
    for (il::int_t i = 0; i < matrix_of_coeff.size(0); ++i) {
      sum = 0;

      for (il::int_t j = 0; j < i; ++j) {
        sum += (matrix_of_coeff(i, j) * solution[j]);
      }

      for (il::int_t k = i + 1; k < matrix_of_coeff.size(1); ++k) {
        sum += (matrix_of_coeff(i, k)) * initial_solution[k];
      }

      solution[i] = (relax_parameter / matrix_of_coeff(i, i)) *
                        (right_hand_side[i] - sum) +
                    (1 - relax_parameter) * initial_solution[i];
      initial_solution[i] = solution[i];
    }

    numb_iterations--;
  }
  return solution;
};
}