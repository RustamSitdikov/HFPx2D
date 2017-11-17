//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FVM_H
#define HFPX2D_FVM_H

// Inclusion from standard library
#include <cmath>

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

// Inclusion from the project
#include "Dilatancy.h"
#include "Permeability.h"
#include "src/Core/Mesh.h"

namespace hfp2d {

struct Parameters_fluid {
  // Fluid viscosity
  double viscosity;
  // Fluid compressibility
  double compressibility;
  // Fluid density matrix
  // {{rho_1left, rho_2right},{rho_2left, rho_2right}, ...} (size -> Nelts + 1)
  il::Array2D<double> density;
};

///
// This function calculates the average between two values for each element
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector that contains the average values of each row (element)
inline il::Array<double> average(const il::Array2D<double> &d) {
  il::Array<double> Average{d.size(0), 0};
  for (il::int_t i = 0; i < d.size(0); ++i) {
    Average[i] = (d(i, 0) + d(i, 1)) / 2;
  }
  return Average;
};
///
// This function calculates the slip/opening at +/- 1/4 -> the control volume is
// centered on the nodes!
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector -> {slip_+1/4 , slip_+3/4}

inline il::Array<double> quarter(const il::Array2D<double> &d) {
  il::Array<double> Quarter(2 * d.size(0), 0);
  for (il::int_t i = 0, j = 0; i < (d.size(0)); ++i, j = j + 2) {
    Quarter[j] = ((3 * d(i, 0)) + d(i, 1)) / 4;
    Quarter[j + 1] = (d(i, 0) + (3 * d(i, 1))) / 4;
  }
  return Quarter;
};

///
// Function to find out the position of a value in a 2D array
// It returns 2x2 array with row&col of the seek value
// It is completely general in a sense that the output can be a vector or a
// matrix (2x2)
inline il::Array2D<int> position_2d_array(const il::Array2D<int> &arr2D, int seek) {
  il::Array2D<int> M{arr2D.size(1) * arr2D.size(0), 2, -1};
  int k = 0;

  for (int i = 0; i < arr2D.size(0); ++i) {
    for (int j = 0; j < arr2D.size(1); ++j) {
      if (arr2D(i, j) == seek) {
        M(k, 0) = i;
        M(k, 1) = j;

        k = k + 1;
      }
    }
  }

  il::Array2D<int> outp{k, 2, 0};

  for (int l = 0; l < k; ++l) {
    for (int j = 0; j < 2; ++j) {
      outp(l, j) = M(l, j);
    }
  }

  return outp;
};

///
// Auxiliary function for assembly process
// It returns a given row (vector - specified by idx) of a 2D array
inline il::Array<int> row_selection(const il::Array2D<int> &arr,
                                    il::int_t idx) {
  il::Array<int> vect{arr.size(1), 0};
  for (il::int_t i = 0; i < vect.size(); ++i) {
    vect[i] = arr(idx, i);
  }

  return vect;
};

il::Array<double> shear_conductivities_newtonian(
    Parameters_fluid &fluid_parameters, Mesh mesh, const il::Array2D<double> &d,
    Parameters_dilatancy &dilat_parameters,
    Parameters_permeability &permeab_parameters);

il::Array2D<double> build_l_matrix(Mesh mesh, const il::Array2D<double> &d,
                                   Parameters_fluid &fluid_parameters,
                                   Parameters_dilatancy &dilat_parameters,
                                   const double &TimeStep,
                                   Parameters_permeability &permeab_parameters);

il::Array2D<double> build_vp_matrix_p1(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array2D<double> &d, il::io_t);

il::Array2D<double> build_vd_matrix_p1(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       il::Array2D<int> Dof,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array2D<double> &d);
}

#endif  // HFPX2D_FVM_H
