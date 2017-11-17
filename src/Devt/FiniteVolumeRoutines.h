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
// Input: vector of slip/opening for each element -> size 2 x Nelts
// Remember: piecewise linear variation over the element
// Output: vector that contains the average values for each element (size Nelts)
inline il::Array<double> average(const il::Array<double> &d) {
  il::Array<double> Average{d.size() / 2, 0};
  for (il::int_t i = 0, k = 0; i < Average.size(); ++i, k = k + 2) {
    Average[i] = (d[k] + d[k + 1]) / 2;
  }
  return Average;
};

///
// This function calculates the slip/opening at +/- 1/4 -> the control volume is
// centered on the nodes!
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector -> {slip_+1/4 , slip_+3/4}

inline il::Array<double> quarter(const il::Array<double> &d) {
  il::Array<double> Quarter(d.size(), 0);
  for (il::int_t i = 0; i < d.size() / 2; ++i, i = i + 2) {
    Quarter[i] = ((3 * d[i]) + d[i + 1]) / 4;
    Quarter[i + 1] = (d[i] + (3 * d[i + 1])) / 4;
  }
  return Quarter;
};

///
// Function to find out the position of a value in a 2D array
// It returns 2x2 array with row&col of the seek value
// It is completely general in a sense that the output can be a vector or a
// matrix (2x2)
inline il::Array2D<il::int_t> position_2d_array(
    const il::Array2D<il::int_t> &arr2D, int seek) {
  il::Array2D<il::int_t> M{arr2D.size(1) * arr2D.size(0), 2, -1};
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

  il::Array2D<il::int_t> outp{k, 2, 0};

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
inline il::Array<il::int_t> row_selection(const il::Array2D<il::int_t> &arr,
                                          il::int_t idx) {
  il::Array<il::int_t> vect{arr.size(1), 0};
  for (il::int_t i = 0; i < vect.size(); ++i) {
    vect[i] = arr(idx, i);
  }
  return vect;
};

il::Array<double> edgeConductivitiesP1Newtonian(
    Mesh &theMesh, FluidProperties &FluidProperties,
    il::Array<double> &permeab_middle, const il::Array<double> &dilat_middle,
    const il::Array<double> &opening_middle);

il::Array<double> shearConductivitiesP1Newtonian(
    Mesh &theMesh, FluidProperties &FluidProperties,
    FractureEvolution &FractureEvolution, const il::Array<double> &slip);

il::Array2D<double> buildLMatrix(Mesh &theMesh, const il::Array<double> &slip,
                                 const il::Array<double> &opening,
                                 FluidProperties &FluidProperties,
                                 FractureEvolution &FractureEvolution,
                                 const double TimeStep);

il::Array2D<double> buildVpMatrix(Mesh &theMesh,
                                  FractureEvolution &FractureEvolution,
                                  FluidProperties &FluidProperties,
                                  il::Array<double> &slip);

il::Array2D<double> buildVdMatrix(Mesh &theMesh,
                                  FractureEvolution &FractureEvolution,
                                  FluidProperties &FluidProperties,
                                  il::Array<double> &slip);
}

#endif  // HFPX2D_FVM_H
