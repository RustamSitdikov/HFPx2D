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
#include "src/core/Mesh.h"

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

il::Array<double> average(const il::Array2D<double> &d);

il::Array<double> quarter(const il::Array2D<double> &d);

il::Array2D<int> position_2d_array(const il::Array2D<int> &arr2D, int seek);

il::Array2D<int> search(const il::Array2D<int> &matrix, int x);

il::Array<int> row_selection(const il::Array2D<int> &arr, il::int_t idx);

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
