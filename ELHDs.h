//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 30.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_ELHDS_H
#define HFPX2D_ELHDS_H

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include "Mesh.h"

namespace hfp2d {

struct Result {

  // vector of friction coefficient for each time step
  il::Array<double> friction;

  // vector of dilatancy for each time step
  il::Array<double> dilatancy;

  // array (matrix) of total stress for each time step
  il::Array2D<double> tot_stress_state;

  // vector of increment of slip
  il::Array<double> incr_d;

  // vector of pore pressure profile at nodal points (size -> 2Nelts)
  il::Array<double> P;

  // number of iterations for each time step
  int iter;

  // vector of final slip points per each time step
  il::Array<il::int_t> final_slip_elements;

  // current time step
  double dt;

  // Length slippage zone
  double slippagezone;
};

void elhds(Result &res, Mesh mesh, const int p,
                    const double Cohes, il::Array2D<double> &kmat,
                    const double Incr_dil, const double d_wd,
                    il::Array2D<double> rho, const double Init_dil,
                    const double CompressFluid, const double TimeStep,
                    const double Visc, il::Array<double> S, const int InjPoint,
                    const int dof_dim, const double Peak_fric,
                    const double Resid_fric, const double d_wf);

il::Array<double> flatten1(il::Array2D<double> &Arr);

int boole_mc(il::Array<double> tau, il::Array<double> sigma_n,
             il::Array<double> fric, const double Cohes);

il::Array<int> flatten2(il::Array2D<int> Arr);

il::Array<int> select(il::Array<int> &arr);

void set_submatrix_non_linear_system(il::Array2D<double> &A, int i0, int i1,
                                     const il::Array2D<double> &B);

void set_subvector_non_linear_system(il::Array<double> &A, int i0,
                                     const il::Array<double> &B);

il::Array2D<double>
take_submatrix_non_linear_system(int i0, int i1, int j0, int j1,
                                 const il::Array2D<double> &A);

il::Array<double> take_subvector_non_linear_system(int i0, int i1,
                                                   const il::Array<double> &A);

double euclidean_distance(double x1, double x2, double y1, double y2);

il::Array<int> delete_duplicates(il::Array<il::int_t> &arr);

il::Array<il::int_t> delete_duplicates2(const il::Array<il::int_t> &arr);
}

#endif // HFPX2D_ELHDS_H
