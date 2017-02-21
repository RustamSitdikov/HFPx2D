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

void elhds(Mesh mesh, int p, double Cohes, const il::Array2D<double> &kmat,
           double Incr_dil, double d_wd, il::Array2D<double> rho,
           double Init_dil, double CompressFluid, double TimeStep, double Visc,
           il::Array<double> S, int InjPoint, int dof_dim, double Peak_fric,
           double Resid_fric, double d_wf, il::Array<double> XColl, il::io_t,
           Result &res);

il::Array<double> flatten1(const il::Array2D<double> &Arr, il::io_t);

int boole_mc(il::Array<double> tau, il::Array<double> sigma_n,
             il::Array<double> fric, double Cohes, il::io_t);

il::Array<int> flatten2(const il::Array2D<int> &Arr, il::io_t);

il::Array<int> select(const il::Array<int> &arr, il::io_t);

void set_submatrix_non_linear_system(il::Array2D<double> &A, int i0, int i1,
                                     const il::Array2D<double> &B);

void set_subvector_non_linear_system(il::Array<double> &A, int i0,
                                     const il::Array<double> &B);

il::Array2D<double>
take_submatrix_non_linear_system(int i0, int i1, int j0, int j1,
                                 const il::Array2D<double> &A);

il::Array<double> take_subvector_non_linear_system(int i0, int i1,
                                                   const il::Array<double> &A);

double euclidean_distance(double x1, double x2, double y1, double y2, il::io_t);

il::Array<int> delete_duplicates(const il::Array<il::int_t> &arr, il::io_t);

il::Array<il::int_t> delete_duplicates2(const il::Array<il::int_t> &arr,
                                        il::io_t);

void sort_ascending_order(const il::Array<il::int_t> &arr, il::io_t);
}

#endif // HFPX2D_ELHDS_H
