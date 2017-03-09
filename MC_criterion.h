//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 02.03.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_MC_CRITERION_H
#define HFPX2D_MC_CRITERION_H

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include "Mesh.h"

namespace hfp2d {

struct Results_one_timeincrement {

  // vector of friction coefficient for each time step
  il::Array<double> friction;

  // array (matrix) of total stress for each time step
  il::Array2D<double> tot_stress_state;

  // vector of total slip
  il::Array<double> d_tot;

  // vector of pore pressure profile at nodal points (size -> 2Nelts)
  il::Array<double> P;

  // vector of pore pressure profile at collocation points (size -> 2Nelts x 2)
  il::Array2D<double> Pcm;

  // number of iterations for each time step
  int iter;

  // vector of active set of collocation points
  il::Array<il::int_t> active_set_collpoints{};

  // current time step
  double dt;

  // Length slippage zone
  double slippagezone;

  // Shear crack velocity
  double crack_velocity;
};

struct MCcheck {

  // Number of collocation points that satisfy the Mohr-Coulomb criterion
  int Ncollpoint_satisfMC = 0;

  // Collocation points that do not satisfy the Mohr-Coulomb criterion
  il::Array<il::int_t> CollPoint_notsatisfMC{};
};

void MC_criterion(Mesh mesh, int p, il::Array<double> cohes,
                  const il::Array2D<double> &kmat,
                  Parameters_friction &fric_parameters,
                  Parameters_dilatancy &dilat_parameters,
                  Parameters_fluid &fluid_parameters, il::Array<double> S,
                  int inj_point, int dof_dim, il::Array<double> XColl,
                  il::Array2D<double> &Fetc, il::Array2D<double> Sigma0,
                  il::io_t, Results_one_timeincrement &res);

il::Array<double> flatten1(const il::Array2D<double> &Arr, il::io_t);

MCcheck boole_mc(il::Array2D<double> &sigma_eff, il::Array<double> fric,
                 il::Array<double> cohes, il::io_t);

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

il::Array<int> delete_duplicates(const il::Array<il::int_t> &arr, il::io_t);

il::Array<il::int_t> delete_duplicates2(const il::Array<il::int_t> &arr,
                                        il::io_t);

void sort_ascending_order(il::Array<il::int_t> &arr, il::io_t);
}

#endif // HFPX2D_MC_CRITERION_H
