//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 03.03.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_ELHDS_H
#define HFPX2D_ELHDS_H

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linear_algebra/dense/norm.h>

// Inclusion from the project
#include "Dilatancy.h"
#include "FVM.h"
#include "Friction.h"
#include "MC_criterion.h"
#include "Mesh.h"

namespace hfp2d {

struct Results_solution_nonlinearsystem {

  // Error on increment in pressure
  double errDp;

  // Error on increment in slip
  double errDd;

  // New matrix of total slip
  il::Array2D<double> new_tot_slip;

  // Final increment of slip after iterations
  il::Array<double> Dd;

  // Final increment of pressure after iterations
  il::Array<double> Dp;

  // Number of iterations to solve non linear system
  int Niterations;
};

Results_solution_nonlinearsystem
ELHDs(Mesh mesh, il::Array2D<double> &kmatd, il::Array2D<double> &Npc,
      Parameters_friction fric_parameters,
      Parameters_dilatancy dilat_parameters, Parameters_fluid fluid_parameters,
      Results_one_timeincrement &SolutionAtTj, il::Array2D<double> sigma_tot,
      il::Array<double> press_prof, il::Array2D<double> Pcm,
      il::Array<double> tot_slip, int itermax, int dof_dim, int p,
      il::Array<double> cohes, il::Status status, il::Norm norm, int inj_point,
      il::Array<double> S, il::Array<int> Dof_slip_coll, il::io_t);
}

#endif // HFPX2D_ELHDS_H
