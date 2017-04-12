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
EHLDs(Mesh mesh, il::Array2D<double> &kmatd, il::Array2D<double> &Npc,
      LayerParameters1 &layer_parameters1, LayerParameters2 &layer_parameters2,
      LayerParameters3 &layer_parameters3, il::Array<il::int_t> id_layers,
      Parameters_dilatancy dilat_parameters, Parameters_fluid fluid_parameters,
      Results_one_timeincrement &SolutionAtTj, il::Array<double> press_prof,
      il::Array<double> tot_slip, int dof_dim, int p, il::Array<double> cohes,
      il::Status &status, il::Norm norm, int inj_point, il::Array<double> S,
      il::Array<il::int_t> Dof_slip_coll, il::Array2D<double> Sigma0,
      il::Array2D<double> sigma_tot,
      hfp2d::simulation_parameters simulation_parameters, double kf, il::io_t);
}

#endif // HFPX2D_EHLDS_H
