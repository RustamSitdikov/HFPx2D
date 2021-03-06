//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 07.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_SIMULATIONPARAMETERS_H
#define HFPX2D_SIMULATIONPARAMETERS_H


#include <il/math.h>

namespace hfp2d {


struct SimulationParameters {
  // placeholder for simulation parameters

  il::int_t ehl_max_its =
      100;  // max number of iterations for elastic - fluid coupled problem

  il::int_t frac_front_max_its =
      30;  // max number of its for fracture front loop.

  double ehl_relaxation =
      0.5;  //  ehl under-relaxation parameters. must be in [0,1]
  double ehl_tolerance = 1.e-6;  // tolerance for stopping criteria for elastic
  // - fluid coupled problem
  double frac_front_tolerance =
      1.e-3;  // tolerance for stopping criteria on fracture front position

  // a string for an option to adapt the under-relaxation
};

}
#endif  // HFPX2_SIMULATIONPARAMETERS_H
