//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 04.05.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <cmath>
#include <il/Array.h>

// Inclusion from the project
#include "Permeability.h"

namespace hfp2d {

il::Array<double> permeability(Parameters_permeability &param,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - param -> structure that contains all the permeability parameters
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> K{d.size(), 0};

  for (il::int_t i = 0; i < K.size(); ++i) {

    K[i] = param.Init_permeab +
           (param.Incr_permeab * (1 - exp(-d[i] / param.d_wd)));
  }

  return K;
};
}