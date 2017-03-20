//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <cmath>

// Inclusion from the project
#include "Dilatancy.h"

namespace hfp2d {

// Function that returns an array that contain dilatancy values (according to
// exponential dilatant hardening law)
il::Array<double> dilatancy(Parameters_dilatancy &param,
                            const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - param -> structure that contains all the dilatancy parameters we need
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> D{d.size(), 0};

  for (il::int_t i = 0; i < D.size(); ++i) {

    D[i] = param.Init_hydr_width + (param.Incr_dil * (1 - exp(-d[i] / param.d_wd)));
  }

  return D;
};

// Function that returns an array that contain the derivative w.r.t slip of
// dilatancy values (according to exponential dilatant hardening law)
il::Array<double> d_dilatancy(Parameters_dilatancy &param,
                              const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - param -> structure that contains all the dilatancy parameters we need
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> DD{d.size(), 0};

  for (il::int_t i = 0; i < DD.size(); ++i) {

    DD[i] = (param.Incr_dil / param.d_wd) * (exp(-d[i] / param.d_wd));
  }

  return DD;
};
}
