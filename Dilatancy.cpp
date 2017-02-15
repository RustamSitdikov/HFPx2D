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
il::Array<double> dilatancy(double Init_dil, double Incr_dil,
                            double d_wd, const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - Init_dil -> Initial value of dilatancy
  //  - Incr_dil -> Increment of dilatancy (difference between residual/peak
  //  dilatancy and initial dilatancy value)
  //  - d_wd -> slip dw for scaling (see dilatancy law in the report)
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> D{d.size(), 0};

  for (il::int_t i = 0; i < D.size(); ++i) {

    D[i] = Init_dil + (Incr_dil * (1 - exp(-d[i] / d_wd)));
  }

  return D;
};

// Function that returns an array that contain the derivative w.r.t slip of
// dilatancy values (according to exponential dilatant hardening law)
il::Array<double> d_dilatancy(double Incr_dil, double d_wd,
                             const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - Incr_dil -> Increment of dilatancy (difference between residual/peak
  //  dilatancy and initial dilatancy value)
  //  - d_wd -> slip dw for scaling (see dilatancy law in the report)
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> DD{d.size(), 0};

  for (il::int_t i = 0; i < DD.size(); ++i) {

    DD[i] = (Incr_dil / d_wd) * (exp(-d[i] / d_wd));
  }

  return DD;
};
}
