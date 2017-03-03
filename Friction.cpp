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
#include "Friction.h"

namespace hfp2d {

// Function that returns an array that contain friction coefficient (according
// to EXPONENTIAL friction weakening law)
il::Array<double> exp_friction(Parameters_friction &param,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - param -> structure that contains all the friction parameters we need
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};

  for (il::int_t i = 0; i < f.size(); ++i) {

    f[i] = param.Peak_fric_coeff -
           ((param.Peak_fric_coeff - param.Resid_fric_coeff) *
            (1 - exp(-(d[i] / param.d_wf))));
  }

  return f;
};

// Function that returns an array that contain friction coefficient (according
// to LINEAR friction weakening law)
il::Array<double> lin_friction(Parameters_friction &param,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - param -> structure that contains all the friction parameters we need
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};
  double_t sl;
  sl = param.Peak_fric_coeff / param.d_wf;

  for (il::int_t i = 0; i < f.size(); ++i) {

    if (d[i] < ((param.Peak_fric_coeff - param.Resid_fric_coeff) / sl)) {

      f[i] = param.Peak_fric_coeff - (sl * d[i]);

    } else {

      f[i] = param.Resid_fric_coeff;
    }
  }

  return f;
};
}
