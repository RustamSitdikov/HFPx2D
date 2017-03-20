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

    f[i] = param.Peak_fric_coeff_layer1 -
           ((param.Peak_fric_coeff_layer1 - param.Resid_fric_coeff_layer1) *
            (1 - exp(-(d[i] / param.d_wf_layer1))));
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
  double_t sl1;
  sl1 = param.Peak_fric_coeff_layer1 / param.d_wf_layer1;
  double_t sl2;
  sl2 = param.Peak_fric_coeff_layer2 / param.d_wf_layer2;
  double_t sl3;
  sl3 = param.Peak_fric_coeff_layer3 / param.d_wf_layer3;

  for (il::int_t i = 0; i < f.size() / 2; ++i) {
    if (d[i] < ((param.Peak_fric_coeff_layer1 - param.Resid_fric_coeff_layer1) / sl1)) {
      f[i] = param.Peak_fric_coeff_layer1 - (sl1 * d[i]);
    } else {
      f[i] = param.Resid_fric_coeff_layer1;
    }
  }

  for (il::int_t i = f.size() / 2; i < (f.size() - (f.size() / 4)); ++i) {
    if (d[i] < ((param.Peak_fric_coeff_layer2 - param.Resid_fric_coeff_layer2) / sl2)) {
      f[i] = param.Peak_fric_coeff_layer2 - (sl2 * d[i]);
    } else {
      f[i] = param.Resid_fric_coeff_layer2;
    }
  }

  for (il::int_t i = (f.size() - (f.size() / 4)); i < f.size(); ++i) {
    if (d[i] < ((param.Peak_fric_coeff_layer3 - param.Resid_fric_coeff_layer3) / sl3)) {
      f[i] = param.Peak_fric_coeff_layer3 - (sl3 * d[i]);
    } else {
      f[i] = param.Resid_fric_coeff_layer3;
    }
  }

  return f;
};
}
