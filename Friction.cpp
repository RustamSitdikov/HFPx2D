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
il::Array<double> exp_friction(double Peak_fric, double Resid_fric, double d_wf,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - Peak_fric -> peak friction coefficient
  //  - Resid_fric -> residual friction coefficient
  //  - d_wf -> slip dw for scaling (see exponential law in the report)
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};

  for (il::int_t i = 0; i < f.size(); ++i) {

    f[i] = Peak_fric - ((Peak_fric - Resid_fric) * (1 - exp(-(d[i] / d_wf))));
  }

  return f;
};

// Function that returns an array that contain friction coefficient (according
// to LINEAR friction weakening law)
il::Array<double> lin_friction(double Peak_fric, double Resid_fric, double d_wf,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - Peak_fric -> peak friction coefficient
  //  - Resid_fric -> residual friction coefficient
  //  - d_wf -> slip dw for scaling (see exponential law in the report)
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};
  double_t sl;
  sl = Peak_fric / d_wf;

  for (il::int_t i = 0; i < f.size(); ++i) {

    //    if (d[i] < d_wf) {
    //
    //      f[i] = Peak_fric - (sl * d[i]);
    //    } else  {
    //
    //      f[i] = Resid_fric;
    //    }

    if (d[i] < 0) {

      f[i] = 1;

    } else if (d[i] < ((Peak_fric-Resid_fric)/sl) ) {

      f[i] = Peak_fric - (sl * d[i]);

    } else {

      f[i] = Resid_fric;
    }
  }

  return f;
};
}
