//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 04.05.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_PERMEABILITY_H
#define HFPX2D_PERMEABILITY_H

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

namespace hfp2d {

struct Parameters_permeability {

  // Initial value of permeability
  double Init_permeab;

  // Increment of permeability (difference between residual permeability (i.e.
  // at large slip) and initial permeability value)
  double Incr_permeab;

  // slip dw for scaling
  double d_wd;
};

il::Array<double> permeability(Parameters_permeability &param,
                               const il::Array<double> &d, il::io_t);
}

#endif // HFPX2D_PERMEABILITY_H
