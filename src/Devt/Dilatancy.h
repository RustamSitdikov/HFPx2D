//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_DILATANCY_H
#define HFPX2D_DILATANCY_H

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

namespace hfp2d {

struct Parameters_dilatancy {

  // Initial value of hydraulic width
  double Init_hydr_width;
  // Increment of dilatancy (difference between residual/peak
  // dilatancy and initial dilatancy value)
  double Incr_dil;
  // slip dw for scaling (see dilatancy law in the report)
  double d_wd;
};

il::Array<double> dilatancy(Parameters_dilatancy &param,
                            const il::Array<double> &d, il::io_t);

il::Array<double> d_dilatancy(Parameters_dilatancy &param,
                              const il::Array<double> &d, il::io_t);
}

#endif // HFPX2D_DILATANCY_H