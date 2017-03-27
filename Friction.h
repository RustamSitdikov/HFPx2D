//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FRICTION_H
#define HFPX2D_FRICTION_H

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

// Inclusion from the project
#include "Mesh.h"

namespace hfp2d {

il::Array<double> exp_friction(LayerParameters1 &layer_parameters1,
                               LayerParameters2 &layer_parameters2,
                               LayerParameters3 &layer_parameters3,
                               const il::Array<il::int_t> &id_layers,
                               il::Array2D<int> Dofw,
                               const il::Array<double> &d, il::io_t);

il::Array<double> lin_friction(LayerParameters1 &layer_parameters1,
                               LayerParameters2 &layer_parameters2,
                               LayerParameters3 &layer_parameters3,
                               const il::Array<il::int_t> &id_layers,
                               il::Array2D<int> Dofw,
                               const il::Array<double> &d, il::io_t);
}

#endif // HFPX2D_FRICTION_H
