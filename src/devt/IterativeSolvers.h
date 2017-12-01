//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 29.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_ITERATIVESOLVERS_H
#define HFPX2D_ITERATIVESOLVERS_H

// Inclusion from the Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include <src/core/Mesh.h>

namespace hfp2d {

il::Array<double> jacobiIterativeSolver(il::Array2D<double> &matrix_of_coeff,
                                        il::Array<double> &right_hand_side,
                                        il::Array<double> &initial_solution,
                                        double relax_parameter,
                                        il::int_t numb_iterations);

il::Array<double> gaussSeidelIterativeSolver(
    il::Array2D<double> &matrix_of_coeff, il::Array<double> &right_hand_side,
    il::Array<double> &initial_solution, double relax_parameter,
    il::int_t numb_iterations);
}

#endif  // HFPX2D_ITERATIVESOLVERS_H
