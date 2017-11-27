//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 06.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_CREATEMESH_H
#define HFPX2DUNITTEST_CREATEMESH_H

// Inclusion from standard library
#include <cmath>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

// Inclusion from the project
#include "src/core/Mesh.h"

namespace hfp2d {

il::Array2D<double> createCustomMesh(
    double x_1,              // First point of the line, x coordinate
    double y_1,              // First point of the line, y coordinate
    double x_2,              // Second point of the line, x coordinate
    double y_2,              // Second point of the line, y coordinate
    il::int_t numElements,   // Number of elements to be generated
    il::int_t interpOrder);  // Element interpolation order

il::Array2D<il::int_t> createAutoConnectivity(il::int_t interpolationOrder,
                                              il::int_t numberOfElements);

il::Array2D<il::int_t> createAutoDisplacementDofHandle(
    il::int_t interpolationOrder, il::int_t numberOfElements);

il::Array2D<il::int_t> createAutoPressureDofHandle(il::int_t interpolationOrder,
                                                   il::int_t numberOfElements);
}
#endif  // HFPX2DUNITTEST_CREATEMESH_H
