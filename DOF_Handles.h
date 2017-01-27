//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_DOF_HANDLES_H
#define HFPX2D_DOF_HANDLES_H

// Inclusion from Inside Loop library
#include <il/Array2D.h>

namespace hfp {

    il::Array2D<int>
    dofhandle_DG2D(const int dof_dim, const int Nelts, const int p);

}

#endif //HFPX2D_DOF_HANDLES_H
