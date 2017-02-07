//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_DOF_HANDLES_H
#define HFPX2D_DOF_HANDLES_H

// Inclusion from Inside Loop library
#include <il/Array2D.h>

namespace hfp2d {

il::Array2D<int> dofhandle_dg_full2d(const int dof_dim, const int Nelts,
                                     const int p);

il::Array2D<int> dofhandle_dg(const int dof_dim, const int Nelts);

il::Array2D<int> dofhandle_cg2d(const int dof_dim, const int Nelts);
}

#endif // HFPX2D_DOF_HANDLES_H
