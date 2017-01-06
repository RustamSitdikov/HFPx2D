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

#include <il/Array2D.h>

void dofhandle_DG2D(il::Array2D<int> &dofhandle, int dof_dim, int nelts, int p);

#endif //HFPX2D_DOF_HANDLES_H
