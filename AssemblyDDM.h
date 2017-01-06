//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_ASSEMBLYDDM_H
#define HFPX2D_ASSEMBLYDDM_H

#include <il/Array2D.h>
#include "Mesh.h"

void BasicAssembly(il::Array2D<double> &Kmat, Mesh mesh, il::Array2D<int> id, int p , double Ep );


#endif //HFPX2D_ASSEMBLYDDM_H
