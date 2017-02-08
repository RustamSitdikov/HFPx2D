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

void set_submatrix(il::Array2D<double>& A, int i0, int i1, const il::StaticArray2D<double,2,4>& B);
void take_submatrix(il::Array2D<double>& sub, int i0, int i1, int j0, int j1, const il::Array2D<double>& A);

void BasicAssembly(il::Array2D<double> &Kmat, Mesh mesh, il::Array2D<int> id, int p , double Ep );


#endif //HFPX2D_ASSEMBLYDDM_H
