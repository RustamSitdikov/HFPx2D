//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_MATRIX_UTILITIES_H
#define HFPX2D_MATRIX_UTILITIES_H

// put here some utilites, especially the take of submatrix...
// different version could be developed depending on efficiency.

namespace hfp2d {

void takeSubMatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
                   const il::Array2D<double> &A);

void setSubMatrix(il::Array2D<double> &A, int i0, int i1,
                  const il::StaticArray2D<double, 2, 4> &B);


//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotationMatrix2D(double theta);


}







#endif //HFPX2D_MATRIX_UTILITIES_H
