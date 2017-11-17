//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>



#include <src/core/Utilities.h>

namespace hfp2d{

// Some utilities //
void takeSubMatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
                   const il::Array2D<double> &A) {
  IL_EXPECT_FAST((i1 - i0 + 1) == sub.size(0));
  IL_EXPECT_FAST((j1 - j0 + 1) == sub.size(1));

  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      sub(i - i0, j - j0) = A(i, j);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

void setSubMatrix(il::Array2D<double> &A, int i0, int i1,
                  const il::StaticArray2D<double, 2, 4> &B) {
  IL_EXPECT_FAST(i0 + B.size(0) <= A.size(0));
  IL_EXPECT_FAST(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
}


//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotationMatrix2D(double theta) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}



}
