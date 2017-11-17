//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from Inside Loop library
#include <il/linear_algebra.h>

// Inclusion from the project
#include <src/Core/Utilities.h>

namespace hfp2d {

////////////////////////////// Some utilities //////////////////////////////////

void take_submatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
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

void set_submatrix(il::Array2D<double> &A, int i0, int i1,
                   const il::StaticArray2D<double, 2, 4> &B) {
  IL_EXPECT_FAST(i0 + B.size(0) <= A.size(0));
  IL_EXPECT_FAST(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}

////////////////////////////////////////////////////////////////////////////////
// Template for the function that return the index ( IndRow,IndColumn ) of a
// given value ("seek") in an array2D.
// arr2D -> array2D in which we want to find the index of  given value
// seek -> value for which we want to find out the index
// It return an array that contain the {N.row, N.col} of the seek value
il::Array<int> find_2d_integer(const il::Array2D<int> &arr2D, int seek) {
  il::Array<int> outp{2};

  for (int i = 0; i < arr2D.size(0); ++i) {
    for (int j = 0; j < arr2D.size(1); ++j) {
      if (arr2D(i, j) == seek) outp[0] = i, outp[1] = j;
    }
  }
  return outp;
}

////////////////////////////////////////////////////////////////////////////////
// Template for the fuunction that delete duplicates in an array of integers
il::Array<int> delete_duplicates_integer(const il::Array<int> &arr) {
  il::Array<int> out{};

  for (il::int_t i = 0; i < arr.size(); ++i) {
    bool already_there = false;
    for (il::int_t j = 0; j < out.size(); ++j) {
      if (arr[i] == out[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
      out.append(arr[i]);
    }
  }

  return out;
};

////////////////////////////////////////////////////////////////////////////////
// This function calculates the 2D euclidean distance between two points
// {x1,y1}, {x2,y2}
// Input -> x- y- coordinates of the two points
// Output -> double precision values that represents the euclidean distance
//           between the two aforementioned points
double euclidean_distance(double x1, double y1, double x2, double y2) {
  double dist;

  double x = x1 - x2;
  double y = y1 - y2;

  dist = (x * x) + (y * y);
  dist = sqrt(dist);

  return dist;
};
}
