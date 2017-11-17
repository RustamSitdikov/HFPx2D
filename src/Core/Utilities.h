//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 25.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_MATRIX_UTILITIES_H
#define HFPX2D_MATRIX_UTILITIES_H

// put here some utilites, especially the take of submatrix...
// different version could be developed depending on efficiency.

namespace hfp2d {

void take_submatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
                    const il::Array2D<double> &A);

void set_submatrix(il::Array2D<double> &A, int i0, int i1,
                   const il::StaticArray2D<double, 2, 4> &B);

//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta);

il::Array<int> find_2d_integer(const il::Array2D<int> &arr2D, int seek);

il::Array<int> delete_duplicates_integer(const il::Array<int> &arr);

double euclidean_distance(double x1, double y1, double x2, double y2);
}

#endif  // HFPX2D_MATRIX_UTILITIES_H
