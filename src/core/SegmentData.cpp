//
// Created by lorenzo on 9/14/17.
//

#include "SegmentData.h"


namespace hfp2d {

//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}


}