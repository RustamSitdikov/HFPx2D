//
// This file is part of WellHDTest.
//
// Created by nikolski on 11/22/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <il/linear_algebra/dense/norm.h>
#include <src/wellbore/WellMesh.h>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
//        well wellMesh class methods
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// get the size of an element
double WellMesh::eltSize(const il::int_t el) {
  // left node
  il::int_t n_l = connectivity_(el, 0);
  // right node
  il::int_t n_r = connectivity_(el, 1);
  // return distance
  return std::fabs(md_[n_r] - md_[n_l]);
}


}