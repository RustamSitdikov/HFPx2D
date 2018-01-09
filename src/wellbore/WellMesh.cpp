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

////////////////////////////////////////////////////////////////////////////////
// method to get elements sharing each node
// (i.e. nodal_connectivity for all nodes with 2 neighbours)
void WellMesh::setNodeAdjElts() {
  node_adj_elt_ = getNodalEltConnectivity(numberOfNodes(), connectivity_);
}

////////////////////////////////////////////////////////////////////////////////
il::Array2D<il::int_t> WellMesh::getNodesSharing2Elts() {
  // case of only 2 nodes for now.....
  // needed for building FD matrix  of WB .....
  // ONLY WORK ON MESH where nodes have maximum of 2 adjacaent elements
  // hardcoded for a continuous wellbore mesh.

  IL_EXPECT_FAST(connectivity_.size(1) == 2);

  il::Array2D<il::int_t> temp(numberOfNodes() - 2, 2, 0);

  il::int_t k = 0;
  for (il::int_t i = 0; i < coordinates_.size(0); i++) {
    if (node_adj_elt_(i, 1) > -1) {
      temp(k, 0) = node_adj_elt_(i, 0);
      temp(k, 1) = node_adj_elt_(i, 1);
      k++;
    }
  }
  return temp;
}
}