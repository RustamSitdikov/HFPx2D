//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 12.12.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_DOMAINMESH_H
#define HFPX2D_DOMAINMESH_H

#include <limits>
#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp2d {

class DomainMesh {
 private:
  il::Array2D<double> nodes_;

  il::Array2D<il::int_t> connectivity_;

  il::Array<il::int_t>
      matid_;  // usually in this kind of domain wellMesh: material
  // is different for each element.. never known

  il::int_t nodes_per_elt_;

 public:
  ///////////////////////////////////////////////////
  // constructor ....
  ///////////////////////////////////////////////////

  DomainMesh(const il::Array2D<double> &nodes,
             const il::Array2D<il::int_t> &connectivity,
             const il::Array<il::int_t> &matid) {
    IL_EXPECT_FAST(connectivity.size(0) == matid.size());

    nodes_ = nodes;
    connectivity_ = connectivity;
    matid_ = matid;

    nodes_per_elt_ = connectivity.size(1);
  };

  ///////////////////////////////////////////////////
  // get functions
  ///////////////////////////////////////////////////

  // returning the matid corresponding to an element.
  il::int_t getmatid(il::int_t k) const { return matid_[k]; };

  il::int_t numberOfElts() const { return connectivity_.size(0); };

  ///////////////////////////////////////////////////
  // Methods
  ///////////////////////////////////////////////////

  il::Array<double> elementCentroid(il::int_t k) {
    // compute the xy coor of  element k centroid.
    il::Array<double> centroid(2, 0.);
    centroid[0] = 0.;
    centroid[1] = 0.;

    for (il::int_t i = 0; i < nodes_per_elt_; i++) {
      centroid[0] += nodes_(connectivity_(k, i), 0) / nodes_per_elt_;
      centroid[1] += nodes_(connectivity_(k, i), 1) / nodes_per_elt_;
    }
    return centroid;
  }

  il::Array2D<double> allElementCentroid() {
    // compute the xy coor of  element k centroid.
    il::Array2D<double> allcentroids{connectivity_.size(0), 2, 0.};
    for (il::int_t k = 0; k < connectivity_.size(0); k++) {
      allcentroids(k, 0) = 0.;
      allcentroids(k, 1) = 0.;
      for (il::int_t i = 0; i < nodes_per_elt_; i++) {
        allcentroids(k, 0) += nodes_(connectivity_(k, i), 0) / nodes_per_elt_;
        allcentroids(k, 1) += nodes_(connectivity_(k, i), 1) / nodes_per_elt_;
      }
    }
    return allcentroids;
  }

  // method locating the element in which the xy point is belonging too
  il::int_t locate(il::Array<double> &xy) {
    IL_EXPECT_FAST(xy.size() == 2);

    il::Array2D<double> centroids = allElementCentroid();

    // compute distance from all centroids.
    double dist;
    double min = std::numeric_limits<double>::max();
    il::int_t min_pos = -1;

    for (il::int_t i = 0; i < connectivity_.size(0); i++) {
      const double dx = centroids(i, 0) - xy[0];
      const double dy = centroids(i, 1) - xy[1];
      dist = dx * dx + dy * dy;
      if (dist < min) {
        min = dist;
        min_pos = i;
      }
    }

    // this will work for a regular Quad wellMesh but NOT necessarily for
    // unstructured triangular wellMesh
    return min_pos;  // let's be bold don t do any more checks  ! aie aie.
    // todo bulletproof this locate routine -> find all elements around via node
    // sharing, then check it is not in one of the element around.
    // this will return the first element found if xy is exactly a node of the
    // background wellMesh.
  }
};

// if Cartesian quad only : we could be more optimized.
//  il::Array<il::int_t>  LocateInQuad(il::Array2D<double> &xy){
//
//  IL_EXPECT_FAST(nodes_per_elt_==4); // we want quad here
//
//
//
//  };
};

#endif  // HFPX2D_DOMAINMESH_H
