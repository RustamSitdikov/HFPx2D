//
// This file is part of WellHDTest.
//
// Created by nikolski on 11/22/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef WELLHDTEST_WELLMESH_H
#define WELLHDTEST_WELLMESH_H

// include std libraries
#include <algorithm>
#include <cmath>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from hfp2d
#include <src/core/Fluid.h>
#include <src/core/SegmentData.h>
#include <src/core/Solution.h>
#include <src/core/Sources.h>
#include <src/core/Mesh.h>

namespace hfp2d {

//////////////////////////////////////////////////////////////////////////
// well wellMesh class; based on hfp2d::Mesh class
//////////////////////////////////////////////////////////////////////////
class WellMesh {

  // A wellMesh object for a well
  // P0 elements:
  // measured depth & true vertical depth are given at the NODES;
  // hydraulic radius & roughness are set at the CENTER of each element

 private:
  // MD / TVD - measured depth (curvilinear) / true vertical depth
  il::Array<double> md_;
  il::Array<double> tvd_;

  // "coordinates_" refer to x-y plane projection of the well
  // oriented as cos(inclination_)*(sin(azimuth), cos(azimuth)) by default
  // "connectivity_" (inherited) is used as is
  il::Array2D<double> coordinates_;
  il::Array2D<il::int_t> connectivity_;

  // elements sharing each node

  il::Array2D<il::int_t> node_adj_elt_;

  il::Array2D<il::int_t> edge_common_;

  // local sine of inclination inclination
  // (used to determine hydrostatic pressure gradient)
  il::Array<double> inclination_;

  double azimuth;  // it may be il::Array<double>

  // HD - hydraulic diameter   at cell center.
  il::Array<double> hd_;

  // dimensionless surface roughness
  il::Array<double> rough_;

  il::int_t  interpolation_order_ = 0;

 public:
  ////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS
  ////////////////////////////////////////////////////////////////////////

  // default constructor
  WellMesh(){};

  // normal constructor
  WellMesh(il::Array<double> &md, il::Array<double> &tvd, il::Array<double> &hd,
           double azimuth,  // todo: make it il::Array<double>
           il::Array<double> &rough) {
    // md:: array of measured depth
    // tvd : array of corresponding true vertical depth
    // hd: array of pipe diameter
    // rough:: array of pipe roughhness
    // azimuth:: double well azimuth w.r. to true north (in the y axis direction by convention)
    //           in degree !

    // check validity of inputs
    IL_EXPECT_FAST(md.size() > 1 && tvd.size() > 1);
    IL_EXPECT_FAST(tvd.size() == md.size());
    IL_EXPECT_FAST(hd.size() == md.size() - 1);

    il::int_t n_pts = md.size();
    il::int_t n_nod = n_pts;
    il::int_t n_elt = n_pts - 1;

    interpolation_order_ = 0;

    if (md[0] == 0.0 && tvd[0] == 0.0) {
      md_ = md;
      tvd_ = tvd;
      hd_ = hd;
    } else {
      ++n_nod;
      ++n_elt;
      md_ = il::Array<double>{n_nod};
      tvd_ = il::Array<double>{n_nod};
      hd_ = il::Array<double>{n_elt};
      // set inlet depth(s) to zero
      md_[0] = 0.0;
      tvd_[0] = 0.0;
      hd_[0] = hd[0];
      // append w. input data
      for (il::int_t n = 0; n < n_pts; ++n) {
        md_[n + 1] = md[n];
        tvd_[n + 1] = tvd[n];
        if (n < n_elt) {
          hd_[n + 1] = hd[n];
        }
      }
    }

    IL_EXPECT_FAST(rough.size() == n_elt);
    rough_ = rough;

    inclination_ = il::Array<double>{n_elt};

    coordinates_ = il::Array2D<double>{n_nod, 2};

    coordinates_(0, 0) = 0.0;
    coordinates_(0, 1) = 0.0;
    connectivity_ = il::Array2D<il::int_t>{n_elt, 2};

    for (il::int_t el = 0; el < n_elt; ++el) {
      // left (top) node
      il::int_t n_l = el;
      // right (bottom) node
      il::int_t n_r = el + 1;

      connectivity_(el, 0) = n_l;
      connectivity_(el, 1) = n_r;

      // cell (element) length
      double el_length = eltSize(el);

      // sine of inclination angle (zero -> horizontal)
      double sin_a = (tvd_[n_r] - tvd_[n_l]) / el_length;
      // decide if sin_a is better to store & use than asin(sin_a)
      inclination_[el] = sin_a;

      double cos_a = std::sqrt(1. - std::pow(sin_a, 2));

      // calculate the well orientation in x-y plane (north is toward y)
      // todo: make azimuth il::Array<double>
      double to_rad = il::pi / 180.;
      coordinates_(n_r, 0) =
          coordinates_(n_l, 0) + std::sin(azimuth*to_rad) * cos_a * el_length;
      coordinates_(n_r, 1) =
          coordinates_(n_l, 1) + std::cos(azimuth*to_rad) * cos_a * el_length;
    }

// build the nodal connected table...
    node_adj_elt_  = getNodalEltConnectivity(coordinates_.size(0), connectivity_);

    edge_common_ =il::Array2D<il::int_t>(md_.size() - 2, 2, 0);

    il::int_t k = 0;
    for (il::int_t i = 0; i < coordinates_.size(0); i++) {
      if (node_adj_elt_(i, 1) > -1) {
        edge_common_(k, 0) = node_adj_elt_(i, 0);
        edge_common_(k, 1) = node_adj_elt_(i, 1);
        k++;
      }
    }

  }

  ////////////////////////////////////////////////////////////////////////
  //        public interfaces
  ////////////////////////////////////////////////////////////////////////

  inline il::int_t numberOfElts() const { return connectivity_.size(0); }
//
  il::int_t numberOfNodes() const { return coordinates_.size(0); }
//
  il::Array2D<double> coordinates() const { return coordinates_; };
//
  double coordinates(il::int_t k, il::int_t i) const {
    return coordinates_(k, i);
  }
//
  il::Array2D<il::int_t> connectivity() const { return connectivity_; };
//
  il::int_t connectivity(il::int_t e, il::int_t i) const {
    return connectivity_(e, i);
  }
//  // nodal connectivity related
//  inline il::Array2D<il::int_t> nodeEltConnectivity() const {
//    return node_adj_elt_;
//  };

  il::Array2D<il::int_t> edgeCommon() const { return edge_common_;};



  il::Array<double> md() { return md_; }

  double md(il::int_t n) { return md_[n]; }

  il::Array<double> tvd() { return tvd_; }

  double tvd(il::int_t n) { return tvd_[n]; }

  il::Array<double> inclination() { return inclination_; }

  double inclination(il::int_t el) { return inclination_[el]; }

  il::Array<double> hd() { return hd_; }

  double hd(il::int_t el) { return hd_[el]; }

  il::Array<double> rough() { return rough_; }

  double rough(il::int_t el) { return rough_[el]; }



  ////////////////////////////////////////////////////////////////////////
  //        set functions
  ////////////////////////////////////////////////////////////////////////

  void setNodeAdjElts();

  ////////////////////////////////////////////////////////////////////////
  //        Methods
  ////////////////////////////////////////////////////////////////////////

  // get the size of a given element
  // (overrides the parent Mesh method
  // since curvilinear coordinate MD is used)
  double eltSize( il::int_t el);

  il::Array2D<il::int_t> getNodesSharing2Elts();


};

};

#endif  // WELLHDTEST_WELLMESH_H
