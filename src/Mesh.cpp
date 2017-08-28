//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "Mesh.h"
#include <il/linear_algebra.h>
#include <il/math.h>
#include <iostream>

namespace hfp2d {

// mesh class
void Mesh::set_values(il::Array2D<double> xy, il::Array2D<int> ien) {
  // check array dimensions ??? -> this is only for 1D mesh so far
  IL_EXPECT_FAST(xy.size(1) == 2);  // check array dimensions ?
  IL_EXPECT_FAST(ien.size(1) == 2);

  node_ = xy;           // list of coordinates of points in the mesh
  connectivity_ = ien;  //  connectivity array -
}
double Mesh::node(il::int_t k, il::int_t i) { return node_(k, i); }

int Mesh::connectivity(il::int_t k, il::int_t i) { return connectivity_(k, i); }

int Mesh::nelts() { return connectivity_.size(0); };

int Mesh::ncoor() { return node_.size(0); };

il::Array2D<double> Mesh::coor() { return node_; };



// needs to add function to add one or more elements ... (needs to have active
// and passive elements to track active/passive fractures etc.)
//
//
// could provide a default constructor for a straight fracture ?




// SOME UTILITIES HERE below -> To be moved in a separate file ??

il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}

// Function returning the segment characteristic from the matrix of the
// coordinates of the end points and the knowledge of the degree of
// interpolation p
//  work for segment mesh
// Inputs
// mesh object
// ne element number in the mesh to get characterstic from
// p order of the interpolation for that mesh
SegmentCharacteristic get_segment_DD_characteristic(Mesh mesh, il::int_t const ne,
                                                    il::int_t const p) {
  //  IL_ASSERT(Xs.size(0) == 2);
  //  IL_ASSERT(Xs.size(1) == 2);

  SegmentCharacteristic segment;

  // compute element size
  il::StaticArray<double, 2> xdiff, s, n, xmean, xaux;
  il::Array2D<double> Xcol{p + 1, 2, 0};
  il::StaticArray2D<double, 2, 2> R;
  il::StaticArray2D<double, 2, 2> Xs;

  Xs(0, 0) = mesh.node(mesh.connectivity(ne, 0), 0);
  Xs(0, 1) = mesh.node(mesh.connectivity(ne, 0), 1);

  Xs(1, 0) = mesh.node(mesh.connectivity(ne, 1), 0);
  Xs(1, 1) = mesh.node(mesh.connectivity(ne, 1), 1);

  xdiff[0] = Xs(1, 0) - Xs(0, 0);
  xdiff[1] = Xs(1, 1) - Xs(0, 1);

  segment.size = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

  // s=xdiff; // tangent vector
  s[0] = xdiff[0] / segment.size;
  s[1] = xdiff[1] / segment.size;
  n[0] = -1. * s[1];
  n[1] = s[0];  // normal vector

  segment.s = s;
  segment.n = n;

  segment.theta = acos(s[0] / sqrt(pow(s[0], 2) + pow(s[1], 2)));
  if (s[1] < 0) {
    segment.theta = -segment.theta;
  };

  xmean[0] = (Xs(1, 0) + Xs(0, 0)) / 2.;
  xmean[1] = (Xs(1, 1) + Xs(0, 1)) / 2.;
  segment.Xmid = xmean;

  switch (p) {
    case 1: {  // linear DD
      Xcol(0, 0) = -1. / sqrt(2.);
      Xcol(0, 1) = 0.;
      Xcol(1, 0) = 1. / sqrt(2.);
      Xcol(1, 1) = 0.;
    }; break;

    case 0: {
      Xcol(0, 0) = 0.;
      Xcol(0, 1) = 0.;
    }; break;
    default:
      std::cout << "error\n";  //  error
      break;
  };

// Returning the collocation point in the global frame

  R = rotation_matrix_2D(segment.theta);

  for (int i = 0; i < p + 1; ++i) {
    xaux[0] = (segment.size) * Xcol(i, 0) / 2.;
    xaux[1] = (segment.size) * Xcol(i, 1) / 2.;

    xaux = il::dot(R, xaux);

    Xcol(i, 0) = xaux[0] + xmean[0];
    Xcol(i, 1) = xaux[1] + xmean[1];
  }

  segment.CollocationPoints = Xcol;

  return segment;  // return structure with all we need on the segment.
}
//----------------------------------------------------

}