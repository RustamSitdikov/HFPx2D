//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_MESH_H
#define HFPX2D_MESH_H

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

 private:

  il::Array2D<double> node_;
  il::Array2D<int> connectivity_;
  il::Array<int> matid_ ;

 public:
  void set_values(il::Array2D<double>, il::Array2D<int>,il::Array<int>);

  double node(il::int_t k, il::int_t i) const;

  int connectivity(il::int_t k, il::int_t i) const;

  int matid(il::int_t k) const;

  int nelts() const;

  int ncoor() const;

  il::Array2D<double> coor() const;

  il::Array2D<int> conn() const;

  il::Array<int> matid() const;

};


////////////////////////////////////////////////////////////////////////////////
struct SegmentData {
  double size;
  double theta;  // angle w.r. to e_1
  // unit normal to segment in global system of coordinates
  il::StaticArray<double, 2> n;
  // unit tangent to segment in global system of coordinates
  il::StaticArray<double, 2> s;
  // segment mid points coordinates.
  il::StaticArray<double, 2> Xmid;
  // collocation points in global system of coordinates
  il::Array2D<double> CollocationPoints;
};

il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta);

SegmentData get_segment_DD_data(const Mesh &mesh,
                                il::int_t ne, il::int_t p);

}
////////////////////////////////////////////////////////////////////////////////

#endif  // HFPX2D_MESH_H
