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
#include <il/container/1d/SmallArray.h>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

 private:

  il::Array2D<double> nodes_;
  il::Array2D<int> connectivity_;
  il::Array<int> matId_;
  il::Array<int> fracId_;

 public:

  //////////////////////////////// Constructors ////////////////////////////////
  // Constructor from data with nodes, connectivity and materials, BUT NOT fracture identifier
  void load1DMesh(il::Array2D<double> nodesCoordinates,
                  il::Array2D<int> elementsConnectivity,
                  il::Array<int> materialIdentifier)
  {
    // Initial assertions to check data consistency
    IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 && elementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(materialIdentifier.size() == elementsConnectivity.size(0));

    // Assignment to the class members
    nodes_ = nodesCoordinates;           // list of coordinates of points in the mesh
    connectivity_ = elementsConnectivity;  //  connectivity array -
    matId_ = materialIdentifier; // material ID array
  };

  // Constructor from data with nodes, connectivity and materials, AND USING fracture identifier
  void load1DMesh(il::Array2D<double> nodesCoordinates,
                  il::Array2D<int> elementsConnectivity,
                  il::Array<int> materialIdentifier,
                  il::Array<int> fractureIdentifier)
  {
    // Initial assertions to check data consistency
    IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 && elementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(materialIdentifier.size() == elementsConnectivity.size(1));
    IL_EXPECT_FAST(fractureIdentifier.size() == elementsConnectivity.size(1));

    // Assignment to the class members
    nodes_ = nodesCoordinates;           // list of coordinates of points in the mesh
    connectivity_ = elementsConnectivity;  //  connectivity array -
    matId_ = materialIdentifier; // material ID array
    fracId_ = fractureIdentifier;
  };


  /// GETTER
  // Read the X coordinate of a node
  double X(il::int_t k) const { return nodes_(k,0); }
  // Read the Y coordinate of a node
  double Y(il::int_t k) const { return nodes_(k,1); }

  // Read a particular element of the node coordinates
  double node(il::int_t k, il::int_t i) const { return nodes_(k, i); }

  il::Array<il::int_t> elemConnectivity(il::int_t k)
  {
    il::Array<il::int_t> temp(connectivity_.size(1));

    for(il::int_t i=0; i<connectivity_.size(1); i++)
    {
      temp[i]=connectivity_(k,i);
    }
  };

  int connectivity(il::int_t k, il::int_t i) const { return connectivity_(k, i); }

  int matid(il::int_t k) const { return matId_[k]; }

  int nelts() const { return connectivity_.size(0); } ;

  int ncoor() const { return nodes_.size(0); };

  il::Array2D<double> coor() const { return nodes_; };

  il::Array2D<int> conn() const { return connectivity_; };

  il::Array<int> matid() const { return matId_; };


//  double node(il::int_t k, il::int_t i) const;
//
//  int connectivity(il::int_t k, il::int_t i) const;
//
//  int matid(il::int_t k) const;
//
//  int nelts() const;
//
//  int ncoor() const;
//
//  il::Array2D<double> coor() const;
//
//  il::Array2D<int> conn() const;
//
//  il::Array<int> matid() const;

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
