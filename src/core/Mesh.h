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

  il::int_t numNodes_;
  il::int_t numElements_;
  il::Array2D<double> nodes_;
  il::Array2D<il::int_t> connectivity_;
  il::Array<il::int_t> matId_;
  il::Array<il::int_t> fracId_;

 public:

  //////////////////////////////// CONSTRUCTORS (a.k.a. initializers) ////////////////////////////////

  // First constructor: it initializes a mesh by using an array with the node coordinates, an
  // array with the connectivity of the elements and the corresponding material identifier.
  // Being the fracture identifier not specified, for this elements it is set as 0.
  void init1DMesh(const il::Array2D<double>& nodesCoordinates,
                  const il::Array2D<il::int_t>& elementsConnectivity,
                  const il::Array<il::int_t>& materialIdentifier)
  {

    // Initial assertions to check data consistency
    IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 && elementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(materialIdentifier.size() == elementsConnectivity.size(0));


    // Assignment to the class members
    nodes_ = nodesCoordinates;           // list of coordinates of points in the mesh
    connectivity_ = elementsConnectivity;  //  connectivity array -
    matId_ = materialIdentifier; // material ID array

    // Assignment of the sizes
    numNodes_=nodes_.size(0);
    numElements_=connectivity_.size(0);

    // Assignment of the fracture Identifier
    fracId_.resize(numElements_);
    for(il::int_t i=0; i<numElements_; i++){
      fracId_[i] = 0;
    }

  };


  // Initialization of a mesh class variable from data with nodes, connectivity and materials,
  // AND USING fracture identifier.
  void init1DMesh(const il::Array2D<double> &nodesCoordinates,
                  const il::Array2D<il::int_t> &elementsConnectivity,
                  const il::Array<il::int_t> &materialIdentifier,
                  const il::Array<il::int_t> &fractureIdentifier)
  {
    // Initial assertions to check data consistency
    IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 && elementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(materialIdentifier.size() == elementsConnectivity.size(0));
    IL_EXPECT_FAST(fractureIdentifier.size() == elementsConnectivity.size(0));

    // Assignment to the class members
    nodes_ = nodesCoordinates;           // list of coordinates of points in the mesh
    connectivity_ = elementsConnectivity;  //  connectivity array -
    matId_ = materialIdentifier; // material ID array
    fracId_ = fractureIdentifier;

    // Assignment of the sizes
    numNodes_=nodes_.size(0);
    numElements_=connectivity_.size(0);

  };


  /// SETTER
  void appendMesh(hfp2d::Mesh &newMesh){

    const il::int_t oldNumNodes = numNodes_;
    const il::int_t oldNumElements = numElements_;

    numNodes_ = numNodes_ + newMesh.numNodes_;
    numElements_ = numElements_ + newMesh.numElements_;

    nodes_.resize(numNodes_,2); // here the 2 is because we are working in 2D
    for(il::int_t i = oldNumNodes; i < numNodes_; i++){
      nodes_(i,0)=newMesh.nodes_(i,0);
      nodes_(i,1)=newMesh.nodes_(i,1);
    }

    // connectivity resize: here the 2 is because the element is 1D with 2 nodes
    // but it should rather be the order of the polynomia on the element plus 1
    connectivity_.resize(numElements_,2);
    matId_.resize(numElements_);
    fracId_.resize(numElements_);

    // add the values of the new mesh
    for(il::int_t i = oldNumElements; i < numElements_; i++){

      connectivity_(i,0)=newMesh.connectivity_(i,0);
      connectivity_(i,1)=newMesh.connectivity_(i,1);

      matId_[i]=newMesh.matId_[i];
      fracId_[i]=newMesh.fracId_[i];

    }

  };

  void appendMesh(const il::Array2D<double>& newNodesCoordinates,
                  const il::Array2D<il::int_t>& newElementsConnectivity,
                  const il::Array<il::int_t>& newMaterialIdentifier)
  {

    // Initial assertions to check data consistency
    IL_EXPECT_FAST(newNodesCoordinates.size(0) > 0 && newNodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(newElementsConnectivity.size(0) > 0 && newElementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(newMaterialIdentifier.size() == newElementsConnectivity.size(0));

    const il::int_t oldNumNodes = numNodes_;
    const il::int_t oldNumElements = numElements_;

    numNodes_ = numNodes_ + newNodesCoordinates.size(0);
    numElements_ = numElements_ + newElementsConnectivity.size(0);

    nodes_.resize(numNodes_,2); // here the 2 is because we are working in 2D
    for(il::int_t i = oldNumNodes; i < numNodes_; i++){
      nodes_(i,0)=newNodesCoordinates(i,0);
      nodes_(i,1)=newNodesCoordinates(i,1);
    }

    // connectivity resize: here the 2 is because the element is 1D with 2 nodes
    // but it should rather be the order of the polynomia on the element plus 1
    connectivity_.resize(numElements_,2);
    matId_.resize(numElements_);
    fracId_.resize(numElements_);

    // add the values of the new mesh
    for(il::int_t i = oldNumElements; i < numElements_; i++){

      connectivity_(i,0)=newElementsConnectivity(i,0);
      connectivity_(i,1)=newElementsConnectivity(i,1);

      matId_[i]=newMaterialIdentifier[i];
      fracId_[i]=0;

    }

  };


  void appendMesh(const il::Array2D<double>& newNodesCoordinates,
                  const il::Array2D<il::int_t>& newElementsConnectivity,
                  const il::Array<il::int_t>& newMaterialIdentifier,
                  const il::Array<il::int_t> &newFractureIdentifier)
  {

    // Initial assertions to check data consistency
    IL_EXPECT_FAST(newNodesCoordinates.size(0) > 0 && newNodesCoordinates.size(1) == 2);
    IL_EXPECT_FAST(newElementsConnectivity.size(0) > 0 && newElementsConnectivity.size(1) == 2);
    IL_EXPECT_FAST(newMaterialIdentifier.size() == newElementsConnectivity.size(0));
    IL_EXPECT_FAST(newFractureIdentifier.size() == newElementsConnectivity.size(0));

    const il::int_t oldNumNodes = numNodes_;
    const il::int_t oldNumElements = numElements_;

    numNodes_ = numNodes_ + newNodesCoordinates.size(0);
    numElements_ = numElements_ + newElementsConnectivity.size(0);

    nodes_.resize(numNodes_,2); // here the 2 is because we are working in 2D
    for(il::int_t i = oldNumNodes; i < numNodes_; i++){
      nodes_(i,0)=newNodesCoordinates(i,0);
      nodes_(i,1)=newNodesCoordinates(i,1);
    }

    // connectivity resize: here the 2 is because the element is 1D with 2 nodes
    // but it should rather be the order of the polynomia on the element plus 1
    connectivity_.resize(numElements_,2);
    matId_.resize(numElements_);
    fracId_.resize(numElements_);

    // add the values of the new mesh
    for(il::int_t i = oldNumElements; i < numElements_; i++){

      connectivity_(i,0)=newElementsConnectivity(i,0);
      connectivity_(i,1)=newElementsConnectivity(i,1);

      matId_[i]=newMaterialIdentifier[i];
      fracId_[i]=newFractureIdentifier[i];

    }

  };

  /// GETTER
  // Read the X coordinate of a node
  double X(il::int_t k) const { return nodes_(k,0); }
  // Read the Y coordinate of a node
  double Y(il::int_t k) const { return nodes_(k,1); }

  // Read a particular element of the node coordinates
  double node(il::int_t k, il::int_t i) const { return nodes_(k, i); }

  // Read sizes of matrices
  il::int_t numberOfNodes() const { return numNodes_; };
  il::int_t numberOfElements() const { return numElements_; };







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

  il::Array2D<il::int_t> conn() const { return connectivity_; };

  il::Array<il::int_t> matid() const { return matId_; };


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
