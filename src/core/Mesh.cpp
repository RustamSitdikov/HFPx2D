//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <src/core/Mesh.h>
#include <il/linear_algebra.h>

namespace hfp2d {

//////////////////////////////// CONSTRUCTORS (a.k.a. initializers) ////////////////////////////////
// Initialization without material properties and boundary conditions,
// as they are considered default (equal to zero) in the fracture mesh
Mesh::Mesh(const il::int_t interpolationOrder,
           const il::Array2D<double> &nodesCoordinates,
           const il::Array2D<il::int_t> &elementsConnectivity,
           const il::Array<il::int_t> &sourceIdentifier) {

  // Initial assertions to check data consistency
  // Non-zero number of nodes and coordinates must be 2D
  IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);

  // Non-zero number of elements and connectivity matrix shall have
  // as many columns as the interpolation order +1
  IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 &&
      elementsConnectivity.size(1) == interpolationOrder + 1);

  // Assignment of node and element sizes
  number_nodes_ = nodesCoordinates.size(0);
  number_elements_ = elementsConnectivity.size(0);

  // Assignment to the class members
  nodes_ = nodesCoordinates;           // list of coordinates of points in the mesh
  connectivity_ = elementsConnectivity;  //  connectivity array -

  // Sources are located at nodes. The source identifier is equal to -1, no source is applied.
  // Otherwise, the index of the source is saved in the vector (starting the numbering from 0)
  IL_EXPECT_FAST(sourceIdentifier.size() == number_nodes_);
  source_id_ = sourceIdentifier;

  //// INSERTING DEFAULT VALUES
  /// This is done for those variables that are required ANYWAY during computation.

  // filling the material_id_ and the fracture_id_ with the standard value
  for (il::int_t i = 0; i < number_elements_; i++) {
    fracture_id_[i] = 0;
    material_id_[i] = 0;
  }

  // filling the dof_handle_displacements (2D discontinuous galerkin)
  for (il::int_t i = 0, j; i < number_elements_; i++) {
    j = i * 2 * (interpolationOrder + 1);

    for (int k = 0; k < 2 * (interpolationOrder + 1); k++) {
      dof_handle_displacement_(i, k) = j + k;
    }

  }

  // fillind the dof_handle_pressure (continuous galerkin)
  for (int i = 0; i < number_elements_; ++i) {
    dof_handle_pressure_(i, 0) = i;
    dof_handle_pressure_(i, 1) = i + 1;
  }

  // is_tip_ default vector - beginning and end are tips
  is_tip_[0] = true;
  for (il::int_t i = 1; i < number_nodes_ - 1; i++) {
    is_tip_[i] = false;
  }
  is_tip_[number_nodes_] = true;

};

// Initialization as before but with dof_handles passed explicitly
Mesh::Mesh(const il::int_t interpolationOrder,
           const il::Array2D<double> &nodesCoordinates,
           const il::Array2D<il::int_t> &elementsConnectivity,
           const il::Array2D<il::int_t> &dofHandleDisplacement,
           const il::Array2D<il::int_t> &dofHandlePressure,
           const il::Array<il::int_t> &sourceIdentifier)
    : Mesh(interpolationOrder,
           nodesCoordinates,
           elementsConnectivity,
           sourceIdentifier) {

  // Let us deal with nodes, connectivity, source
  // dof_handles, fracture and material IDs are set to default values
//  this->Mesh(interpolationOrder,
//             nodesCoordinates,
//             elementsConnectivity,
//             sourceIdentifier);

  // dofHandle for displacement has to have as many rows as elements and
  // as many columns as nodes per element x displacements dofs per node
  IL_EXPECT_FAST(dofHandleDisplacement.size(0) == number_elements_ &&
      dofHandleDisplacement.size(1) == (interpolationOrder + 1) * 2);

  // dofHandle for pressure has to have as many rows as elements and
  // as many columns as nodes per element x pressure dofs per node
  IL_EXPECT_FAST(dofHandlePressure.size(0) == number_elements_ &&
      dofHandlePressure.size(1) == interpolationOrder + 1);

  dof_handle_displacement_ = dofHandleDisplacement;
  dof_handle_pressure_ = dofHandlePressure;

};

// Initialization with everything but tip identifier, for constant meshes
// which have the tip at the beginning and the end (default construction
// of the tip)
Mesh::Mesh(const il::int_t interpolationOrder,
           const il::Array2D<double> &nodesCoordinates,
           const il::Array2D<il::int_t> &elementsConnectivity,
           const il::Array2D<il::int_t> &dofHandleDisplacement,
           const il::Array2D<il::int_t> &dofHandlePressure,
           const il::Array<il::int_t> &fractureIdentifier,
           const il::Array<il::int_t> &materialIdentifier,
           const il::Array<il::int_t> &farStressCondId,
           const il::Array<il::int_t> &porePressCondId,
           const il::Array<il::int_t> &sourceIdentifier)
    : Mesh(interpolationOrder,
           nodesCoordinates,
           elementsConnectivity,
           dofHandleDisplacement,
           dofHandlePressure,
           sourceIdentifier) {

//  this->Mesh(interpolationOrder,
//             nodesCoordinates,
//             elementsConnectivity,
//             dofHandleDisplacement,
//             dofHandlePressure,
//             sourceIdentifier);

  // fracture identifier is an integer that determines to which "group" the element pertains
  IL_EXPECT_FAST(fractureIdentifier.size() == number_elements_);
  // material identifier is an integer that determines which rock failure
  // and flow/transport properties are applied to the element
  IL_EXPECT_FAST(materialIdentifier.size() == number_elements_);

  fracture_id_ = fractureIdentifier; // fracture ID vector
  material_id_ = materialIdentifier; // material ID vector


  // far field stress conditions are applied to the collocation points.
  // consequently, its identifier has the length equal to the number of collocation points,
  // that is number of elements x (interpolation order+1)
  IL_EXPECT_FAST(farStressCondId.size() == number_nodes_);
  // pore pressure condition is applied to the nodes, so this vector length is the number of nodes
  IL_EXPECT_FAST(porePressCondId.size() == number_nodes_);

  far_field_stress_condition_id_ = farStressCondId;
  pressure_condition_id_ = porePressCondId;

};

// Full mesh initialization
Mesh::Mesh(const il::int_t interpolationOrder,
           const il::Array2D<double> &nodesCoordinates,
           const il::Array2D<il::int_t> &elementsConnectivity,
           const il::Array2D<il::int_t> &dofHandleDisplacement,
           const il::Array2D<il::int_t> &dofHandlePressure,
           const il::Array<il::int_t> &fractureIdentifier,
           const il::Array<il::int_t> &materialIdentifier,
           const il::Array<il::int_t> &farStressCondId,
           const il::Array<il::int_t> &porePressCondId,
           const il::Array<il::int_t> &sourceIdentifier,
           const il::Array<bool> &tipIdentifier)
    : Mesh(interpolationOrder,
           nodesCoordinates,
           elementsConnectivity,
           dofHandleDisplacement,
           dofHandlePressure,
           fractureIdentifier,
           materialIdentifier,
           farStressCondId,
           porePressCondId,
           sourceIdentifier) {

//  this->Mesh(interpolationOrder,
//             nodesCoordinates,
//             elementsConnectivity,
//             dofHandleDisplacement,
//             dofHandlePressure,
//             fractureIdentifier,
//             materialIdentifier,
//             farStressCondId,
//             porePressCondId,
//             sourceIdentifier);

  // isTip vector is a boolean vector which value is true if the node is part is a tip of the fracture
  IL_EXPECT_FAST(tipIdentifier.size() == number_nodes_);
  is_tip_ = tipIdentifier;

};

Mesh::Mesh(const il::int_t interpolationOrder,
           const il::Array2D<double> &nodesCoordinates,
           const il::Array2D<il::int_t> &elementsConnectivity,
           const il::Array2D<il::int_t> &dofHandleDisplacement,
           const il::Array2D<il::int_t> &dofHandlePressure,
           const il::int_t fractureIdentifier,
           const il::int_t materialIdentifier,
           const il::int_t farStressCondId,
           const il::int_t porePressCondId,
           const il::String &injectionLocation) {

  // Initial assertions to check data consistency
  // Non-zero number of nodes and coordinates must be 2D
  IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);

  // Non-zero number of elements and connectivity matrix shall have
  // as many columns as the interpolation order +1
  IL_EXPECT_FAST(elementsConnectivity.size(0) > 0 &&
      elementsConnectivity.size(1) == interpolationOrder + 1);

  nodes_ = nodesCoordinates;
  connectivity_ = elementsConnectivity;

  number_nodes_ = nodesCoordinates.size(0);
  number_elements_ = elementsConnectivity.size(0);

  // dofHandle for displacement has to have as many rows as elements and
  // as many columns as nodes per element x displacements dofs per node
  IL_EXPECT_FAST(dofHandleDisplacement.size(0) == number_elements_ &&
      dofHandleDisplacement.size(1) == (interpolationOrder + 1) * 2);

  // dofHandle for pressure has to have as many rows as elements and
  // as many columns as nodes per element x pressure dofs per node
  IL_EXPECT_FAST(dofHandlePressure.size(0) == number_elements_ &&
      dofHandlePressure.size(1) == interpolationOrder + 1);

  dof_handle_displacement_ = dofHandleDisplacement;
  dof_handle_pressure_ = dofHandlePressure;

  for (il::int_t i = 0; i < number_nodes_; i++) {
    fracture_id_[i] = fractureIdentifier;
    material_id_[i] = materialIdentifier;
    far_field_stress_condition_id_[i] = farStressCondId;
    pressure_condition_id_[i] = porePressCondId;
    source_id_[i] = -1;
  }

  // dealing with the injection location
  if (injectionLocation == "start") { //place injection at the beginning
    source_id_[0] = 0
  } else if (injectionLocation == "end") { //place injection at the end
    source_id_[number_nodes_] = 0
  } else if (injectionLocation == "center") { //place injection in the center

    if (number_elements_ % 2 == 1) { // odd nodes, even number of elements, 1 node will be source
      source_id_[(number_elements_ - 1) / 2] = 0; // localize middle node and set source
    } else { // even nodes, odd elements
      source_id_[(number_elements_ - 1) / 2] = 0;
      source_id_[(number_elements_) / 2] = 0;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////




// mesh class
//void Mesh::loadMesh(il::Array2D<double> xy, il::Array2D<int> ien, il::Array<int> mat) {
//  // check array dimensions ?? -> this is only for 1D mesh so far
//  IL_EXPECT_FAST(xy.size(1) == 2);  // check array dimensions ?
//  IL_EXPECT_FAST(ien.size(1) == 2);
//
//  IL_EXPECT_FAST(ien.size(0) == mat.size());
//
//  nodes_ = xy;           // list of coordinates of points in the mesh
//  connectivity_ = ien;  //  connectivity array -
//  material_id_ = mat; // material ID array
//
//  // one could think of having a FracID ...
//}

//double Mesh::node(il::int_t k, il::int_t i) const { return nodes_(k, i); }
//
//int Mesh::connectivity(il::int_t k, il::int_t i) const { return connectivity_(k, i); }
//
//int Mesh::matid(il::int_t k) const { return material_id_[k]; }
//
//int Mesh::nelts() const { return connectivity_.size(0); } ;
//
//int Mesh::ncoor() const { return nodes_.size(0); };
//
//il::Array2D<double> Mesh::coor() const { return nodes_; };
//
//il::Array2D<int> Mesh::conn() const { return connectivity_; };
//
//il::Array<int> Mesh::matid() const { return material_id_; };


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
SegmentData get_segment_DD_data(const Mesh &mesh, il::int_t ne,
                                il::int_t p) {
  //  IL_ASSERT(Xs.size(0) == 2);
  //  IL_ASSERT(Xs.size(1) == 2);

  SegmentData segment;

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

  // mid point of the element
  xmean[0] = (Xs(1, 0) + Xs(0, 0)) / 2.;
  xmean[1] = (Xs(1, 1) + Xs(0, 1)) / 2.;
  segment.Xmid = xmean;

  switch (p) {
  case 1: {  // linear DD
    Xcol(0, 0) = -1. / sqrt(2.);
    Xcol(0, 1) = 0.;
    Xcol(1, 0) = 1. / sqrt(2.);
    Xcol(1, 1) = 0.;
  };
    break;

  case 0: {
    Xcol(0, 0) = 0.;
    Xcol(0, 1) = 0.;
  };
    break;
  default:std::cout << "error\n";  //  error
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