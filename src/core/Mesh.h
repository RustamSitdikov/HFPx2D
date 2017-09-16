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
#include <cmath>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/container/1d/SmallArray.h>
#include <il/String.h>
#include <iostream>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

private:

  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> nodes_;
  // Connectivity matrix - size: number of elements x (order interpolation + 1)
  il::Array2D<il::int_t> connectivity_;

  // Interpolation order
  il::Array<il::int_t> interpolation_order_;

  // Dof handle matrices
  // for displacements - size: number of elements x 2dofs per node x (order interpolation + 1)
  il::Array2D<il::int_t> dof_handle_displacement_;
  // for pressure - size: number of nodes x 1dof per node x (order interpolation + 1)
  il::Array2D<il::int_t> dof_handle_pressure_;

  // Identifier number of the fracture - size: number of elements
  il::Array<il::int_t> fracture_id_;
  // Material identifier - size: number of elements
  il::Array<il::int_t> material_id_;

  // Condition identifiers for displacement and pressure, since stress is per collocation point
  // and pressure is for nodes
  il::Array<il::int_t> far_field_stress_condition_id_; // size: number of collocation nodes -----------------
  il::Array<il::int_t> pressure_condition_id_; // size: number of nodes ----------------

  // Source identifier, per node. It substitutes is_source.
  il::Array<il::int_t> source_id_;

  // Tip identifier
  il::Array<bool> is_tip_;

  // Characteristic "dimensions" of the class members
  // Number of nodes (default initialized at zero)
  il::int_t number_nodes_ = 0;
  // Number of elements (default initialized at zero)
  il::int_t number_elements_ = 0;
  // Number of interpolation order types
  il::int_t number_interpolations_=0;

  // Number of fractures
  il::int_t number_fractures_=0;
  // Number of material properties
  il::int_t number_materials_=0;

  // Number of displacement dofs per element
  il::int_t number_displacement_dofs_=0;
  // Number of pressure dofs per element
  il::int_t number_pressure_dofs_=0;

  // Number of far field stress conditions
  //il::int_t number_far_field_stress_conditions_=0;
  // Number of pore pressure conditions
  //il::int_t number_pore_pressure_conditions_=0;
  // Number of sources
  //il::int_t number_sources_=0;

  // Number of tips
  //il::int_t number_tips_=0;


public:

  //////////////////////////////// CONSTRUCTORS (a.k.a. initializers) ////////////////////////////////

  // Initialization without material properties and boundary conditions,
  // as they are considered default (equal to zero) in the fracture mesh
  Mesh(il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array<il::int_t> &sourceIdentifier);

  // Initialization as before but with dof_handles passed explicitly
  Mesh(il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t> &sourceIdentifier);

  // Initialization with everything but tip identifier, for constant meshes
  // which have the tip at the beginning and the end (default construction
  // of the tip)
  Mesh(il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t> &fractureIdentifier,
       const il::Array<il::int_t> &materialIdentifier,
       const il::Array<il::int_t> &farStressCondId,
       const il::Array<il::int_t> &porePressCondId,
       const il::Array<il::int_t> &sourceIdentifier);

  // Full mesh initialization
  Mesh(il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t> &fractureIdentifier,
       const il::Array<il::int_t> &materialIdentifier,
       const il::Array<il::int_t> &farStressCondId,
       const il::Array<il::int_t> &porePressCondId,
       const il::Array<il::int_t> &sourceIdentifier,
       const il::Array<bool> &tipIdentifier);

  Mesh(il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       il::int_t fractureIdentifier,
       il::int_t materialIdentifier,
       il::int_t farStressCondId,
       il::int_t porePressCondId,
       const il::String &sourceIdentifier);

  ////////////////////////////////////////////////////////////////////////////////////////////

/// SETTER
  void appendMesh(const Mesh &newMesh,
                        bool isJoined);


  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
                  const il::Array2D<il::int_t> &newElementsConnectivity,
                  const il::Array<il::int_t> &newMaterialIdentifier);

  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
                  const il::Array2D<il::int_t> &newElementsConnectivity,
                  const il::Array<il::int_t> &newMaterialIdentifier,
                  const il::Array<il::int_t> &newFractureIdentifier);

  void appendNodeToMeshTip(il::int_t mesh_node,
                              double x_new,
                              double y_new);

  /// GETTER
  // Read sizes of matrices
  il::int_t numberOfNodes() const { return number_nodes_; };
  il::int_t numberOfElements() const { return number_elements_; };
  il::int_t numberOfInterpolations() const { return number_interpolations_; };
  il::int_t numberOfFractures() const { return number_fractures_; };
  il::int_t numberOfMaterials() const { return number_materials_; };
  il::int_t numberOfDisplacementDofsPerElement() const { return number_displacement_dofs_; };
  il::int_t numberOfPressureDofsPerElement() const { return number_pressure_dofs_; };
  il::int_t numberOfFarFieldStressConditions() const { return number_far_field_stress_conditions_; };
  il::int_t numberOfPorePressureConditions() const { return number_pore_pressure_conditions_; };
  il::int_t numberOfSources() const { return number_sources_; };
  il::int_t numberOfTips() const { return number_tips_; };

  // Read the X coordinate of a node
  double X(il::int_t k) const { return nodes_(k, 0); }
  // Read the Y coordinate of a node
  double Y(il::int_t k) const { return nodes_(k, 1); }

  // Read a particular element of the node coordinates
  double node(il::int_t k, il::int_t i) const { return nodes_(k, i); }

  il::Array<il::int_t> elemConnectivity(il::int_t k) {
    il::Array<il::int_t> temp(connectivity_.size(1));

    for (il::int_t i = 0; i < connectivity_.size(1); i++) {
      temp[i] = connectivity_(k, i);
    }
  };

  il::int_t connectivity(il::int_t k, il::int_t i) const { return connectivity_(k, i); }

  il::int_t matid(il::int_t k) const { return material_id_[k]; }

  il::int_t nelts() const { return connectivity_.size(0); };

  il::int_t ncoor() const { return nodes_.size(0); };

  il::Array2D<double> coor() const { return nodes_; };

  il::Array2D<il::int_t> conn() const { return connectivity_; };

  il::Array<il::int_t> matid() const { return material_id_; };


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

}

#endif  // HFPX2D_MESH_H
