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

// include std libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/String.h>
#include <il/base.h>
#include <il/container/1d/SmallArray.h>
#include <il/linear_algebra.h>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

 private:
  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> nodes_;

  // Connectivity matrix - size: number of elements x (order interpolation + 1)
  il::Array2D<il::int_t> connectivity_;

  // Interpolation order
  il::int_t interpolation_order_;

  // Dof handle matrices
  // for displacements - size: number of elements x 2dofs per node x (order
  // interpolation + 1)
  il::Array2D<il::int_t> dof_handle_displacement_;
  // for pressure - size: number of nodes x 1dof per node x (order interpolation
  // + 1)
  il::Array2D<il::int_t> dof_handle_pressure_;

  // Identifier number of the fracture - size: number of elements
  il::Array<il::int_t> fracture_id_;

  // Material identifier - size: number of elements
  il::Array<il::int_t> material_id_;
  // Material identifier - size: number of elements
  il::Array<il::int_t> condition_id_;

 public:
  /////////////////////          CONSTRUCTORS (a.k.a. initializers)
  ///////////////////////////////////
// todo: naming of the different entities are not consistent AND TOO LONG

  Mesh(){};  // TODO: remove empty initialization of mesh class variables if
             // possible.
  // In any case, we need always to have temporal variable to construct a mesh.

  // Constructor with only nodes and elements. TODO: to be substituted in the
  // Griffith test examples
  Mesh(const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity) {
    nodes_ = nodesCoordinates;
    connectivity_ = elementsConnectivity;
  };

  Mesh(const il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &displacementDOFHandle) {
    interpolation_order_ = interpolationOrder;
    nodes_ = nodesCoordinates;
    connectivity_ = elementsConnectivity;
    dof_handle_displacement_ = displacementDOFHandle;
  };

  Mesh(const il::int_t interpolationOrder,
       const il::Array2D<double> &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity, // connectivity is enough
       const il::Array2D<il::int_t> &displ_dof_handle, // should be dd not displ.
       const il::Array2D<il::int_t> &press_dof_handle,
       const il::Array<il::int_t> &fractureID,
       const il::Array<il::int_t> &materialID,
       const il::Array<il::int_t> &conditionID);

  ////////////////////////////////////////////////////////////////////////////////////////////
// todo : do comments each methods !
// todo: change interpolationOrder to DDOrder as it relates to the interpolation of Displacement Discontinuities
  /// SETTER
  void appendMesh(const Mesh &newMesh, bool isJoined);

  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
                  const il::Array2D<il::int_t> &newElementsConnectivity,
                  const il::Array<il::int_t> &newMaterialIdentifier);

  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
                  const il::Array2D<il::int_t> &newElementsConnectivity,
                  const il::Array<il::int_t> &newMaterialIdentifier,
                  const il::Array<il::int_t> &newFractureIdentifier);

  void appendNodeToMeshTip(il::int_t mesh_node, double x_new, double y_new);

  /// GETTER FUNCTIONS

  // please avoid long names ... number-> Nn
  // Displ should be DD
  il::int_t numberOfNodes() const { return nodes_.size(0); }
  il::int_t numberOfElements() const { return connectivity_.size(0); }
  il::int_t interpolationOrder() const { return interpolation_order_; }
  il::int_t numberOfFractures() const {
    return (*std::max_element(fracture_id_.begin(), fracture_id_.end()) + 1);
  }
  il::int_t numberOfMaterials() const {
    return (*std::max_element(material_id_.begin(), material_id_.end()) + 1);
  }
  il::int_t numberOfConditions() const {
    return (*std::max_element(condition_id_.begin(), condition_id_.end()) + 1);
  }
  il::int_t numberOfDisplDofsPerElement() const {
    return dof_handle_displacement_.size(1);
  }
  il::int_t numberOfPressDofsPerElement() const {
    return dof_handle_pressure_.size(1);
  }

  // ??? what is this ?
  il::int_t numberOfPressDofs() const {
    return (numberOfElements() * interpolationOrder() + numberOfFractures()); // this is not necessarily true
  }

  il::int_t numberOfDisplDofs() const {
    return (numberOfElements() * (interpolation_order_ + 1) * 2);
  }

  //  what is this mess below ?

  // Read the X coordinate of a node
  double X(il::int_t k) const { return nodes_(k, 0); }
  // Read the Y coordinate of a node
  double Y(il::int_t k) const { return nodes_(k, 1); }

  // Read a particular element of the node coordinates
  double node(il::int_t k, il::int_t i) const { return nodes_(k, i); }
  // todo: restructure "node" method in order to have mesh.node(i).X for the x
  // coordinate of node i ??  NO NO please. let's work with array

  il::Array<il::int_t> elemConnectivity(il::int_t k) {
    il::Array<il::int_t> temp(connectivity_.size(1));

    for (il::int_t i = 0; i < connectivity_.size(1); i++) {
      temp[i] = connectivity_(k, i);
    }

    return temp;
  };

  // todo: should be dofDD not dofDispl
  il::int_t connectivity(il::int_t k, il::int_t i) const {
    return connectivity_(k, i);
  }
  il::int_t dofPress(il::int_t k, il::int_t i) const {
    return dof_handle_pressure_(k, i);
  }
  il::int_t dofDispl(il::int_t k, il::int_t i) const {
    return dof_handle_displacement_(k, i);
  }
  il::int_t fracID(il::int_t k) const { return fracture_id_[k]; }
  il::int_t matID(il::int_t k) const { return material_id_[k]; }
  il::int_t condID(il::int_t k) const { return condition_id_[k]; }

  // il::int_t matid(il::int_t k) const { return material_id_[k]; }

  il::int_t nelts() const { return connectivity_.size(0); };

  il::int_t ncoor() const { return nodes_.size(0); };

  il::Array2D<double> coor() const { return nodes_; };

  il::Array2D<il::int_t> conn() const { return connectivity_; };

  il::Array<il::int_t> matid() const { return material_id_; };

// a method to output

  il::Array<il::int_t> eltinFrac(il::int_t k) ;

};
}

#endif  // HFPX2D_MESH_H
