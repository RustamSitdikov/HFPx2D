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

// Inclusion from standard library
#include <algorithm>
#include <cmath>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include <src/Core/SegmentData.h>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> coordinates_;

  // Connectivity matrix - size: number of elements x (order interpolation + 1)
  il::Array2D<il::int_t> connectivity_;

  // Interpolation order
  il::int_t interpolation_order_;

  // Dof handle matrices
  // for displacements  Discontinuities - size: number of elements x 2dofs per
  // coordinates x (order
  // interpolation + 1)
  il::Array2D<il::int_t> dof_handle_dd_;

  // for pressure - size: number of nodes x 1dof per coordinates x (order
  // interpolation  )
  il::Array2D<il::int_t> dof_handle_pressure_;

  // Identifier number of the fracture - size: number of elements
  //  il::Array<il::int_t> fracture_id_; //   not needed.....

  // Material identifier - size: number of elements
  il::Array<int> material_id_;

  // a structure  with nodes and corresponding adjacent elements .....
  // row node #, columms element sharing that nodes, if  entry is -1 then
  // no more connected elt.  todo: switch to a sparse matrix and use
  // smallArrays?
  il::Array2D<il::int_t> node_adj_elt_;

  //  2 arrays containing the tipnodes and the corresponding tipelts (could be a
  //  matrix)
  il::Array<il::int_t> tipnodes_;
  il::Array<il::int_t> tipelts_;

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////

 public:
  Mesh() = default;  // Default constructor
  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity,
       const il::int_t interpolationOrder);  // Constructor 1
  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity, const il::Array<int> &MatID,
       const il::int_t interpolationOrder);  // Constructor 2

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

 public:
  inline il::int_t numberOfElts() const { return connectivity_.size(0); };

  inline il::int_t numberOfNodes() const { return coordinates_.size(0); };

  // nodal coordinates related.
  inline il::Array2D<double> coordinates() const { return coordinates_; };
  // Read a particular element of the coordinates coordinates
  inline double coordinates(il::int_t k, il::int_t i) const {
    return coordinates_(k, i);
  }
  // Read the X coordinate of a coordinates
  inline double X(il::int_t k) const { return coordinates_(k, 0); }
  // Read the Y coordinate of a coordinates
  inline double Y(il::int_t k) const { return coordinates_(k, 1); }

  inline il::StaticArray<double, 2> coordinates(il::int_t k) const {
    il::StaticArray<double, 2> temp;
    temp[0] = coordinates_(k, 0);
    temp[1] = coordinates_(k, 1);
    return temp;
  };

  // connectivity related
  inline il::Array2D<il::int_t> connectivity() const { return connectivity_; };
  // get the connectivity of an element -> A StaticArray of size 2 here !
  il::StaticArray<il::int_t, 2> connectivity(il::int_t k) const {
    il::StaticArray<il::int_t, 2> temp;
    for (il::int_t i = 0; i < connectivity_.size(1); i++) {
      temp[i] = connectivity_(k, i);
    }
    return temp;
  };

  inline il::int_t connectivity(il::int_t e, il::int_t i) const {
    // element e, local coordinates i - return global nodes
    return connectivity_(e, i);
  }

  // nodal connectivity related
  inline il::Array2D<il::int_t> nodeEltConnectivity() const {
    return node_adj_elt_;
  };

  inline il::int_t nodeEltConnectivity(il::int_t k, il::int_t l) const {
    return node_adj_elt_(k, l);
  };

  inline il::Array<il::int_t> nodeEltConnectivity(il::int_t k) const {
    il::Array<il::int_t> temp(node_adj_elt_.size(1));
    for (il::int_t i = 0; i < node_adj_elt_.size(1); i++) {
      temp[i] = node_adj_elt_(k, i);
    }
    return temp;
  };

  // get Tip nodes
  inline il::Array<il::int_t> tipNodes() const { return tipnodes_; };
  inline il::int_t tip_nodes(il::int_t k) const { return tipnodes_[k]; };

  // get Tip elts
  inline il::Array<il::int_t> tipElts() const { return tipelts_; };
  inline il::int_t tip_elts(il::int_t k) const { return tipelts_[k]; };

  // material ID related
  inline il::Array<int> matid() const { return material_id_; };
  inline il::int_t matid(il::int_t k) const { return material_id_[k]; }

  inline il::int_t numberOfMaterials() const {
    return (*std::max_element(material_id_.begin(), material_id_.end()) + 1);
  }

  // interpolation order
  inline il::int_t interpolationOrder() const { return interpolation_order_; }

  // dofs related.....
  inline il::int_t numberDDDofsPerElt() const { return dof_handle_dd_.size(1); }

  inline il::int_t numberPressDofsPerElt() const {
    return dof_handle_pressure_.size(1);
  }

  inline il::int_t numberPressDofs() {
    il::int_t aux;
    switch (interpolation_order_) {
      case 0: {
        aux = dof_handle_pressure_.size(0);
      }
      case 1: {
        aux = coordinates_.size(0);
      }
    }
    return aux;
  }

  inline il::int_t numberDDDofs() const {
    return (numberOfElts() * (interpolation_order_ + 1) * 2);
  }

  inline il::int_t dofPress(il::int_t k, il::int_t i) const {
    // coordinates k, dof i -> return global equation iD
    return dof_handle_pressure_(k, i);
  }

  inline il::int_t dofDD(il::int_t k, il::int_t i) const {
    // coordinates k, dof i -> return global equation iD
    return dof_handle_dd_(k, i);  // element , dof dim.
  }

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  // a method to get the size of a given element.
  double eltSize(il::int_t &e);

  il::Array<double> allEltSize();

  // method to get nodes sharing 2 elts. (i.e. the nodal_connectivity for all
  // nodes with 2 neighbours)
  il::Array2D<il::int_t> getNodesSharing2Elts();

  // method to get the ribbon elements  - of a given mesh. (elements nearest to
  // a Tip element)
  il::Array<il::int_t> getRibbonElements();

  hfp2d::SegmentData getElementData(const il::int_t ne);

  // method to add N element ahead of a Tip node of a Tip element at a given
  // kick angle
  void addNTipElements(const il::int_t t_e, const il::int_t the_tip_node,
                       const il::int_t n_add, double kink_angle);
};
}

#endif  // HFPX2D_MESH_H
