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

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linear_algebra/dense/norm.h>

// Inclusion from hfp2d
#include <src/core/SegmentData.h>

namespace hfp2d {

il::Array2D<il::int_t> getNodalEltConnectivity(
    const il::int_t nt_nodes, const il::Array2D<il::int_t> &connectivity);

il::Array<il::int_t> buildTipNodes(
    const il::Array2D<il::int_t> &node_connectivity);

il::Array<il::int_t> buildTipElts(
    const il::Array2D<il::int_t> &node_connectivity,
    const il::Array<il::int_t> &tipnodes);

///// 1D mesh class
class Mesh {  // class for 1D wellMesh of 1D segment elements ?

 private:
  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> coordinates_;

  // Connectivity matrix - size: number of elements x  2
  il::Array2D<il::int_t> connectivity_;

  // Interpolation order
  il::int_t interpolation_order_;

  // Dof handle matrices
  // for displacements  Discontinuities - size: number of elements x 2 dofs per
  // coordinates x (order interpolation + 1)
  il::Array2D<il::int_t> dof_handle_dd_;

  // for pressure - size: number of nodes x 1dof per coordinates x (order
  // interpolation  )
  il::Array2D<il::int_t> dof_handle_pressure_;

  // Identifier number of the fracture - size: number of elements
  il::Array<il::int_t> fracture_id_;  //

  // Material identifier - size: number of elements
  il::Array<il::int_t> material_id_;

  // a structure  with nodes and corresponding adjacent elements .....
  // row node #, columms element sharing that nodes, if  entry is -1 then
  // no more connected elt.  todo: switch to a sparse matrix and use
  // smallArrays?
  il::Array2D<il::int_t> node_adj_elt_;

  //  2 arrays containing the tipnodes and the corresponding tipelts (could be a
  //  matrix)
  il::Array<il::int_t> tipnodes_;
  il::Array<il::int_t> tipelts_;

 public:
  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////////

  // todo: naming of the different entities are not consistent AND TOO LONG

  //   Mesh()default;
  Mesh(){};

  // Basic constructor with  coordinates and connectivity array and
  // interpolation order
  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity,
       const il::int_t interpolationOrder) {
    // check validity of inputs

    IL_EXPECT_FAST(Coordinates.size(0) > 1 && Coordinates.size(1) == 2);
    // P0 and P1 elements only for now
    IL_EXPECT_FAST(interpolationOrder == 0 || interpolationOrder == 1);
    // check connectivity and coordinates consistency ??? currently no ->
    // they should be properly compatible

    coordinates_ = Coordinates;
    connectivity_ = Connectivity;
    interpolation_order_ = interpolationOrder;

    // matid_ was not passed as input, so we assume the material is homogeneous
    il::Array<il::int_t> material_id_(connectivity_.size(0), 1);

    il::int_t nelts = connectivity_.size(0);
    il::int_t p = interpolation_order_;

    /// Discontinuous Polynomial DOF handles
    il::Array2D<il::int_t> id_dd{nelts, 2 * (p + 1), 0};
    for (il::int_t i = 0; i < nelts; i++) {
      for (il::int_t j = 0; j < 2 * (p + 1); j++) {
        id_dd(i, j) = i * 2 * (p + 1) + j;
      }
    }
    dof_handle_dd_ = id_dd;  /// dof

    /// //    dof(element, local nnodes number)
    // actually this is the connectivity_ array for  p =1 and
    // a simple elt number of P0
    if (interpolation_order_ == 0) {
      il::Array2D<il::int_t> id_press{nelts, 1};
      for (il::int_t e = 0; e < nelts; e++) {
        id_press(e, 0) = e;
      };
      dof_handle_pressure_ = id_press;
    } else {  // case 1
      dof_handle_pressure_ = connectivity_;
    };

    // build the nodal connected table...
    node_adj_elt_ =
        getNodalEltConnectivity(coordinates_.size(0), connectivity_);

    // built tip nodes table...
    tipnodes_ = buildTipNodes(node_adj_elt_);
    tipelts_ = buildTipElts(node_adj_elt_, tipnodes_);
    // built fracture_id_ ?...
    // default just one fracture if not passed ?-> todo generalize the
    // construction to automatic detection of # fracture ?
    il::Array<il::int_t> fractureID{nelts, 0};
    fracture_id_ = fractureID;
  };

  // case where the matid vector and the fracture_ID vector are provided
  // constructor with interpolation order and coordinates and connectivity array
  Mesh(const il::Array2D<double> &Coordinates,
       const il::Array2D<il::int_t> &Connectivity,
       const il::Array<il::int_t> &MatID,
       const il::Array<il::int_t> &FractureID,
       const il::int_t interpolationOrder) {
    // check validity of inputs

    IL_EXPECT_FAST(Coordinates.size(0) > 1 && Coordinates.size(1) == 2);
    IL_EXPECT_FAST(Connectivity.size(0) == MatID.size());
    IL_EXPECT_FAST(Connectivity.size(0) == FractureID.size());

    // P0 and P1 elements only for now
    IL_EXPECT_FAST(interpolationOrder == 0 || interpolationOrder == 1);
    // check connectivity and coordinates consistency ??? currently no ->
    // they should be properly compatible

    coordinates_ = Coordinates;
    connectivity_ = Connectivity;
    interpolation_order_ = interpolationOrder;

    // check that the entries of FractureID are within a range 0 to N by step
    // of one without blank
    il::Array<il::int_t> checkMat = MatID;
    std::sort(checkMat.begin(), checkMat.end());
    auto last = std::unique(checkMat.begin(), checkMat.end());
    IL_EXPECT_FAST(checkMat[0] == 0);
    il::int_t nmat = 0;
    for (il::int_t i = 0; i < checkMat.size(); i++) {
      if (checkMat[i] == i) {
        nmat = nmat + 1;
      } else {
        break;
      }
    }
    auto nmat2 = (*std::max_element(MatID.begin(), MatID.end()) + 1);
    IL_EXPECT_FAST(nmat == nmat2);
    material_id_ = MatID;

    // check that the entries of FractureID are within a range0 to N by step
    // of one without blank
    il::Array<il::int_t> checkFrac = FractureID;
    std::sort(checkFrac.begin(), checkFrac.end());
    last = std::unique(checkFrac.begin(), checkFrac.end());
    IL_EXPECT_FAST(checkFrac[0] == 0);
    il::int_t nfracs = 0;
    for (il::int_t i = 0; i < checkFrac.size(); i++) {
      if (checkFrac[i] == i) {
        nfracs = nfracs + 1;
      } else {
        break;
      }
    }
    auto nfrac2 = (*std::max_element(FractureID.begin(), FractureID.end()) + 1);
    IL_EXPECT_FAST(nfracs == nfrac2);
    fracture_id_ = FractureID;

    il::int_t nelts = connectivity_.size(0);
    il::int_t p = interpolation_order_;

    // Discontinuous Polynomial DOF handle for DD
    il::Array2D<il::int_t> id_dd{nelts, 2 * (p + 1), 0};
    for (il::int_t i = 0; i < nelts; i++) {
      for (il::int_t j = 0; j < 2 * (p + 1); j++) {
        id_dd(i, j) = i * 2 * (p + 1) + j;
      }
    }
    dof_handle_dd_ = id_dd;  /// dof

    ///  Continous scalar polynomial for Pressure
    // actually this is the connectivity_ array for  p =1 and
    // a simple elt number of P0
    if (interpolation_order_ == 0) {
      il::Array2D<il::int_t> id_press{nelts, 1};
      for (il::int_t e = 0; e < nelts; e++) {
        id_press(e, 0) = e;
      };
      dof_handle_pressure_ = id_press;
    } else {  // case 1
      dof_handle_pressure_ = connectivity_;
    };

    // build the nodal connected table...
    node_adj_elt_ =
        getNodalEltConnectivity(coordinates_.size(0), connectivity_);

    // built tip nodes table...
    tipnodes_ = buildTipNodes(node_adj_elt_);
    tipelts_ = buildTipElts(node_adj_elt_, tipnodes_);
    IL_EXPECT_FAST(tipnodes_.size()==nfracs*2);
    IL_EXPECT_FAST(tipelts_.size()==nfracs*2);

    //order them as function of the fractureID such that it is easier
    // to track the tips
    il::int_t  jf=0;
    il::int_t  nl=0;// local node
    il::Array<il::int_t> orderedTipNodes{2*nfracs,-1}; // init with -1
    il::Array<il::int_t> orderedTipElts{2*nfracs,-1};
    il::StaticArray<double,2> tipcoor1;
    il::StaticArray<double,2> tipcoor2;
    double dist1=0; double dist2=0.;
    il::int_t iaux; il::int_t nl2=1;
    for (il::int_t i=0;i<tipnodes_.size();i++ ){
      tipcoor1=coordinates(tipnodes_[i]);
      dist1 = il::norm(tipcoor1,il::Norm::L2);
      jf=fracture_id_[tipelts_[i]];
      nl=0;
      if ((orderedTipElts[jf*2]!=-1))
      {  // already a node stored
        tipcoor2=coordinates(orderedTipElts[jf*2]);
        dist2 = il::norm(tipcoor2,il::Norm::L2);
        if (dist1>dist2){
          nl=1;
        } else {
          if (dist1==dist2){
            if (tipcoor1[0]>tipcoor2[0]){
              nl=1;
            }
          }
        }
        if (nl==0) {
          iaux = orderedTipElts[jf*2];
          orderedTipElts[jf*2+1]=iaux;
          iaux = orderedTipNodes[jf*2];
          orderedTipNodes[jf*2+1]=iaux;
        }
      };

      orderedTipElts[jf*2+nl]=tipelts_[i];
      orderedTipNodes[jf*2+nl]=tipnodes_[i];
    }
    tipnodes_=orderedTipNodes;
    tipelts_=orderedTipElts;
  };

  //////////////////////////////////////////////////////////////////////////
  //        get functions  - i.e. public interfaces
  //////////////////////////////////////////////////////////////////////////

  il::int_t numberOfElts() const { return connectivity_.size(0); };

  il::int_t numberOfNodes() const { return coordinates_.size(0); };

  // nodal coordinates related.
  il::Array2D<double> coordinates() const { return coordinates_; };

  // Read a particular element of the coordinates coordinates
  double coordinates(il::int_t k, il::int_t i) const {
    return coordinates_(k, i);
  }

  // Read the X coordinate of a coordinates
  double X(il::int_t k) const { return coordinates_(k, 0); }
  // Read the Y coordinate of a coordinates
  double Y(il::int_t k) const { return coordinates_(k, 1); }

  il::StaticArray<double, 2> coordinates(il::int_t k) const {
    il::StaticArray<double, 2> temp;
    temp[0] = coordinates_(k, 0);
    temp[1] = coordinates_(k, 1);
    return temp;
  };

  // connectivity related
  il::Array2D<il::int_t> connectivity() const { return connectivity_; };
  // get the connectivity of an element -> A StaticArray of size 2 here !
  il::StaticArray<il::int_t, 2> connectivity(il::int_t k) const {
    il::StaticArray<il::int_t, 2> temp;
    for (il::int_t i = 0; i < connectivity_.size(1); i++) {
      temp[i] = connectivity_(k, i);
    }
    return temp;
  };

  //
  il::int_t connectivity(il::int_t e, il::int_t i) const {
    // element e, local coordinates i - return global nodes
    return connectivity_(e, i);
  }

  // nodal connectivity related
  il::Array2D<il::int_t> nodeEltConnectivity() const { return node_adj_elt_; };
  il::int_t nodeEltConnectivity(il::int_t k, il::int_t l) const {
    return node_adj_elt_(k, l);
  };

  il::Array<il::int_t> nodeEltConnectivity(il::int_t k) const {
    il::Array<il::int_t> temp(node_adj_elt_.size(1));
    for (il::int_t i = 0; i < node_adj_elt_.size(1); i++) {
      temp[i] = node_adj_elt_(k, i);
    }
    return temp;
  };

  // get tip nodes
  il::Array<il::int_t> tipNodes() const { return tipnodes_; };
  il::int_t tipNodes(il::int_t k) const { return tipnodes_[k]; };

  // get tip elts
  il::Array<il::int_t> tipElts() const { return tipelts_; };
  il::int_t tipElts(il::int_t k) const { return tipelts_[k]; };

  // material ID related
  il::Array<il::int_t> matid() const { return material_id_; };
  il::int_t matid(il::int_t k) const { return material_id_[k]; }

  // todo check that material_id_ and frature_id_ contains entries in the whole
  // range 0 to N (no blank)
  il::int_t numberOfMaterials() const {  // assume entries from 0 to N
    return (*std::max_element(material_id_.begin(), material_id_.end()) + 1);
  }

  // fracture ID related
  il::Array<il::int_t> fracid() const { return fracture_id_; };
  il::int_t fracid(il::int_t k) const { return fracture_id_[k]; }

  il::int_t numberOfFractures() const {  // assume entries from 0 to N
    return (*std::max_element(fracture_id_.begin(), fracture_id_.end()) + 1);
  }

  // interpolation order
  il::int_t interpolationOrder() const { return interpolation_order_; }

  // dofs related.....
  il::int_t numberDDDofsPerElt() const { return dof_handle_dd_.size(1); }

  il::int_t numberPressDofsPerElt() const {
    return dof_handle_pressure_.size(1);
  }

  il::int_t numberPressDofs() const {
    il::int_t aux = 0;
    aux = dof_handle_pressure_.size(0);  // p0
    if (interpolation_order_ == 1) {     // p1
      aux = coordinates_.size(0);
    }
    return aux;
  }

  il::int_t numberDDDofs() const {
    return (numberOfElts() * (interpolation_order_ + 1) * 2);
  }

  il::int_t dofPress(il::int_t k, il::int_t i) const {
    // coordinates k, dof i -> return global equation iD
    return dof_handle_pressure_(k, i);
  }

  il::int_t dofDD(il::int_t k, il::int_t i) const {
    // coordinates k, dof i -> return global equation iD
    return dof_handle_dd_(k, i);  // element , dof dim.
  }

  ////////////////////////////////////////////////////////////////////////////////////////////
  //   Methods
  ////////////////////////////////////////////////////////////////////////////////////////////

  hfp2d::SegmentData getElementData(il::int_t ne);

  // a method to get the size of a given element.
  double eltSize(il::int_t &e);

  il::Array<double> allEltSize();

  // method to get nodes sharing 2 elts. (i.e. the nodal_connectivity for all
  // nodes with 2 neighbours)
  // todo rename
  il::Array2D<il::int_t> getNodesSharing2Elts();

  // method to get the ribbon elements  - of a given wellMesh. (elements nearest
  // to a tip element)
  il::Array<il::int_t> getRibbonElements();

  // method to add N element ahead of a tip node of a tip element at a given
  // kick angle
  void addNTipElts(il::int_t t_e, il::int_t the_tip_node, il::int_t n_add,
                   double kink_angle);

  void reorderTip() {  // todo
    il::Array<il::int_t> tip_fid{tipelts_.size(), 0};
    for (il::int_t i = 0; i < tipelts_.size(); i++) {
      tip_fid[i] = fracture_id_[tipelts_[i]];
    }
  };
};
}

#endif  // HFPX2D_MESH_H
