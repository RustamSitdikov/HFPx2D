//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//
#include <il/linear_algebra/dense/norm.h>

#include <src/core/Mesh.h>
#include <src/core/SegmentData.h>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
//   FUNCTIONS REQUIRED IN THE CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
il::Array2D<il::int_t> getNodalEltConnectivity(
    const il::int_t nt_nodes, const il::Array2D<il::int_t> &connectivity) {
  // get element sharing common edges

  // loop over all nodes of the wellMesh
  // find corresponding elements
  // if n of elt sharing that coordinates ==2 -> put in array...
  //

  // maximum number of vertex with 2 elements is
  // this will not work for the case of more than 2 elements sharing the
  // coordinates.
  // we don t care of that case for now.

  // format is col1: element1, col2: element2 etc.   (note we don t store the
  // corresponding coordinates here....)

  // we should return a sparse matrix of integers....

  il::int_t n_elts = connectivity.size(0);

  il::Array<il::int_t> neConnec{nt_nodes};

  int j = 0;

  // finding  the number of elements that connect on node i
  for (il::int_t i = 0; i < nt_nodes; i++) {
    j = 0;
    for (il::int_t e = 0; e < n_elts; e++) {
      if (connectivity(e, 0) == i) {
        j++;
      } else if (connectivity(e, 1) == i) {
        j++;
      };
    };
    neConnec[i] = j;
  }

  // now get the nodes with the maximum of elements connected to it. // should
  // be a low number....
  il::int_t maxConnected =
      (*std::max_element(neConnec.begin(), neConnec.end()));

  // here we would like to have a sparse matrix in a sense...

  il::Array2D<il::int_t> node_connectivity{nt_nodes, maxConnected, -1};

  for (il::int_t i = 0; i < nt_nodes; i++) {
    j = 0;
    for (il::int_t e = 0; e < n_elts; e++) {
      if (connectivity(e, 0) == i) {
        node_connectivity(i, j) = e;
        j++;
      } else if (connectivity(e, 1) == i) {
        node_connectivity(i, j) = e;
        j++;
      };
      if (j == neConnec[i]) {
        break;
      }
    }
  }

  return node_connectivity;
}

////////////////////////////////////////////////////////////////////////////////
// function to build tip nodes vector from nodal_connectivity table
il::Array<il::int_t> buildTipNodes(
    const il::Array2D<il::int_t> &node_connectivity) {
  // find all nodes who do have only  1 adjacent element
  il::Array<il::int_t> temp(node_connectivity.size(0), 0);
  il::int_t k = 0;
  for (il::int_t i = 0; i < node_connectivity.size(0); i++) {
    if (node_connectivity(i, 1) == -1) {
      temp[k] = i;
      k++;
    }
  };

  il::Array<il::int_t> tipnodes(k, 0);
  for (il::int_t i = 0; i < k; i++) {
    tipnodes[i] = temp[i];
  }

  return tipnodes;
}
////////////////////////////////////////////////////////////////////////////////
// function to build tip element vector from nodal_connectivity table and the
// tipnodes array
il::Array<il::int_t> buildTipElts(
    const il::Array2D<il::int_t> &node_connectivity,
    const il::Array<il::int_t> &tipnodes) {
  il::Array<il::int_t> tipelt(tipnodes.size());

  for (il::int_t i = 0; i < tipnodes.size(); i++) {
    tipelt[i] = node_connectivity(tipnodes[i], 0);
  }
  return tipelt;
}
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////
//          METHODS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// get size of an element
double Mesh::eltSize(il::int_t &e) {
  il::StaticArray<double, 2> xdiff;
  xdiff[0] = coordinates_(connectivity_(e, 1), 0) -
             coordinates_(connectivity_(e, 0), 0);
  xdiff[1] = coordinates_(connectivity_(e, 1), 1) -
             coordinates_(connectivity_(e, 0), 1);
  double hx = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

  return hx;
};

////////////////////////////////////////////////////////////////////////////////
// get size of a all elements - output is an array of double
il::Array<double> Mesh::allEltSize() {
  il::Array<double> temp(connectivity_.size(0));
  for (il::int_t i = 0; i < numberOfElts(); i++) {
    temp[i] = eltSize(i);
  };
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
il::Array2D<il::int_t> Mesh::getNodesSharing2Elts() {
  // case of only 2 nodes for now.....
  // needed for building FD matrix  without fracture intersection.....
  // ONLY WORK ON MESH where nodes have maximum of 2 adjacaent elements
  // ..... no fracture intersection ...
  // todo generalize to account for fracture intersection

  IL_EXPECT_FAST(connectivity_.size(1) == 2);

  il::Array<il::int_t> tip = tipnodes_;

  il::Array2D<il::int_t> temp(numberOfNodes() - tip.size(), 2, 0);

  il::int_t k = 0;
  for (il::int_t i = 0; i < coordinates_.size(0); i++) {
    if (nodeEltConnectivity(i, 1) > -1) {
      temp(k, 0) = nodeEltConnectivity(i, 0);
      temp(k, 1) = nodeEltConnectivity(i, 1);
      k++;
    }
  }
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
// method to get next to tip elements - so-called ribbons element in an
// ILSA-like scheme
il::Array<il::int_t> Mesh::getRibbonElements() {
  il::Array<il::int_t> ribbon_elts{tipelts_.size(), 0};
  il::int_t nt_1;

  for (il::int_t e = 0; e < ribbon_elts.size(); e++) {
    // local connectivity of the tip elements
    if (connectivity(tipelts_[e], 0) == tipnodes_[e]) {
      nt_1 = connectivity(tipelts_[e], 1);
    } else {
      if (connectivity(tipelts_[e], 1) == tipnodes_[e]) {
        nt_1 = connectivity(tipelts_[e], 0);
      } else {
        il::abort();
      }
    }

    if (nodeEltConnectivity(nt_1, 0) == tipelts_[e]) {
      ribbon_elts[e] = nodeEltConnectivity(nt_1, 1);
    } else {
      if (nodeEltConnectivity(nt_1, 1) == tipelts_[e]) {
        ribbon_elts[e] = nodeEltConnectivity(nt_1, 0);
      } else {
        il::abort();
      }
    }
  }

  return ribbon_elts;
}

////////////////////////////////////////////////////////////////////////////////
// could do a method returning the element and nodes of a given fracture

////////////////////////////////////////////////////////////////////////////////
// Function returning the segment characteristic object for element ne
hfp2d::SegmentData Mesh::getElementData(const il::int_t ne) {
  il::StaticArray2D<double, 2, 2> Xs;
  Xs(0, 0) = coordinates_(connectivity_(ne, 0), 0);
  Xs(0, 1) = coordinates_(connectivity_(ne, 0), 1);

  Xs(1, 0) = coordinates_(connectivity_(ne, 1), 0);
  Xs(1, 1) = coordinates_(connectivity_(ne, 1), 1);

  return SegmentData(Xs, interpolation_order_);
};

////////////////////////////////////////////////////////////////////////////////
// adding elements metod...
void Mesh::addNTipElts(const il::int_t t_e, const il::int_t the_tip_node,
                       const il::int_t n_add, double kink_angle) {
  // add n_add elements in the Mesh object ahead of the nodes the_tip_node
  // (global numbering) of element t_e
  //  with a kink_angle   with respect to element t_e
  //  The size of the added elements are equal to the size of element of t_e
  //  the kick_angle should be given in the local tip coordinates system.....
  //
  //  API could be changed ...

  // we need the segment data .....
  hfp2d::SegmentData tipEltData = Mesh::getElementData(t_e);

  il::int_t fracNumber = fracid(t_e);

  il::StaticArray<il::int_t, 2> tipEltConn = connectivity(t_e);
  double h = tipEltData.size();
  double global_prop_angle = 0.;

  il::int_t local_tip_node = 0;
  il::StaticArray<double, 2> local_dir;

  // we need to know how if the tip nodes is the first or second nodes of the
  // element to know the direction of propagation

  if (the_tip_node == tipEltConn[0]) {
    local_tip_node = 0;
    global_prop_angle = tipEltData.theta() + il::pi - kink_angle;
  } else {
    if (the_tip_node == tipEltConn[1]) {
      local_tip_node = 1;
      global_prop_angle = tipEltData.theta() + kink_angle;
    } else {
      std::cout << "ERROR in AddNtipElements \n ";
      il::abort();  // must throw an exception here
    }
  };

  il::StaticArray<double, 2> prop_dir;  // in global coordinates....
  prop_dir[0] = h * cos(global_prop_angle);
  prop_dir[1] = h * sin(global_prop_angle);

  il::Array2D<double> new_nodes(n_add, 2);
  il::StaticArray<double, 2> tip_nodes_coor = coordinates(the_tip_node);

  // coordinates of the new nodes to be added
  for (il::int_t i = 0; i < n_add; i++) {
    if (i == 0) {
      new_nodes(i, 0) = tip_nodes_coor[0] + prop_dir[0];
      new_nodes(i, 1) = tip_nodes_coor[1] + prop_dir[1];
    } else {
      new_nodes(i, 0) = new_nodes(i - 1, 0) + prop_dir[0];
      new_nodes(i, 1) = new_nodes(i - 1, 1) + prop_dir[1];
    }
  }

  // reconstructing the whole coordinates array (sub-optimal)
  // could use resize()
  il::int_t nelts_old = numberOfElts();
  il::int_t nnodes_old = numberOfNodes();

  il::Array2D<double> new_all_coor(nnodes_old + n_add, 2);

  for (il::int_t i = 0; i < nnodes_old; i++) {
    new_all_coor(i, 0) = coordinates(i, 0);
    new_all_coor(i, 1) = coordinates(i, 1);
  }
  //
  for (il::int_t i = 0; i < n_add; i++) {
    new_all_coor(i + nnodes_old, 0) = new_nodes(i, 0);
    new_all_coor(i + nnodes_old, 1) = new_nodes(i, 1);
  }

  // duplicating the old connectivity table....
  il::Array2D<il::int_t> new_conn(nelts_old + n_add, 2);
  for (il::int_t i = 0; i < nelts_old; i++) {
    new_conn(i, 0) = connectivity(i, 0);
    new_conn(i, 1) = connectivity(i, 1);
  }
  // adding the new connectivity at the end (so old element numbers are still
  // the same).
  for (il::int_t i = 0; i < n_add; i++) {
    if (i == 0) {
      new_conn(i + nelts_old, 0) = connectivity(t_e, local_tip_node);
      new_conn(i + nelts_old, 1) = nnodes_old + i;
    } else {
      new_conn(i + nelts_old, 0) = nnodes_old + i - 1;
      new_conn(i + nelts_old, 1) = nnodes_old + i;
    }
  };

  il::int_t new_tip_elt = nelts_old+n_add - 1;
  il::int_t new_tip_node = nnodes_old+n_add -1 ;

  // now UPDATE THE MESH....
  coordinates_ = new_all_coor;
  connectivity_ = new_conn;

  il::int_t nelts = connectivity_.size(0);
  il::int_t p = interpolation_order_;

  // Material ID UPdate
  // this is for uniform material only
  // todo: handling heterogeneous media (mat_ID)

  material_id_.resize(nelts);
  for (il::int_t i = 0; i < n_add; i++) {
    material_id_[i + nelts_old] = 0;
  };

  // fracture id update...
  fracture_id_.resize(nelts);
  for (il::int_t i = 0; i < n_add; i++) {
    fracture_id_[i + nelts_old] = fracNumber;
  };

  // DOF HANDLES
  // COULD  be optimized below.... here we re-built everything from scratch...
  // anyway copy would be needed....
  /// Discontinuous Polynomial DOF handles
  il::Array2D<il::int_t> id_dd{nelts, 2 * (p + 1), 0};
  for (il::int_t i = 0; i < nelts; i++) {
    for (il::int_t j = 0; j < 2 * (p + 1); j++) {
      id_dd(i, j) = i * 2 * (p + 1) + j;
    }
  }
  dof_handle_dd_ = id_dd;  /// dof

  //   dof(element, local nnodes number)
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

  // rebuild the nodal connected table...
  node_adj_elt_ = getNodalEltConnectivity(coordinates_.size(0), connectivity_);

  // look

  // rebuilt tip nodes table...
  // find the index of this tip in the tipelts_
  //
  il::int_t j=0;
  for (il::int_t i=0;i<tipelts_.size();i++){
      if (tipelts_[i]==t_e){
        j=i;
        break;
      }
  }
  tipnodes_[j] = new_tip_node;
  tipelts_[j] = new_tip_elt;

//  tipnodes_ = buildTipNodes(node_adj_elt_);
//  tipelts_ = buildTipElts(node_adj_elt_, tipnodes_);
//
//  auto nfracs = numberOfFractures();
//  IL_EXPECT_FAST(tipnodes_.size()==nfracs*2);
//  IL_EXPECT_FAST(tipelts_.size()==nfracs*2);
//  //order them as function of the fractureID such that it is easier
//  // to track the tips
//  il::int_t  jf=0;
//  il::int_t  nl=0;// local node
//  il::Array<il::int_t> orderedTipNodes{2*nfracs,-1}; // init with -1
//  il::Array<il::int_t> orderedTipElts{2*nfracs,-1};
//  il::StaticArray<double,2> tipcoor1;
//  il::StaticArray<double,2> tipcoor2;
//  double dist1=0; double dist2=0.;
//  il::int_t iaux; il::int_t nl2=1;
//  for (il::int_t i=0;i<tipnodes_.size();i++ ){
//    tipcoor1=coordinates(tipnodes_[i]);
//    dist1 = il::norm(tipcoor1,il::Norm::L2);
//    jf=fracture_id_[tipelts_[i]];
//    nl=0;
//    if ((orderedTipElts[jf*2]!=-1))
//    {  // already a node stored
//      tipcoor2=coordinates(orderedTipElts[jf*2]);
//      dist2 = il::norm(tipcoor2,il::Norm::L2);
//      if (tipcoor1[0]>tipcoor2[0]){
//        nl=1;
//      } else {
//        if (tipcoor1[0]==tipcoor2[0]){
//          if (tipcoor1[1]>tipcoor2[1]){
//            nl=1;
//          }
//        }
//      }
//      if (nl==0) {
//        iaux = orderedTipElts[jf*2];
//        orderedTipElts[jf*2+1]=iaux;
//        iaux = orderedTipNodes[jf*2];
//        orderedTipNodes[jf*2+1]=iaux;
//      }
//    };
//
//    orderedTipElts[jf*2+nl]=tipelts_[i];
//    orderedTipNodes[jf*2+nl]=tipnodes_[i];
//  }
//  tipnodes_=orderedTipNodes;
//  tipelts_=orderedTipElts;
}
}