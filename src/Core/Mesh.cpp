//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <src/Core/Mesh.h>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
//        CONSTRUCTORS/DESTRUCTORS
////////////////////////////////////////////////////////////////////////////////

il::Array2D<il::int_t> GetNodalEltConnectivity(
    const il::int_t nt_nodes, const il::Array2D<il::int_t> &connectivity) {
  // get element sharing common edges

  // loop over all nodes of the mesh
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

// function to build Tip nodes vector from nodal_connectivity table
il::Array<il::int_t> BuildTipNodes(
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

// function to build Tip element vector from nodal_connectivity table and the
// tipnodes array
il::Array<il::int_t> BuildTipElts(
    const il::Array2D<il::int_t> &node_connectivity,
    const il::Array<il::int_t> &tipnodes) {
  il::Array<il::int_t> tipelt(tipnodes.size());

  for (il::int_t i = 0; i < tipnodes.size(); i++) {
    tipelt[i] = node_connectivity(tipnodes[i], 0);
  }
  return tipelt;
}

Mesh::Mesh(const il::Array2D<double> &Coordinates,
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

  // matid_ was not passed as Input, so we assume the material is
  // homogeneous
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
  switch (interpolation_order_) {
    case 0: {
      il::Array2D<il::int_t> id_press{nelts, 1, 0};
      for (il::int_t e = 0; e < nelts; e++) {
        id_press(e, 0) = e;
      };
      dof_handle_pressure_ = id_press;
    }
    case 1:
      dof_handle_pressure_ = connectivity_;  // 1 unknowns per nodes ....
  };

  // build the nodal connected table...
  node_adj_elt_ = GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

  // built Tip nodes table...
  tipnodes_ = BuildTipNodes(node_adj_elt_);
  tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
};

Mesh::Mesh(const il::Array2D<double> &Coordinates,
           const il::Array2D<il::int_t> &Connectivity,
           const il::Array<int> &MatID,
           const il::int_t interpolationOrder) {
  // check validity of inputs

  IL_EXPECT_FAST(Coordinates.size(0) > 1 && Coordinates.size(1) == 2);
  IL_EXPECT_FAST(Connectivity.size(0) == MatID.size());

  // P0 and P1 elements only for now
  IL_EXPECT_FAST(interpolationOrder == 0 || interpolationOrder == 1);
  // check connectivity and coordinates consistency ??? currently no ->
  // they should be properly compatible

  coordinates_ = Coordinates;
  connectivity_ = Connectivity;
  interpolation_order_ = interpolationOrder;
  material_id_ = MatID;

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
  switch (interpolation_order_) {
    case 0: {
      il::Array2D<il::int_t> id_press{nelts, 1, 0};
      for (il::int_t e = 0; e < nelts; e++) {
        id_press(e, 0) = e;
      };
      dof_handle_pressure_ = id_press;
    }
    case 1:
      dof_handle_pressure_ = connectivity_;
  };

  // build the nodal connected table...
  node_adj_elt_ = GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

  // built Tip nodes table...
  tipnodes_ = BuildTipNodes(node_adj_elt_);
  tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
};

////////////////////////////////////////////////////////////////////////////////
//          METHODS
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

///
// get size of a all elements - output is an array of double
il::Array<double> Mesh::allEltSize() {
  il::Array<double> temp(connectivity_.size(0));
  for (il::int_t i = 0; i < numberOfElts(); i++) {
    temp[i] = eltSize(i);
  };
  return temp;
}

///
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

///
// method to get next to Tip elements - so-called ribbons element in an
// ILSA-like scheme

il::Array<il::int_t> Mesh::getRibbonElements() {
  il::Array<il::int_t> ribbon_elts{tipelts_.size(), 0};
  il::int_t nt_1;

  for (il::int_t e = 0; e < ribbon_elts.size(); e++) {
    // local connectivity of the Tip elements
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

///
// Function returning the segment characteristic object for element ne
hfp2d::SegmentData Mesh::getElementData(const il::int_t ne) {
  il::StaticArray2D<double, 2, 2> Xs;
  Xs(0, 0) = coordinates_(connectivity_(ne, 0), 0);
  Xs(0, 1) = coordinates_(connectivity_(ne, 0), 1);

  Xs(1, 0) = coordinates_(connectivity_(ne, 1), 0);
  Xs(1, 1) = coordinates_(connectivity_(ne, 1), 1);

  return SegmentData(Xs, interpolation_order_);
};

///
// adding elements function...
void Mesh::addNTipElements(const il::int_t t_e, const il::int_t the_tip_node,
                           const il::int_t n_add, double kink_angle) {
  // add n_add elements in the Mesh object ahead of the nodes the_tip_node
  // (global numbering) of
  // element t_e
  // with a kink_angle   with respect to element t_e
  //  The size of the added elements are equal to the size of element of t_e
  // the kick_angle should be given in the local Tip coordinates system.....
  //
  // api could be changed ...

  // we need the segment data .....
  hfp2d::SegmentData tipEltData = Mesh::getElementData(t_e);

  il::StaticArray<il::int_t, 2> tipEltConn = connectivity(t_e);
  double h = tipEltData.size(), global_prop_angle;

  il::int_t local_tip_node;
  il::StaticArray<double, 2> local_dir;

  // we need to know how if the Tip nodes is the first or second nodes of the
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
  il::Array2D<double> new_all_coor(numberOfNodes() + n_add, 2);

  for (il::int_t i = 0; i < numberOfNodes(); i++) {
    new_all_coor(i, 0) = coordinates(i, 0);
    new_all_coor(i, 1) = coordinates(i, 1);
  }
  //
  for (il::int_t i = 0; i < n_add; i++) {
    new_all_coor(i + numberOfNodes(), 0) = new_nodes(i, 0);
    new_all_coor(i + numberOfNodes(), 1) = new_nodes(i, 1);
  }

  // duplicating the old connectivity table....
  il::Array2D<il::int_t> new_conn(numberOfElts() + n_add, 2);
  for (il::int_t i = 0; i < numberOfElts(); i++) {
    new_conn(i, 0) = connectivity(i, 0);
    new_conn(i, 1) = connectivity(i, 1);
  }
  // adding the new connectivity at the end (so old element numbers are still
  // the same).
  for (il::int_t i = 0; i < n_add; i++) {
    if (i == 0) {
      new_conn(i + numberOfElts(), 0) = connectivity(t_e, local_tip_node);
      new_conn(i + numberOfElts(), 1) = numberOfNodes() + i;
    } else {
      new_conn(i + numberOfElts(), 0) = numberOfNodes() + i - 1;
      new_conn(i + numberOfElts(), 1) = numberOfNodes() + i;
    }
  };

  // now UPDATE THE MESH....
  coordinates_ = new_all_coor;
  connectivity_ = new_conn;

  //  the other changes....
  // this is for uniform material only
  il::Array<il::int_t> material_id_(connectivity_.size(0), 1);

  il::int_t nelts = connectivity_.size(0);
  il::int_t p = interpolation_order_;

  // COULD to be optimized below.... here we re-built everything from scratch...
  // anyway copy would be needed....

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
  switch (interpolation_order_) {
    case 0: {
      il::Array2D<il::int_t> id_press{nelts, 1, 0};
      for (il::int_t e = 0; e < nelts; e++) {
        id_press(e, 0) = e;
      };
      dof_handle_pressure_ = id_press;
    }
    case 1:
      dof_handle_pressure_ = connectivity_;  // 1 unknowns per nodes ....
  };

  // rebuild the nodal connected table...
  node_adj_elt_ = GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

  // rebuilt Tip nodes table...
  tipnodes_ = BuildTipNodes(node_adj_elt_);
  tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
}
}