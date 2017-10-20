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

#include <src/core/SegmentData.h>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
//   FUNCTIONS REQUIRED IN THE CONSTRUCTOR
////////////////////////////////////////////////////////////////////////////////

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

  // finding for the number of elements that connect on node i
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
  // be a low number.
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
// build tip  nodes vector
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
////////////////////////////////////////////////////////////////////////////////
// build tip element vector
il::Array<il::int_t> BuildTipElts(
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

////////////////////////////////////////////////////////////////////////////////
///  METHODS

// appendMesh is not tested and NEEDS TO BE MODIFIED following the other
// changes...

void Mesh::appendMesh(const Mesh &newMesh, const bool isJoined) {
  // IL_EXPECT_FAST(interpolation_order_==newMesh.interpolation_order_);

  // Strategy for appendMesh method:
  // newMesh arrives and it can be either connected or not.
  // isJoined determines if the mesh will share some nodes with the previous
  // one.
  //

  if (isJoined) {
    /* PLACEHOLDER
      // If isJoined == true, then the two meshes are connected VIA THE TIPS.
      // The is_tip vector determines which nodes in the mesh are tips. The
      following
      // step requires to check which nodes are in contacts, by computing the
      distance
      // of the two tip nodes and finding the pair with the minimum distance

      il::int_t connectedOldNode;
      il::int_t connectedNewNode;

      double distTol = 10 ^-6;

      for (il::int_t i = 0; i < number_nodes_; i++) {

        if (is_tip_[i]) {

          for (il::int_t j = 0; j < newMesh.number_nodes_; j++) {

            if (newMesh.is_tip_[j]) {

              double distance =
                  sqrt(pow((coordinates_(i, 0) - newMesh.coordinates_(j, 0)), 2)
      +
      pow((coordinates_(i, 1) - newMesh.coordinates_(j, 1)), 2));

              if (distance < distTol) {
                connectedOldNode = i;
                connectedNewNode = j;
              }
            }
          }
        }
      }

      // Now, the number of nodes are summed -1 (one node must be collapsed)
      // The number of elements is summed.
      // Nodes coordinates array is modified as following. The node to be
      collapsed is
      // stored in the variable connectedOldNode. That is maintained.
      // In the new mesh, the node that is collapsed is called connectedNewNode.
      That is
      // removed from the list and the other nodes are appended to the old list
      of
      // nodal coordinates. Then, in the connectivity matrix, the nodes with
      number less
      // than connectedNewNode are increased by oldNumNodes. connectedNewNode
      itself is
      // substituted by connectedOldNodes and the nodes with number greater than
      // connectedNewNodes are increased by oldNumNodes-1
      //
      //
      // 0,1,2,3....15,16 old nodes \
        //    common node is 13 and 7 ----> 0-16 (17 total old nodes),
      0+17,1+17,...6+17,8+17-1,...11+17-1,12+17-1
      // 0,1,2,.....11,12 new nodes /
      //
      //
      // The node number 7 (connected new node) in the new mesh is omitted and
      substituted
      // in the connectivity by the number 13 (connected old node).
      // Handle matrices for displacements are appended increased by
      oldNumElements*2*(interpolation+1)
      // Handle matrices for pressure are appended with NumNodes offset before
      the connectedNewNode.
      // Connected new coordinates value is substituted with the old connected
      value.
      The subsequent nodes
      // are increased by the old number of nodes -1
      //
      //
      //
      //
      //


      //

      const il::int_t oldNumNodes = number_nodes_;
      const il::int_t oldNumElements = number_elements_;

      number_nodes_ = number_nodes_ + newMesh.number_nodes_;
      number_elements_ = number_elements_ + newMesh.number_elements_;

      coordinates_.resize(number_nodes_, 2); // here the 2 is because we are
      working
      in 2D
      for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
        coordinates_(i, 0) = newMesh.coordinates_(i, 0);
        coordinates_(i, 1) = newMesh.coordinates_(i, 1);
      }

      ///// APPENDING CONNECTIVITY AND DOF HANDLES
      // se dobbiamo aggiungere una mesh che é collegata a quella giá salvata
      // possiamo:
      // - per i dof_handle_displacement basta concatenare i dofs aggiuntivi
      //   visto che gli elementi per gli spostamenti sono indipendenti
      // - per i dof_handle_pressure, bisogna avere lo stesso valore del nodo
      //   libero. quindi bisogna sapere da quale nodo ripartire
      // - per la connettivita' basta aggiungere i nuovi nodi con una nuova
      //   numerazione (numerazione locale + shift degli elementi che sono già
      //   inseriti nella mesh => number_elements)




      // connectivity resize: here the 2 is because the element is 1D with 2
      nodes
      // but it should rather be the order of the polynomia on the element plus
      1
      connectivity_.resize(number_elements_, 2);
      material_id_.resize(number_elements_);
      fracture_id_.resize(number_elements_);

      // add the values of the new mesh
      for (il::int_t i = oldNumElements; i < number_elements_; i++) {

        connectivity_(i, 0) = newMesh.connectivity_(i, 0);
        connectivity_(i, 1) = newMesh.connectivity_(i, 1);

        material_id_[i] = newMesh.material_id_[i];
        fracture_id_[i] = newMesh.fracture_id_[i];

      }*/

  } else {
    // If isJoined == false, then the two meshes are not connected.
    // Hence, the total number of nodes and elements are summed together.
    // Connectivity matrices are concatenated. The new connectivity matrix
    // contains
    // values from 0 to the new number of elements. Consequently, all the
    // connectivity
    // will be shifted by the old number of elements.
    // Similarly, dof_handle matrices are shifted respectively by old number of
    // nodes
    // (in the case of the pressure one) and
    // num_elements*2*(interpolation_order+1)
    // (in the case of the displacement one).
    // Fracture ID and Material ID are saved as such.

    // Computing the new number of nodes (being 2 independent fractures, they
    // are just summed)
    const il::int_t old_number_of_nodes = this->numberOfNodes();
    const il::int_t new_number_of_nodes =
        this->numberOfNodes() + newMesh.numberOfNodes();

    // Resize the array of nodes coordinates accordingly and fill it with the
    // new values
    coordinates_.resize(new_number_of_nodes, 2);

    for (il::int_t i = 0; i < newMesh.numberOfNodes(); i++) {
      coordinates_(i + old_number_of_nodes, 0) = newMesh.coordinates_(i, 0);
      coordinates_(i + old_number_of_nodes, 1) = newMesh.coordinates_(i, 1);
    }

    // Compute the new number of elements (simply summed as well)
    const il::int_t old_number_of_elements = this->numberOfElements();
    const il::int_t new_number_of_elements =
        this->numberOfElements() + newMesh.numberOfElements();

    // Resize connectivity accordingly and save the connectivity of new
    // elements, which nodes identifier
    // have with a shift of old_number_of_elements w.r.t. the newMesh values

    // also materialID and fractureID are element related and concatenated here

    const il::int_t columns_conn_mtx =
        (interpolation_order_ == 0) ? 2 : (interpolation_order_ + 1);

    connectivity_.resize(new_number_of_elements, columns_conn_mtx);

    //    fracture_id_.resize(new_number_of_elements);

    material_id_.resize(new_number_of_elements);

    for (il::int_t i = 0; i < newMesh.numberOfElements(); i++) {
      for (il::int_t j = 0; j < columns_conn_mtx; j++) {
        connectivity_(i + old_number_of_elements, j) =
            newMesh.connectivity_(i, j) + old_number_of_nodes;
      }

      material_id_[i + old_number_of_elements] = newMesh.material_id_[i];
      //      fracture_id_[i + old_number_of_elements] =
      //      newMesh.fracture_id_[i];
    }

    std::cout << "got connectivity" << std::endl;

    // Displacement dof handles: concatenated and shifted by old_number_of_nodes
    // * 2 * (interpolationOrder + 1)
    const il::int_t columns_displ_dofh_mtx =
        ((interpolation_order_ == 0) ? (2 * 2)
                                     : (2 * (interpolation_order_ + 1)));
    const il::int_t displ_dof_handle_shift =
        old_number_of_elements * columns_displ_dofh_mtx;

    dof_handle_dd_.resize(new_number_of_elements, columns_displ_dofh_mtx);

    for (il::int_t i = 0; i < newMesh.numberOfElements(); i++) {
      for (il::int_t j = 0; j < columns_displ_dofh_mtx; j++) {
        dof_handle_dd_(i + old_number_of_elements, j) =
            newMesh.dof_handle_dd_(i, j) + displ_dof_handle_shift;
      }
    }

    std::cout << "got displ dof" << std::endl;

    // Pressure dof handles are concatenated and shifted by
    // old_number_of_elements + 1
    const il::int_t columns_press_dofh_mtx =
        (interpolation_order_ == 0) ? 2 : (interpolation_order_ + 1);

    // WRONG below
    const il::int_t press_dof_handle_shift =
        old_number_of_elements * (columns_press_dofh_mtx - 1) + 1;
    //        numberOfFractures();

    dof_handle_pressure_.resize(new_number_of_elements, columns_press_dofh_mtx);

    for (il::int_t i = 0; i < newMesh.numberOfElements(); i++) {
      for (il::int_t j = 0; j < columns_press_dofh_mtx; j++) {
        dof_handle_pressure_(i + old_number_of_elements, j) =
            newMesh.dof_handle_pressure_(i, j) + press_dof_handle_shift;
      }
    }

    std::cout << "got press dof" << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
//          METHODS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// get size of an element
double Mesh::elt_size(il::int_t &e) {
  il::StaticArray<double, 2> xdiff;
  xdiff[0] = coordinates_(connectivity_(e, 1), 0) -
             coordinates_(connectivity_(e, 0), 0);
  xdiff[1] = coordinates_(connectivity_(e, 1), 1) -
             coordinates_(connectivity_(e, 0), 1);
  double hx = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

  return hx;
};

////////////////////////////////////////////////////////////////////////////////
// get size of a all elements
il::Array<double> Mesh::All_elt_size() {
  il::Array<double> temp(connectivity_.size(0));
  for (il::int_t i = 0; i < numberOfElements(); i++) {
    temp[i] = elt_size(i);
  };
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
il::Array2D<il::int_t> Mesh::GetNodesSharing2Elts() {
  // needed for building FD matrix  without fracture intersection.....
  // ONLY WORK ON MESH where nodes have maximum of 2 adjacaent elements
  // ..... no fracture intersection .....

  IL_EXPECT_FAST(connectivity_.size(1) == 2);

  il::Array<il::int_t> tip = tipnodes_;

  il::Array2D<il::int_t> temp(numberOfNodes() - tip.size(), 2, 0);

  il::int_t k = 0;
  for (il::int_t i = 0; i < coordinates_.size(0); i++) {
    if (node_elt_connectivity(i, 1) > -1) {
      temp(k, 0) = node_elt_connectivity(i, 0);
      temp(k, 1) = node_elt_connectivity(i, 1);
      k++;
    }
  }
  // case of only 2 nodes for now.....
  return temp;
}

////////////////////////////////////////////////////////////////////////////////
// could do a method returning the element and nodes of a given fracture
/// this is honestly not really needed.....

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
// add element function...
void Mesh::AddNTipElements(const il::int_t t_e, const il::int_t the_tip_node,
                           const il::int_t n_add, double kink_angle) {
  // add n_add elements in the Mesh object ahead of the nodes the_tip_node of
  // element t_w
  //  with size
  // equal to element of t_e
  // with a kink_angle   with respect to element t_e
  // kick_angle in local tip coordinates system.....
  // api could be changed....

  // we need the segment data .....
  hfp2d::SegmentData tipEltData = Mesh::getElementData(t_e);

  il::StaticArray<il::int_t, 2> tipEltConn = connectivity(t_e);
  double h = tipEltData.size(), global_prop_angle;

  il::int_t local_tip_node;
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

  for (il::int_t i = 0; i < n_add; i++) {
    if (i == 0) {
      new_nodes(i, 0) = tip_nodes_coor[0] + prop_dir[0];
      new_nodes(i, 1) = tip_nodes_coor[1] + prop_dir[1];
    } else {
      new_nodes(i, 0) = new_nodes(i - 1, 0) + prop_dir[0];
      new_nodes(i, 1) = new_nodes(i - 1, 1) + prop_dir[1];
    }
  }

  il::Array2D<double> new_all_coor(numberOfNodes() + n_add, 2);

  for (il::int_t i = 0; i < numberOfNodes(); i++) {
    new_all_coor(i, 0) = coordinates(i, 0);
    new_all_coor(i, 1) = coordinates(i, 1);
  }
  for (il::int_t i = 0; i < n_add; i++) {
    new_all_coor(i + numberOfNodes(), 0) = new_nodes(i, 0);
    new_all_coor(i + numberOfNodes(), 1) = new_nodes(i, 1);
  }

  il::Array2D<il::int_t> new_conn(numberOfElements() + n_add, 2);
  for (il::int_t i = 0; i < numberOfElements(); i++) {
    new_conn(i, 0) = connectivity(i, 0);
    new_conn(i, 1) = connectivity(i, 1);
  }
  for (il::int_t i = 0; i < n_add; i++) {
    if (i == 0) {
      new_conn(i + numberOfElements(), 0) = connectivity(t_e, local_tip_node);
      new_conn(i + numberOfElements(), 1) = numberOfNodes() + i;
    } else {
      new_conn(i + numberOfElements(), 0) = numberOfNodes() + i - 1;
      new_conn(i + numberOfElements(), 1) = numberOfNodes() + i;
    }
  };

  coordinates_ = new_all_coor;
  connectivity_ = new_conn;

  //  the other changes....
  // this is for uniform material only
  il::Array<il::int_t> material_id_(connectivity_.size(0), 1);
  il::Array<il::int_t> condition_id_(connectivity_.size(0), 1);

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

  // build the nodal connected table...
  node_adj_elt_ = GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

  // built tip nodes table...
  tipnodes_ = BuildTipNodes(node_adj_elt_);
  tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
}


}