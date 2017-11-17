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
// function to build tip element vector from nodal_connectivity table and the
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
////////////////////////////////////////////////////////////////////////////////




//////////////////////////////// CONSTRUCTORS (a.k.a. initializers) ////////////////////////////////
// Complete initialization of the mesh class
//Mesh::Mesh(const il::int_t interpolationOrder,
//           const il::Array2D<double> &nodesCoordinates,
//           const il::Array2D<il::int_t> &elementsConnectivity,
//           const il::Array2D<il::int_t> &displ_dof_handle,
//           const il::Array2D<il::int_t> &press_dof_handle,
//           const il::Array<il::int_t> &fractureID,
//           const il::Array<il::int_t> &materialID,
//           const il::Array<il::int_t> &conditionID) {
//
//  // Initial assertions to check data consistency
//  // Non-zero number of nodes and coordinates must be 2D
//  IL_EXPECT_FAST(nodesCoordinates.size(0) > 0 && nodesCoordinates.size(1) == 2);
//
//  // Non-zero number of elements and connectivity matrix shall have
//  // as many columns as the interpolation order +1
//  IL_EXPECT_FAST(elementsConnectivity.size(0) > 0);
//  IL_EXPECT_FAST((interpolationOrder==0 && elementsConnectivity.size(1)==2)
//                      || (interpolationOrder>0  && elementsConnectivity.size
//                          (1)==interpolationOrder + 1) );
//
//  // Check of dof handles sizes
//  // - same number of elements for pressure and displacements
//  IL_EXPECT_FAST(elementsConnectivity.size(0)==displ_dof_handle.size(0));
//  IL_EXPECT_FAST(elementsConnectivity.size(0)==press_dof_handle.size(0));
//  // - 2 displacement dofs per node x number of nodes in an element == dof handle per element
//  IL_EXPECT_FAST((interpolationOrder==0 && displ_dof_handle.size(1)==2) ||
//                     elementsConnectivity.size(1)*2 == displ_dof_handle.size(1));
//  // - 1 pressure dof per node x number of nodes in an element == dof handle per element
//  IL_EXPECT_FAST((interpolationOrder==0 && press_dof_handle.size(1)==2) ||
//                     (interpolationOrder + 1 == press_dof_handle.size(1)));
//
//  // Check of the size of fractureID and materialID
//  IL_EXPECT_FAST(elementsConnectivity.size(0) == fractureID.size());
//  IL_EXPECT_FAST(elementsConnectivity.size(0) == materialID.size());
//
//  // Assignment to the class members
//  coordinates_ = nodesCoordinates;
//  connectivity_ = elementsConnectivity;
//  dof_handle_dd_ = displ_dof_handle;
//  dof_handle_pressure_ = press_dof_handle;
//
//  fracture_id_=fractureID;
//  material_id_=materialID;
//  //condition_id_=conditionID;
//
//  interpolation_order_=interpolationOrder;
//
//};


// here constructor loading the

////////////////////////////////////////////////////////////////////////////////////////////

/// SETTER
void Mesh::appendMesh(const Mesh &newMesh, const bool isJoined) {

  //IL_EXPECT_FAST(interpolation_order_==newMesh.interpolation_order_);

  // Strategy for appendMesh method:
  // newMesh arrives and it can be either connected or not.
  // isJoined determines if the mesh will share some nodes with the previous one.
  //

  if(isJoined){
/* PLACEHOLDER
  // If isJoined == true, then the two meshes are connected VIA THE TIPS.
  // The is_tip vector determines which nodes in the mesh are tips. The following
  // step requires to check which nodes are in contacts, by computing the distance
  // of the two tip nodes and finding the pair with the minimum distance

  il::int_t connectedOldNode;
  il::int_t connectedNewNode;

  double distTol = 10 ^-6;

  for (il::int_t i = 0; i < number_nodes_; i++) {

    if (is_tip_[i]) {

      for (il::int_t j = 0; j < newMesh.number_nodes_; j++) {

        if (newMesh.is_tip_[j]) {

          double distance =
              sqrt(pow((coordinates_(i, 0) - newMesh.coordinates_(j, 0)), 2) + pow((coordinates_(i, 1) - newMesh.coordinates_(j, 1)), 2));

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
  // Nodes coordinates array is modified as following. The node to be collapsed is
  // stored in the variable connectedOldNode. That is maintained.
  // In the new mesh, the node that is collapsed is called connectedNewNode. That is
  // removed from the list and the other nodes are appended to the old list of
  // nodal coordinates. Then, in the connectivity matrix, the nodes with number less
  // than connectedNewNode are increased by oldNumNodes. connectedNewNode itself is
  // substituted by connectedOldNodes and the nodes with number greater than
  // connectedNewNodes are increased by oldNumNodes-1
  //
  //
  // 0,1,2,3....15,16 old nodes \
    //    common node is 13 and 7 ----> 0-16 (17 total old nodes), 0+17,1+17,...6+17,8+17-1,...11+17-1,12+17-1
  // 0,1,2,.....11,12 new nodes /
  //
  //
  // The node number 7 (connected new node) in the new mesh is omitted and substituted
  // in the connectivity by the number 13 (connected old node).
  // Handle matrices for displacements are appended increased by oldNumElements*2*(interpolation+1)
  // Handle matrices for pressure are appended with NumNodes offset before the connectedNewNode.
  // Connected new node value is substituted with the old connected value. The subsequent nodes
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

  coordinates_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
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




  // connectivity resize: here the 2 is because the element is 1D with 2 nodes
  // but it should rather be the order of the polynomia on the element plus 1
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
    // Connectivity matrices are concatenated. The new connectivity matrix contains
    // values from 0 to the new number of elements. Consequently, all the connectivity
    // will be shifted by the old number of elements.
    // Similarly, dof_handle matrices are shifted respectively by old number of nodes
    // (in the case of the pressure one) and num_elements*2*(interpolation_order+1)
    // (in the case of the displacement one).
    // Fracture ID and Material ID are saved as such.

    // Computing the new number of nodes (being 2 independent fractures, they are just summed)
    const il::int_t old_number_of_nodes = this->numNodes();
    const il::int_t new_number_of_nodes = this->numNodes() + newMesh.numNodes();

    // Resize the array of nodes coordinates accordingly and fill it with the new values
    coordinates_.resize(new_number_of_nodes,2);

    for(il::int_t i = 0; i < newMesh.numNodes(); i++){

      coordinates_(i + old_number_of_nodes,0) = newMesh.coordinates_(i,0);
      coordinates_(i + old_number_of_nodes,1) = newMesh.coordinates_(i,1);

    }

    // Compute the new number of elements (simply summed as well)
    const il::int_t old_number_of_elements = this->numElems();
    const il::int_t new_number_of_elements = this->numElems() + newMesh.numElems();

    // Resize connectivity accordingly and save the connectivity of new elements, which nodes identifier
    // have with a shift of old_number_of_elements w.r.t. the newMesh values

    // also materialID and fractureID are element related and concatenated here

    const il::int_t columns_conn_mtx = (interpolation_order_ == 0) ? 2 : (interpolation_order_+1);

    connectivity_.resize(new_number_of_elements, columns_conn_mtx);

    fracture_id_.resize(new_number_of_elements);

    material_id_.resize(new_number_of_elements);

    for(il::int_t i = 0; i < newMesh.numElems(); i++){

      for(il::int_t j = 0; j < columns_conn_mtx; j++ ) {

        connectivity_(i + old_number_of_elements, j) = newMesh.connectivity_(i, j) + old_number_of_nodes;

      }

      material_id_[i+old_number_of_elements]=newMesh.material_id_[i];
      fracture_id_[i+old_number_of_elements]=newMesh.fracture_id_[i];
    }

    //std::cout << "got connectivity" <<std::endl;

    // Displacement dof handles: concatenated and shifted by old_number_of_nodes * 2 * (interpOrd + 1)
    const il::int_t columns_displ_dofh_mtx = ((interpolation_order_ == 0) ? (2*2) : (2*(interpolation_order_+1)));
    const il::int_t displ_dof_handle_shift = old_number_of_elements * columns_displ_dofh_mtx;

    dof_handle_dd_.resize( new_number_of_elements, columns_displ_dofh_mtx );

    for (il::int_t i = 0; i < newMesh.numElems(); i++) {

      for(il::int_t j=0; j< columns_displ_dofh_mtx; j++) {

        dof_handle_dd_(i + old_number_of_elements, j) = newMesh.dof_handle_dd_(i, j) + displ_dof_handle_shift;

      }

    }

    //std::cout << "got displ dof" <<std::endl;


    // Pressure dof handles are concatenated and shifted by old_number_of_elements + 1
    const il::int_t columns_press_dofh_mtx = (interpolation_order_ == 0) ? 2 : (interpolation_order_+1);
    const il::int_t press_dof_handle_shift = old_number_of_elements * ( columns_press_dofh_mtx -1 ) + numFracs();

//    std::cout << "Here, concatenation of pressure" << std::endl;
//    std::cout << "Number of fractures is " << numFracs() << std::endl;

    dof_handle_pressure_.resize(new_number_of_elements,columns_press_dofh_mtx);

    for (il::int_t i=0; i< newMesh.numElems(); i++){

      for(il::int_t j=0; j< columns_press_dofh_mtx; j++){

        dof_handle_pressure_(i + old_number_of_elements, j) = newMesh.dof_handle_pressure_(i, j) + press_dof_handle_shift;

      }
    }

    //std::cout << "got press dof" <<std::endl;
  }
}


/*void Mesh::appendMesh(const il::Array2D<double> &newNodesCoordinates,
                const il::Array2D<il::int_t> &newElementsConnectivity,
                const il::Array<il::int_t> &newMaterialIdentifier) {

  // Initial assertions to check data consistency
  IL_EXPECT_FAST(newNodesCoordinates.size(0) > 0 && newNodesCoordinates.size(1) == 2);
  IL_EXPECT_FAST(newElementsConnectivity.size(0) > 0 && newElementsConnectivity.size(1) == 2);
  IL_EXPECT_FAST(newMaterialIdentifier.size() == newElementsConnectivity.size(0));

  const il::int_t oldNumNodes = number_nodes_;
  const il::int_t oldNumElements = number_elements_;

  number_nodes_ = number_nodes_ + newNodesCoordinates.size(0);
  number_elements_ = number_elements_ + newElementsConnectivity.size(0);

  coordinates_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
  for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
    coordinates_(i, 0) = newNodesCoordinates(i, 0);
    coordinates_(i, 1) = newNodesCoordinates(i, 1);
  }

  // connectivity resize: here the 2 is because the element is 1D with 2 nodes
  // but it should rather be the order of the polynomia on the element plus 1
  connectivity_.resize(number_elements_, 2);
  material_id_.resize(number_elements_);
  fracture_id_.resize(number_elements_);

  // add the values of the new mesh
  for (il::int_t i = oldNumElements; i < number_elements_; i++) {

    connectivity_(i, 0) = newElementsConnectivity(i, 0);
    connectivity_(i, 1) = newElementsConnectivity(i, 1);

    material_id_[i] = newMaterialIdentifier[i];
    fracture_id_[i] = 0;

  }

};


void Mesh::appendMesh(const il::Array2D<double> &newNodesCoordinates,
                const il::Array2D<il::int_t> &newElementsConnectivity,
                const il::Array<il::int_t> &newMaterialIdentifier,
                const il::Array<il::int_t> &newFractureIdentifier) {

  // Initial assertions to check data consistency
  IL_EXPECT_FAST(newNodesCoordinates.size(0) > 0 && newNodesCoordinates.size(1) == 2);
  IL_EXPECT_FAST(newElementsConnectivity.size(0) > 0 && newElementsConnectivity.size(1) == 2);
  IL_EXPECT_FAST(newMaterialIdentifier.size() == newElementsConnectivity.size(0));
  IL_EXPECT_FAST(newFractureIdentifier.size() == newElementsConnectivity.size(0));

  const il::int_t oldNumNodes = number_nodes_;
  const il::int_t oldNumElements = number_elements_;

  number_nodes_ = number_nodes_ + newNodesCoordinates.size(0);
  number_elements_ = number_elements_ + newElementsConnectivity.size(0);

  coordinates_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
  for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
    coordinates_(i, 0) = newNodesCoordinates(i, 0);
    coordinates_(i, 1) = newNodesCoordinates(i, 1);
  }

  // connectivity resize: here the 2 is because the element is 1D with 2 nodes
  // but it should rather be the order of the polynomia on the element plus 1
  connectivity_.resize(number_elements_, 2);
  material_id_.resize(number_elements_);
  fracture_id_.resize(number_elements_);

  // add the values of the new mesh
  for (il::int_t i = oldNumElements; i < number_elements_; i++) {

    connectivity_(i, 0) = newElementsConnectivity(i, 0);
    connectivity_(i, 1) = newElementsConnectivity(i, 1);

    material_id_[i] = newMaterialIdentifier[i];
    fracture_id_[i] = newFractureIdentifier[i];

  }

};

void Mesh::appendNodeToMeshTip(const il::int_t mesh_node, const double x_new, const double y_new) {
  //// JUST FOR LINEAR ELEMENT FOR THE MOMENT WITH UNIFORM CONDITIONS

  // check if mesh_node is tip
  if (is_tip_[mesh_node]) {

    il::int_t oldTipElement;
    // if the old mesh node is tip, there is only one element which has the old mesh node
    // Consequently it will appear once in the connectivity matrix.
    il::int_t i = 0;
    while (i < number_elements_ &&
        (connectivity_(i, 0) != mesh_node || connectivity_(i, 1) != mesh_node)) {
      i++; // goes to the next row in the matrix if mesh_node is not found
    }
    // once the mesh node is found, we exit from the while loop
    // then, i is the element of the old mesh node
    oldTipElement = i;

    // increase number of nodes by 1 and add x_new, y_new
    number_nodes_ = number_nodes_ + 1;

    coordinates_.resize(number_nodes_, 2);
    coordinates_(number_nodes_, 0) = x_new;
    coordinates_(number_nodes_, 1) = y_new;

    // increase number of elements by 1 and create the connectivity (mesh_node, newNode)
    number_elements_ = number_elements_ + 1;

    connectivity_.resize(number_elements_, 2);
    connectivity_(number_elements_, 0) = mesh_node;
    connectivity_(number_elements_, 1) = number_nodes_;

    // dof_handle_displacement: for the new element (a new line at the end of the dof table),
    // take the total size of the matrix and add 4 dofs
    dof_handle_dd_.resize(number_elements_, 4);
    dof_handle_dd_(number_elements_, 0) = number_elements_ * 4;
    dof_handle_dd_(number_elements_, 1) = number_elements_ * 4 + 1;
    dof_handle_dd_(number_elements_, 2) = number_elements_ * 4 + 2;
    dof_handle_dd_(number_elements_, 3) = number_elements_ * 4 + 3;

    // dof_handle_pressure: for the new element (nodal!) add a new line with the old mesh node
    // and new mesh node
    dof_handle_pressure_.resize(number_nodes_, 2);

    dof_handle_pressure_(number_nodes_, 0) = mesh_node;
    dof_handle_pressure_(number_nodes_, 1) = number_nodes_;

    // for the new element: same fractureID, same materialID;
    fracture_id_.append(fracture_id_[oldTipElement]);
    material_id_.append(material_id_[oldTipElement]);

    // for the new collocation points: far_field_stress_condition_id_ same as before [to be discussed how to change that]
    // far_field_stress_condition_id_ is defined as
    for (il::int_t i = 0; i < 4; i++) {
      far_field_stress_condition_id_.append(far_field_stress_condition_id_[oldTipElement]);  /// ONLY FOR UNIFORM STRESS
    }

    // for the new nodes: pressure_condition_id_  same as before [to be discussed how to change that]
    pressure_condition_id_.append(pressure_condition_id_[mesh_node]);

    // for the source id vector: add a line with -1 (no source considered)
    source_id_.append(-1);

    // assign the new tip to the new node
    is_tip_[mesh_node] = false;
    is_tip_[number_nodes_] = true;

  } else {
    std::cerr << "Selected node is not tip" << std::endl;
    exit(4);
  }

};*/



// mesh class
//void Mesh::loadMesh(il::Array2D<double> xy, il::Array2D<int> ien, il::Array<int> mat) {
//  // check array dimensions ?? -> this is only for 1D mesh so far
//  IL_EXPECT_FAST(xy.size(1) == 2);  // check array dimensions ?
//  IL_EXPECT_FAST(ien.size(1) == 2);
//
//  IL_EXPECT_FAST(ien.size(0) == mat.size());
//
//  coordinates_ = xy;           // list of coordinates of points in the mesh
//  connectivity_ = ien;  //  connectivity array -
//  material_id_ = mat; // material ID array
//
//  // one could think of having a FracID ...
//}

//double Mesh::node(il::int_t k, il::int_t i) const { return coordinates_(k, i); }
//
//int Mesh::connectivity(il::int_t k, il::int_t i) const { return connectivity_(k, i); }
//
//int Mesh::matid(il::int_t k) const { return material_id_[k]; }
//
//int Mesh::nelts() const { return connectivity_.size(0); } ;
//
//int Mesh::ncoor() const { return coordinates_.size(0); };
//
//il::Array2D<double> Mesh::coor() const { return coordinates_; };
//
//il::Array2D<int> Mesh::conn() const { return connectivity_; };
//
//il::Array<int> Mesh::matid() const { return material_id_; };


// needs to add function to add one or more elements ... (needs to have active
// and passive elements to track active/passive fractures etc.)
//
//
// could provide a default constructor for a straight fracture ?



/*// SOME UTILITIES HERE below -> To be moved in a separate file ??
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
}*/
//----------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
//          METHODS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// get size of an element
//double Mesh::elt_size(const il::int_t e) {
//  il::StaticArray<double, 2> xdiff;
//  xdiff[0] = coordinates_(connectivity_(e, 1), 0) -
//      coordinates_(connectivity_(e, 0), 0);
//  xdiff[1] = coordinates_(connectivity_(e, 1), 1) -
//      coordinates_(connectivity_(e, 0), 1);
//  double hx = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));
//
//  return hx;
//};

////////////////////////////////////////////////////////////////////////////////
// get size of a all elements - output is an array of double
//il::Array<double> Mesh::All_elt_size() {
//  il::Array<double> temp(connectivity_.size(0));
//  for (il::int_t i = 0; i < numElems(); i++) {
//    temp[i] = elt_size(i);
//  };
//  return temp;
//}

////////////////////////////////////////////////////////////////////////////////
//il::Array2D<il::int_t> Mesh::GetNodesSharing2Elts() {
//  // case of only 2 nodes for now.....
//  // needed for building FD matrix  without fracture intersection.....
//  // ONLY WORK ON MESH where nodes have maximum of 2 adjacaent elements
//  // ..... no fracture intersection ...
//  // todo generalize to account for fracture intersection
//
//  IL_EXPECT_FAST(connectivity_.size(1) == 2);
//
//  il::Array<il::int_t> tip = tipnodes_;
//
//  il::Array2D<il::int_t> temp(numNodes() - tip.size(), 2, 0);
//
//  il::int_t k = 0;
//  for (il::int_t i = 0; i < coordinates_.size(0); i++) {
//    if (node_elt_connectivity(i, 1) > -1) {
//      temp(k, 0) = node_elt_connectivity(i, 0);
//      temp(k, 1) = node_elt_connectivity(i, 1);
//      k++;
//    }
//  }
//  return temp;
//}

////////////////////////////////////////////////////////////////////////////////
// method to get next to tip elements - so-called ribbons element in an
// ILSA-like scheme

//il::Array<il::int_t> Mesh::getRibbonElements(){
//
//
//  il::Array<il::int_t> ribbon_elts{tipelts_.size(),0};
//  il::int_t  nt_1;
//
//  for (il::int_t e=0; e< ribbon_elts.size(); e++){
//
//    // local connectivity of the tip elements
//    if( connectivity(tipelts_[e],0) == tipnodes_[e] ){
//      nt_1 = connectivity(tipelts_[e],1);
//    } else {
//      if (connectivity(tipelts_[e],1) == tipnodes_[e] ) {
//        nt_1 = connectivity(tipelts_[e],0);
//      } else
//      {
//        il::abort();
//      }
//    }
//
//    if (node_elt_connectivity(nt_1,0)== tipelts_[e]){
//      ribbon_elts[e]=node_elt_connectivity(nt_1,1);
//    }
//    else
//    {
//      if (node_elt_connectivity(nt_1,1)== tipelts_[e]) {
//        ribbon_elts[e]=node_elt_connectivity(nt_1,0);
//      } else
//      {il::abort();
//      }
//    }
//
//
//  }
//
//  return ribbon_elts;
//
//}

////////////////////////////////////////////////////////////////////////////////
// could do a method returning the element and nodes of a given fracture
/// this is honestly not really needed for now.....

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
// adding elements function...
//void Mesh::AddNTipElements(const il::int_t t_e, const il::int_t the_tip_node,
//                           const il::int_t n_add, double kink_angle) {
//  // add n_add elements in the Mesh object ahead of the nodes the_tip_node
//  // (global numbering) of
//  // element t_e
//  // with a kink_angle   with respect to element t_e
//  //  The size of the added elements are equal to the size of element of t_e
//  // the kick_angle should be given in the local tip coordinates system.....
//  //
//  // api could be changed ...
//
//  // we need the segment data .....
//  hfp2d::SegmentData tipEltData = Mesh::getElementData(t_e);
//
//  il::StaticArray<il::int_t, 2> tipEltConn = connectivity(t_e);
//  double h = tipEltData.size(), global_prop_angle;
//
//  il::int_t local_tip_node;
//  il::StaticArray<double, 2> local_dir;
//
//  // we need to know how if the tip nodes is the first or second nodes of the
//  // element to know the direction of propagation
//
//  if (the_tip_node == tipEltConn[0]) {
//    local_tip_node = 0;
//    global_prop_angle = tipEltData.theta() + il::pi - kink_angle;
//  } else {
//    if (the_tip_node == tipEltConn[1]) {
//      local_tip_node = 1;
//      global_prop_angle = tipEltData.theta() + kink_angle;
//    } else {
//      std::cout << "ERROR in AddNtipElements \n ";
//      il::abort();  // must throw an exception here
//    }
//  };
//
//  il::StaticArray<double, 2> prop_dir;  // in global coordinates....
//  prop_dir[0] = h * cos(global_prop_angle);
//  prop_dir[1] = h * sin(global_prop_angle);
//
//  il::Array2D<double> new_nodes(n_add, 2);
//  il::StaticArray<double, 2> tip_nodes_coor = coordinates(the_tip_node);
//
//  // coordinates of the new nodes to be added
//  for (il::int_t i = 0; i < n_add; i++) {
//    if (i == 0) {
//      new_nodes(i, 0) = tip_nodes_coor[0] + prop_dir[0];
//      new_nodes(i, 1) = tip_nodes_coor[1] + prop_dir[1];
//    } else {
//      new_nodes(i, 0) = new_nodes(i - 1, 0) + prop_dir[0];
//      new_nodes(i, 1) = new_nodes(i - 1, 1) + prop_dir[1];
//    }
//  }
//
//  // reconstructing the whole coordinates array (sub-optimal)
//  il::Array2D<double> new_all_coor(numNodes() + n_add, 2);
//
//  for (il::int_t i = 0; i < numNodes(); i++) {
//    new_all_coor(i, 0) = coordinates(i, 0);
//    new_all_coor(i, 1) = coordinates(i, 1);
//  }
//  //
//  for (il::int_t i = 0; i < n_add; i++) {
//    new_all_coor(i + numNodes(), 0) = new_nodes(i, 0);
//    new_all_coor(i + numNodes(), 1) = new_nodes(i, 1);
//  }
//
//  // duplicating the old connectivity table....
//  il::Array2D<il::int_t> new_conn(numElems() + n_add, 2);
//  for (il::int_t i = 0; i < numElems(); i++) {
//    new_conn(i, 0) = connectivity(i, 0);
//    new_conn(i, 1) = connectivity(i, 1);
//  }
//  // adding the new connectivity at the end (so old element numbers are still
//  // the same).
//  for (il::int_t i = 0; i < n_add; i++) {
//    if (i == 0) {
//      new_conn(i + numElems(), 0) = connectivity(t_e, local_tip_node);
//      new_conn(i + numElems(), 1) = numNodes() + i;
//    } else {
//      new_conn(i + numElems(), 0) = numNodes() + i - 1;
//      new_conn(i + numElems(), 1) = numNodes() + i;
//    }
//  };
//
//  // now UPDATE THE MESH....
//  coordinates_ = new_all_coor;
//  connectivity_ = new_conn;
//
//  //  the other changes....
//  // this is for uniform material only
//  il::Array<il::int_t> material_id_(connectivity_.size(0), 1);
//
//  il::int_t nelts = connectivity_.size(0);
//  il::int_t p = interpolation_order_;
//
//  // COULD to be optimized below.... here we re-built everything from scratch...
//  // anyway copy would be needed....
//
//  /// Discontinuous Polynomial DOF handles
//  il::Array2D<il::int_t> id_dd{nelts, 2 * (p + 1), 0};
//  for (il::int_t i = 0; i < nelts; i++) {
//    for (il::int_t j = 0; j < 2 * (p + 1); j++) {
//      id_dd(i, j) = i * 2 * (p + 1) + j;
//    }
//  }
//  dof_handle_dd_ = id_dd;  /// dof
//
//  /// //    dof(element, local nnodes number)
//  // actually this is the connectivity_ array for  p =1 and
//  // a simple elt number of P0
//  switch (interpolation_order_) {
//    case 0: {
//      il::Array2D<il::int_t> id_press{nelts, 1, 0};
//      for (il::int_t e = 0; e < nelts; e++) {
//        id_press(e, 0) = e;
//      };
//      dof_handle_pressure_ = id_press;
//    }
//    case 1:
//      dof_handle_pressure_ = connectivity_;  // 1 unknowns per nodes ....
//  };
//
////  // rebuild the nodal connected table...
////  node_adj_elt_ = GetNodalEltConnectivity(coordinates_.size(0), connectivity_);
////
////  // rebuilt tip nodes table...
////  tipnodes_ = BuildTipNodes(node_adj_elt_);
////  tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
//}

}