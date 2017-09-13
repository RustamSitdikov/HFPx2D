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
#include <iostream>

namespace hfp2d {

///// 1D mesh class
class Mesh {  // class for 1D mesh of 1D segment elements ?

private:

  // Number of nodes (default initialized at zero)
  il::int_t number_nodes_ = 0;
  // Number of elements (default initialized at zero)
  il::int_t number_elements_ = 0;
  // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
  il::Array2D<double> nodes_;
  // Connectivity matrix - size: number of elements x (order interpolation + 1)
  il::Array2D<il::int_t> connectivity_;

  // Dof handle matrices
  // for displacements - size: number of elements x 2dofs per node x (order interpolation + 1)
  il::Array2D<il::int_t> dof_handle_displacement_;
  // for pressure - size: number of elements x 1dof per node x (order interpolation + 1)
  il::Array2D<il::int_t> dof_handle_pressure_;

  // Identifier number of the fracture - size: number of elements
  il::Array<il::int_t> fracture_id_;
  // Material identifier - size: number of elements
  il::Array<il::int_t> material_id_;

  // Condition identifiers for displacement and pressure, since stress is per collocation point
  // and pressure is for nodes
  il::Array<il::int_t> far_field_stress_condition_id_; // size: number of collocation nodes
  il::Array<il::int_t> pressure_condition_id_; // size: number of nodes

  // Source identifier, per node. It substitutes is_source.
  il::Array<il::int_t> source_id_;

  // Tip identifier
  il::Array<bool> is_tip_;

public:

  //////////////////////////////// CONSTRUCTORS (a.k.a. initializers) ////////////////////////////////

  Mesh(      il::int_t              interpolationOrder,
       const il::Array2D<double>    &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array<il::int_t>   &sourceIdentifier) ;


  Mesh(      il::int_t              interpolationOrder,
       const il::Array2D<double>    &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t>   &sourceIdentifier);

  // Initialization without tip identifier, for static mesh case
  Mesh(      il::int_t              interpolationOrder,
       const il::Array2D<double>    &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t>   &fractureIdentifier,
       const il::Array<il::int_t>   &materialIdentifier,
       const il::Array<il::int_t>   &farStressCondId,
       const il::Array<il::int_t>   &porePressCondId,
       const il::Array<il::int_t>   &sourceIdentifier);

  // First constructor: it initializes a mesh by using all available information on the mesh
  Mesh(      il::int_t              interpolationOrder,
       const il::Array2D<double>    &nodesCoordinates,
       const il::Array2D<il::int_t> &elementsConnectivity,
       const il::Array2D<il::int_t> &dofHandleDisplacement,
       const il::Array2D<il::int_t> &dofHandlePressure,
       const il::Array<il::int_t>   &fractureIdentifier,
       const il::Array<il::int_t>   &materialIdentifier,
       const il::Array<il::int_t>   &farStressCondId,
       const il::Array<il::int_t>   &porePressCondId,
       const il::Array<il::int_t>   &sourceIdentifier,
       const il::Array<bool>        &tipIdentifier);


  ////////////////////////////////////////////////////////////////////////////////////////////




  /// SETTER
  void appendMesh(hfp2d::Mesh &newMesh, bool isJoined) {

    // Strategy for appendMesh method:
    // newMesh arrives and it can be either connected or not.
    // isJoined determines if the mesh will share some nodes with the previous one.
    //
    // If isJoined == false, then the two meshes are not connected.
    // Hence, the total number of nodes and elements are summed together.
    // Connectivity matrices are concatenated. The new connectivity matrix contains
    // values from 0 to the new number of elements. Consequently, all the connectivity
    // will be shifted by the old number of elements.
    // Similarly, dof_handle matrices are shifted respectively by num_nodes (in the case
    // of the pressure one) and num_elements*2*(interpolation_error+1) (in the case
    // of the displacement one).
    // Fracture ID is the old one +1. Material ID, condition ID, far field stress ID,
    // pressure condition ID, source ID are maintained as is and concatenated.
    // The vector is_tip is concatenated as well.
    //
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
                sqrt(pow((nodes_(i, 0) - newMesh.nodes_(j, 0)), 2) + pow((nodes_(i, 1) - newMesh.nodes_(j, 1)), 2));

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

    nodes_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
    for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
      nodes_(i, 0) = newMesh.nodes_(i, 0);
      nodes_(i, 1) = newMesh.nodes_(i, 1);
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

    }

  };

  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
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

    nodes_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
    for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
      nodes_(i, 0) = newNodesCoordinates(i, 0);
      nodes_(i, 1) = newNodesCoordinates(i, 1);
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

  void appendMesh(const il::Array2D<double> &newNodesCoordinates,
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

    nodes_.resize(number_nodes_, 2); // here the 2 is because we are working in 2D
    for (il::int_t i = oldNumNodes; i < number_nodes_; i++) {
      nodes_(i, 0) = newNodesCoordinates(i, 0);
      nodes_(i, 1) = newNodesCoordinates(i, 1);
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

  void appendNodeToMeshTip(il::int_t mesh_node, double &x_new, double &y_new) {
    //// JUST FOR LINEAR ELEMENT FOR THE MOMENT

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

      nodes_.resize(number_nodes_, 2);
      nodes_(number_nodes_, 0) = x_new;
      nodes_(number_nodes_, 1) = y_new;

      // increase number of elements by 1 and create the connectivity (mesh_node, newNode)
      number_elements_ = number_elements_ + 1;

      connectivity_.resize(number_elements_, 2);
      connectivity_(number_elements_, 0) = mesh_node;
      connectivity_(number_elements_, 1) = number_nodes_;

      // dof_handle_displacement: for the new element (a new line at the end of the dof table),
      // take the total size of the matrix and add 4 dofs
      dof_handle_displacement_.resize(number_elements_, 4);
      dof_handle_displacement_(number_elements_, 0) = number_elements_ * 4;
      dof_handle_displacement_(number_elements_, 1) = number_elements_ * 4 + 1;
      dof_handle_displacement_(number_elements_, 2) = number_elements_ * 4 + 2;
      dof_handle_displacement_(number_elements_, 3) = number_elements_ * 4 + 3;

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
        far_field_stress_condition_id_.append(far_field_stress_condition_id_[oldTipElement]);
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

  }

  /// GETTER
  // Read sizes of matrices
  il::int_t numberOfNodes() const { return number_nodes_; };
  il::int_t numberOfElements() const { return number_elements_; };

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
