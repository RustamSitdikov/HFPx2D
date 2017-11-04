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
#include <cmath>
#include <iostream>
#include <algorithm>

// Inclusion from Inside Loop library
#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>
#include <il/container/1d/SmallArray.h>
#include <il/String.h>
#include <il/linear_algebra.h>

// Inclusion from hfp2d
#include <src/core/SegmentData.h>

namespace hfp2d
{

il::Array2D<il::int_t> GetNodalEltConnectivity(
    const il::int_t nt_nodes, const il::Array2D<il::int_t> &connectivity);

il::Array<il::int_t> BuildTipNodes(
    const il::Array2D<il::int_t> &node_connectivity);

il::Array<il::int_t> BuildTipElts(
    const il::Array2D<il::int_t> &node_connectivity,
    const il::Array<il::int_t> &tipnodes);

///// 1D mesh class
class Mesh
{  // class for 1D mesh of 1D segment elements ?

private:

    // Coordinates of the nodes - size: number of nodes x problem dimension (2D)
    il::Array2D<double> coordinates_;

    // Connectivity matrix - size: number of elements x (order interpolation + 1)
    il::Array2D<il::int_t> connectivity_;

    // Interpolation order
    il::int_t interpolation_order_;

    // Dof handle matrices
    // for displacements - size: number of elements x 2dofs per node x (order interpolation + 1)
    il::Array2D<il::int_t> dof_handle_dd_;
    // for pressure - size: number of nodes x 1dof per node x (order interpolation + 1)
    il::Array2D<il::int_t> dof_handle_pressure_;

    // Identifier number of the fracture - size: number of elements
    il::Array<il::int_t> fracture_id_;

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
    Mesh(){};  // TODO: remove empty initialization of mesh class variables if
    // possible.

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
        node_adj_elt_ =
            GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

        // built tip nodes table...
        tipnodes_ = BuildTipNodes(node_adj_elt_);
        tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
    };

    // case where matid vector is provided
    // constructor with interpolation order and coordinates and connectivity array
    Mesh(const il::Array2D<double> &Coordinates,
         const il::Array2D<il::int_t> &Connectivity,
         const il::Array<il::int_t> &MatID, const il::int_t interpolationOrder) {
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
        node_adj_elt_ =
            GetNodalEltConnectivity(coordinates_.size(0), connectivity_);

        // built tip nodes table...
        tipnodes_ = BuildTipNodes(node_adj_elt_);
        tipelts_ = BuildTipElts(node_adj_elt_, tipnodes_);
    };



//    // Constructor with only nodes and elements.
//    Mesh(const il::Array2D<double> &nodesCoordinates,
//         const il::Array2D<il::int_t> &elementsConnectivity)
//    {
//
//        coordinates_ = nodesCoordinates;
//        connectivity_ = elementsConnectivity;
//
//    };
//
//    Mesh(const il::int_t interpolationOrder,
//         const il::Array2D<double> &nodesCoordinates,
//         const il::Array2D<il::int_t> &elementsConnectivity,
//         const il::Array2D<il::int_t> &displacementDOFHandle)
//    {
//
//        interpolation_order_ = interpolationOrder;
//        coordinates_ = nodesCoordinates;
//        connectivity_ = elementsConnectivity;
//        dof_handle_dd_ = displacementDOFHandle;
//
//    };
//
    Mesh(const il::int_t interpolationOrder,
         const il::Array2D<double> &nodesCoordinates,
         const il::Array2D<il::int_t> &elementsConnectivity,
         const il::Array2D<il::int_t> &displ_dof_handle,
         const il::Array2D<il::int_t> &press_dof_handle,
         const il::Array<il::int_t> &fractureID,
         const il::Array<il::int_t> &materialID,
         const il::Array<il::int_t> &conditionID);

    ////////////////////////////////////////////////////////////////////////////////////////////

/// SETTER
    void
    appendMesh(const Mesh &newMesh,
               bool isJoined);

    void
    appendMesh(const il::Array2D<double> &newNodesCoordinates,
               const il::Array2D<il::int_t> &newElementsConnectivity,
               const il::Array<il::int_t> &newMaterialIdentifier);

    void
    appendMesh(const il::Array2D<double> &newNodesCoordinates,
               const il::Array2D<il::int_t> &newElementsConnectivity,
               const il::Array<il::int_t> &newMaterialIdentifier,
               const il::Array<il::int_t> &newFractureIdentifier);

    void
    appendNodeToMeshTip(il::int_t mesh_node,
                        double x_new,
                        double y_new);

    //////////////////////////////////////////////////////////////////////////
    //        get functions  - i.e. public interfaces
    //////////////////////////////////////////////////////////////////////////

    // number of nodes
    il::int_t numNodes() const { return coordinates_.size(0); }
    // number of elements
    il::int_t numElems() const { return connectivity_.size(0); }
    // interpolation order
    il::int_t interpOrd() const { return interpolation_order_; }

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


    il::int_t
    numFracs() const
    {

        auto thePosition =
            std::max_element(fracture_id_.begin(), fracture_id_.end());
        il::int_t theValue = *thePosition;

//    std::cout << "Position " << thePosition << std::endl;
//    std::cout << "Begin " << fracture_id_.begin() << " End " << fracture_id_.end() << std::endl;
//    std::cout << "Value " << theValue << " +1 " << theValue+1 << std::endl;

        return (theValue + 1);

    }

    il::int_t
    numMats() const
    {
        return (*std::max_element(material_id_.begin(), material_id_.end()))
            + 1;
    }

    il::int_t
    numDDDofsPerElem() const { return dof_handle_dd_.size(1); }
    il::int_t
    numPressDofsPerElem() const { return dof_handle_pressure_.size(1); }

    il::int_t
    numPressDofs() const
    {
        return (numElems() * interpolation_order_ + numFracs());
    }

    il::int_t
    numDDDofs() const
    {
        return (numElems() * (interpolation_order_ + 1) * 2);
    }


    // Read a particular element of the node coordinates
    double
    node(il::int_t k, il::int_t i) const { return coordinates_(k, i); }
    // todo: restructure "node" method in order to have mesh.node(i).X for the x coordinate of node i ??

    il::Array<il::int_t>
    elemConnectivity(il::int_t k)
    {
        il::Array<il::int_t> temp(connectivity_.size(1));

        for (il::int_t i = 0; i < connectivity_.size(1); i++)
        {
            temp[i] = connectivity_(k, i);
        }

        return temp;
    };

    il::int_t
    connectivity(il::int_t k, il::int_t i) const { return connectivity_(k, i); }
    il::int_t
    dofPress(il::int_t k, il::int_t i) const
    {
        return dof_handle_pressure_(k,
                                    i);
    }
    il::int_t
    dofDD(il::int_t k,
          il::int_t i) const { return dof_handle_dd_(k, i); }
    il::int_t
    fracID(il::int_t k) const { return fracture_id_[k]; }
    il::int_t
    matID(il::int_t k) const { return material_id_[k]; }


    //il::int_t matid(il::int_t k) const { return material_id_[k]; }

    il::int_t
    nelts() const { return connectivity_.size(0); };

    il::int_t
    ncoor() const { return coordinates_.size(0); };

    il::Array2D<double>
    coor() const { return coordinates_; };

    il::Array2D<il::int_t>
    conn() const { return connectivity_; };

    il::Array<il::int_t>
    matid() const { return material_id_; };

    double
    eltsize(il::int_t e)
    {
        // NOTE: size of element only in the case of constant or linear elements

        il::StaticArray<double, 2> xdiff;
        xdiff[0] =
            coordinates_(connectivity_(e, 1), 0) - coordinates_(connectivity_(e, 0), 0);
        xdiff[1] =
            coordinates_(connectivity_(e, 1), 1) - coordinates_(connectivity_(e, 0), 1);
        double hx = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

        return hx;

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
//    il::int_t connectivity(il::int_t e, il::int_t i) const {
//        // element e, local coordinates i - return global nodes
//        return connectivity_(e, i);
//    }

    // nodal connectivity related
    il::Array2D<il::int_t> node_elt_connectivity() const {
        return node_adj_elt_;
    };
    il::int_t node_elt_connectivity(il::int_t k, il::int_t l) const {
        return node_adj_elt_(k, l);
    };

    il::Array<il::int_t> node_elt_connectivity(il::int_t k) const {
        il::Array<il::int_t> temp(node_adj_elt_.size(1));
        for (il::int_t i = 0; i < node_adj_elt_.size(1); i++) {
            temp[i] = node_adj_elt_(k, i);
        }
        return temp;
    };

    // get tip nodes
    il::Array<il::int_t> tip_nodes() const { return tipnodes_; };
    il::int_t tip_nodes(il::int_t k) const { return tipnodes_[k]; };

    // get tip elts
    il::Array<il::int_t> tip_elts() const { return tipelts_; };
    il::int_t tip_elts(il::int_t k) const { return tipelts_[k]; };

    // material ID related
//    il::Array<il::int_t> matid() const { return material_id_; };
    il::int_t matid(il::int_t k) const { return material_id_[k]; }

    il::int_t numberOfMaterials() const {
        return (*std::max_element(material_id_.begin(), material_id_.end()) + 1);
    }


    // interpolation order
    il::int_t interpolationOrder() const { return interpolation_order_; }

    // dofs related.....
    il::int_t DDDofsPerElement() const { return dof_handle_dd_.size(1); }

    il::int_t PressDofsPerElement() const { return dof_handle_pressure_.size(1); }

    il::int_t numberOfPressDofs() const {
        return dof_handle_pressure_.size(0);
    };  // this is a scalar continuous polynomial.

    il::int_t numberOfDDDofs() const {
        return (numElems() * (interpolation_order_ + 1) * 2);
    }

//    il::int_t dofPress(il::int_t k, il::int_t i) const {
//        // coordinates k, dof i -> return global equation iD
//        return dof_handle_pressure_(k, i);
//    }

//    il::int_t dofDD(il::int_t k, il::int_t i) const {
//        // coordinates k, dof i -> return global equation iD
//        return dof_handle_dd_(k, i); // element , dof dim.
//    }


    ////////////////////////////////////////////////////////////////////////////////////////////
    //   Methods
    ////////////////////////////////////////////////////////////////////////////////////////////

    hfp2d::SegmentData getElementData(const il::int_t ne);

    // a method to get the size of a given element.
    double elt_size(il::int_t &e);

    il::Array<double> All_elt_size();

    // method to get nodes sharing 2 elts. (i.e. the nodal_connectivity for all nodes with 2 neighbours)
    // todo rename
    il::Array2D<il::int_t> GetNodesSharing2Elts();

    // method to get the ribbon elements  - of a given mesh.
    il::Array<il::int_t>  getRibbonElements() ;


    // method to add N element ahead of a tip node of a tip element at a given kick angle
    void AddNTipElements(const il::int_t t_e, const il::int_t the_tip_node,
                         const il::int_t n_add, double kink_angle);

    // todo : do comments each methods better !


    // TODO: remove all methods that are not needed

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
