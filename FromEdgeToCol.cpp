//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "FromEdgeToCol.h"
#include "DOF_Handles.h"
#include <cmath>
#include <il/linear_algebra.h>
#include <iostream>


// Function that allow us to switch from end points (two values per node -> dof_dim = 2) to collocation points
// Remember: the elasticity is evaluated at collocation points
// The collocation points are located in the reference element at location {-1/sqrt(2) , +1/sqrt(2)}
// It is general ina sense that it includes both opening & shear

il::Array<double> FromEdgeToCol(il::Array<double> &d_edge, il::Array2D<int> id, int Nelts, int dof_dim){

    il::Array<double> Fetc{d_edge.size(),0.};
    il::Array2D<double> ShapeFunction{4,4,.0};

    ShapeFunction(0,0) = (1+(1/sqrt(2)))/2;
    ShapeFunction(0,1) = (1+(1/sqrt(2)))/2;
    ShapeFunction(0,2) = (1-(1/sqrt(2)))/2;
    ShapeFunction(0,3) = (1-(1/sqrt(2)))/2;

    ShapeFunction(1,0) = ShapeFunction(0,0);
    ShapeFunction(1,1) = ShapeFunction(0,1);
    ShapeFunction(1,2) = ShapeFunction(0,2);
    ShapeFunction(1,3) = ShapeFunction(0,3);

    ShapeFunction(2,0) = (1-(1/sqrt(2)))/2;
    ShapeFunction(2,1) = (1-(1/sqrt(2)))/2;
    ShapeFunction(2,2) = (1+(1/sqrt(2)))/2;
    ShapeFunction(2,3) = (1+(1/sqrt(2)))/2;

    ShapeFunction(3,0) = ShapeFunction(2,0);
    ShapeFunction(3,1) = ShapeFunction(2,1);
    ShapeFunction(3,2) = ShapeFunction(2,2);
    ShapeFunction(3,3) = ShapeFunction(2,3);


    for (il::int_t i = 0, k =0; i < Nelts; ++i) {

        for (int j = 0; j < 2*dof_dim; j = j+2) {


            Fetc[id(i,0)] = Fetc[id(i,0)] + ShapeFunction(0, j) * d_edge[id(i, j)];
            Fetc[id(i,1)] = Fetc[id(i,1)] + ShapeFunction(1, j) * d_edge[id(i, j)];
            Fetc[id(i,2)] = Fetc[id(i,2)] + ShapeFunction(2, j) * d_edge[id(i, j)];
            Fetc[id(i,3)] = Fetc[id(i,3)] + ShapeFunction(3, j) * d_edge[id(i, j)];

        }

    }

    return Fetc;
};