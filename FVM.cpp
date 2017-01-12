//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "FVM.h"
#include "Mesh.h"
#include "Dilatancy.h"
#include <iostream>

/////// Some utilities ///////

// This function calculates the average between two values
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector that contain the average values of each row (element) of the input matrix

il::Array<double> Average(il::Array2D<double> &d){

    il::Array<double> Average{d.size(0),0.};

    for (il::int_t i = 0; i < d.size(0); ++i) {

            Average[i] = (d(i,0) + d(i,1))/2;

    }

    return Average;

};


// This function calculates the slip/opening at +/- 1/4 -> the control value is centered on the nodes!
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector -> {slip_+1/4 , slip_+3/4}

il::Array<double> Quarter(il::Array2D<double> &d){

  il::Array<double> Quarter(2*d.size(0),0.);

    for (il::int_t i = 0, j=0; i < (d.size(0)); ++i, j=j+2) {

        Quarter[j] = ((3*d(i,0))+d(i,1))/4;
        Quarter[j+1] = (d(i,0)+(3*d(i,1)))/4;

    }

    return Quarter;

};

// Function to find out the position of a value in a 2D array
// It returns 2x2 array with row&col of the seek value
il::Array2D<int> FindPosit_2DArray(il::Array2D<int> &arr2D, double_t seek) {

    il::Array2D<int> outp{2,2,0.};

    for (il::int_t i = 0; i < arr2D.size(0); ++i) {

        for (il::int_t j = 0; j < arr2D.size(1); ++j) {

            if (arr2D(i, j) == seek){

                outp(j,0) = i;
                outp(j,1) = j;

            }
        }
    }

    return outp;
}


// Auxiliary function for assembly process
// It returns the a given row (vector - specified by idx) of a 2D array
il::Array<int> Auxiliary(il::Array2D<int> &arr, il::int_t idx){


    il::Array<int> vect{arr.size(1),0.};

    for (il::int_t i = 0; i < vect.size(); ++i) {

        vect[i] = arr(idx,i);

    }

    return vect;

};


//Function for the coefficients of the Finite Difference Matrix "L"
// Output: array (vector) that contains all the coefficients for each element

il::Array<double> ConductivitiesNewtonian(const int Visc, Mesh mesh, il::Array2D<double> rho, il::Array2D<double> &d, const double Incr_dil, const double d_wd, const double Init_dil){


    il::Array<double> d_mid, wh_mid, rho_mid;

    d_mid = Average(d);
    wh_mid = Dilatancy(Init_dil,Incr_dil,d_wd,d_mid);
    rho_mid = Average(rho);

    //create the array of element size
    il::Array<double> EltSizes{(mesh.conn).size(0), 0.};

    for (il::int_t i = 0; i < EltSizes.size(); ++i) {

        for (il::int_t j = 0; j < (mesh.conn).size(1); ++j) {

            EltSizes[i] =  fabs(fabs(EltSizes[i]) - fabs(mesh.Coor(mesh.conn(i,j),0)));

        }

    }


    il::Array<double> Res{(mesh.conn).size(0),0.};

    for (il::int_t k = 0; k < Res.size(); ++k) {

        Res[k] = ((rho_mid[k]*(pow(wh_mid[k],3)))/EltSizes[k])*(1/(12*Visc));

    }

    return Res;
};


// Function that assemble the Finite Difference matrix "L"

il::Array2D<double> Build_L_matrix(Mesh mesh, il::Array2D<double> &d, il::Array2D<double> rho, const int Visc, const double Incr_dil, const double d_wd, const double Init_dil, const double &TimeStep){


    il::Array<double> Kk;
    Kk = ConductivitiesNewtonian(Visc,mesh,rho,d,Incr_dil,d_wd,Init_dil);


    il::Array2D<double> LL{(mesh.conn).size(0)+1,(mesh.conn).size(0)+1,0.};
    il::Array2D<int> ed;
    il::int_t ej;
    il::int_t dofj;
    il::Array<int> t;


    for (il::int_t i = 1; i < (mesh.conn).size(0); ++i) {

        ed = FindPosit_2DArray(mesh.conn,((double)i));

        for (il::int_t j = 0; j < 2; ++j) {

            ej = ed(j,0);
            t = Auxiliary(mesh.conn,ej);

            for (il::int_t k = 0; k < t.size(); ++k) {

                if(t[k]!= i) dofj = t[k];

            }

            LL(i,i) = LL(i,i) - Kk[ej];
            LL(i,dofj) = LL(i,dofj) + Kk[ej];

        }

    }

    // Set the values at the boundary
    LL(0,0) = Kk[0];
    LL(0,1) = Kk[1];

    LL((mesh.conn).size(0)+1,(mesh.conn).size(0)) = Kk[(mesh.conn).size(0)-1];
    LL((mesh.conn).size(0)+1,(mesh.conn).size(0)+1) = Kk[(mesh.conn).size(0)];


    il::Array2D<double> T{(mesh.conn).size(0)+1,(mesh.conn).size(0)+1,TimeStep};

    il::Array2D<double> L{(mesh.conn).size(0)+1,(mesh.conn).size(0)+1,0.};

    for (il::int_t k = 0; k < (mesh.conn).size(0)+1; ++k) {

        for (il::int_t i = 0; i < (mesh.conn).size(0)+1; ++i) {

            L(k,i) = LL(k,i)*T(k,i);

        }

    }

    return L;
};