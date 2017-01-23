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
// It is not completely general -> the output is always a matrix (2x2)
il::Array2D<int> FindPosit_2DArray(il::Array2D<int> &arr2D, double seek) {

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

// Function to find out the position of a value in a 2D array
// It returns 2x2 array with row&col of the seek value
// It is completely general in a sense that the output can be a vector or a matrix (2x2)
il::Array2D<int> Position_2DArray(il::Array2D<int> &arr2D, double seek) {

    il::Array2D<il::int_t> M{2 * arr2D.size(0), 2, -1};
    il::int_t k = 0;

    for (il::int_t i = 0; i < arr2D.size(0); ++i) {

        for (il::int_t j = 0; j < arr2D.size(1); ++j) {

            if (arr2D(i, j) == seek) {

                M(k, 0) = i;
                M(k, 1) = j;

                k = k + 1;

            }

        }
    }

    il::Array2D<int> outp{k, 2, 0};

    for (il::int_t l = 0; l < k; ++l) {

        for (il::int_t j = 0; j < 2; ++j) {

            outp(l, j) = M(l, j);

        }

    }

    return outp;
};


// Auxiliary function for assembly process
// It returns a given row (vector - specified by idx) of a 2D array
il::Array<int> Auxiliary(il::Array2D<int> &arr, il::int_t idx){


    il::Array<int> vect{arr.size(1),0.};

    for (il::int_t i = 0; i < vect.size(); ++i) {

        vect[i] = arr(idx,i);

    }

    return vect;

};


//// FVM routines ////

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
    il::int_t ni;
    il::int_t ej;
    il::int_t dofj;
    il::Array<int> t;

    // Loop over all the "inner" nodes (the boundary nodes are not included)
    for (il::int_t i = 0; i < LL.size(0); ++i) {

        ed = Position_2DArray(mesh.conn,((double)i));
        ni = ed.size(0);

        for (il::int_t j = 0; j < ni; ++j) {

            ej = ed(j,0);
            t = Auxiliary(mesh.conn,ej);

            for (il::int_t k = 0; k < t.size(); ++k) {

                if(t[k]!= i) dofj = t[k];

            }

            LL(i,i) = LL(i,i) - Kk[ej];
            LL(i,dofj) = LL(i,dofj) + Kk[ej];

        }

    }

    il::Array2D<double> T{(mesh.conn).size(0)+1,(mesh.conn).size(0)+1,TimeStep};

    il::Array2D<double> L{(mesh.conn).size(0)+1,(mesh.conn).size(0)+1,0.};

    for (il::int_t k = 0; k < (mesh.conn).size(0)+1; ++k) {

        for (il::int_t i = 0; i < (mesh.conn).size(0)+1; ++i) {

            L(k,i) = LL(k,i)*T(k,i);

        }

    }

    return L;
};


// Function for assembling the Pressure matrix "Vp"
il::Array2D<double> BuildVpMatrix(Mesh mesh, const double Incr_dil, const double Init_dil, const double CompressFluid, il::Array2D<double> &d, const double d_wd){


    //Create an auxiliary vector for the assembling
    il::Array2D<int> h{2*(mesh.conn).size(0) + 1, 2, 0};

    for (il::int_t i = 0; i < h.size(0); ++i) {

        h(i,0) = i;
        h(i,1) = i+1;

    }

    // mesh nodes {0,1,2,3..}
    il::Array<int> vertices{(mesh.conn).size(0)+1,0};
    for (il::int_t j = 0; j < vertices.size(); ++j) {

        vertices[j] = j;

    }

    //create the array of element size
    il::Array<double> EltSizes{(mesh.conn).size(0), 0.};

    for (il::int_t i = 0; i < EltSizes.size(); ++i) {

        for (il::int_t j = 0; j < (mesh.conn).size(1); ++j) {

            EltSizes[i] =  fabs(fabs(EltSizes[i]) - fabs(mesh.Coor(mesh.conn(i,j),0)));

        }

    }

    // Create the all the vectors for compressibility of fluid
    il::Array<double> Cf{vertices.size(),CompressFluid}; // Vector of compressibility of the fluid at nodal points
    il::Array<double> Cfmid{(mesh.conn).size(0),CompressFluid}; // Vector of compressibility of the fluid at the midpoints of each element
    il::Array<double> Cfquart{2*(mesh.conn).size(0),CompressFluid}; // Vector of compressibility of the fluid at +/- 1/4 of each element


    il::Array<double> d_left{(mesh.conn).size(0),0.};
    for (il::int_t k = 0; k < d_left.size(); ++k) {

        d_left[k] = d(k,0);
    }

    il::Array<double> d_right{(mesh.conn).size(0),0.};
    for (il::int_t j = 0; j < d_right.size(); ++j) {

        d_right[j] = d(j,1);
    }

    il::Array<double> whi_left{(mesh.conn).size(0),0.} , whi_right{(mesh.conn).size(0),0.};

    whi_left = Dilatancy(Init_dil,Incr_dil,d_wd,d_left);
    whi_right = Dilatancy(Init_dil,Incr_dil,d_wd,d_right);

    il::Array2D<double> whi{(mesh.conn).size(0),2,0.};

    for (il::int_t l = 0; l < whi.size(0); ++l) {

        whi(l,0) = whi_left[l];
        whi(l,1) = whi_right[l];

    }

    il::Array<double> whi_mid{(mesh.conn).size(0),0.};
    whi_mid = Average(whi);

    il::Array<double> wquart{2*(mesh.conn).size(0),0.};
    wquart = Quarter(whi);


    // Assembling the matrix
    il::Array2D<double> Vp{(mesh.conn).size(0)+1, (mesh.conn).size(0)+1, 0.};
    il::Array2D<int> ed;
    il::Array2D<int> hi;
    il::int_t ni;
    il::int_t ej;
    il::int_t hj;
    il::int_t dofj;
    il::Array<int> t;

    il::Array<double> Whi{2*(mesh.conn).size(0),0.}; //Vector that we need for assembling process
    for (il::int_t n = 0, q = 0; n < whi_left.size(); ++n, q = q+2) {

        Whi[q] = whi_left[n];
        Whi[q+1] = whi_right[n];

    }


    for (il::int_t m = 0; m < Vp.size(0); ++m) {

        ed = Position_2DArray(mesh.conn,((double)m));
        hi = Position_2DArray(h,((double)(2*m)));
        ni = ed.size(0);

        for (il::int_t i = 0; i < ni; ++i) {

            ej = ed(i,0);
            hj = hi(i,0);
            t = Auxiliary(mesh.conn,ej);

            for (il::int_t j = 0; j < t.size(); ++j) {

                if(t[j]!= vertices[m]) dofj = t[j];

            }

            Vp(vertices[m],vertices[m]) = Vp(vertices[m],vertices[m]) + (EltSizes[ej]/12)*((Whi[hj]*Cf[vertices[m]]) + (0.5*whi_mid[ej]*Cfmid[ej]) + (3*wquart[hj]*Cfquart[hj]));
            Vp(vertices[m], dofj) = Vp(vertices[m],dofj) + (EltSizes[ej]/12)*((0.5*whi_mid[ej]*Cfmid[ej]) + (wquart[hj]*Cfquart[hj]));

        }

    }

    return Vp;

};


// Function for assembling the Mass matrix "Vd"
il::Array2D<double> BuildVdMatrix(Mesh mesh, const double Incr_dil, const double d_wd, il::Array2D<int> id, il::Array2D<double> rho, il::Array2D<double> &d){


    //Create an auxiliary vector for the assembling
    il::Array2D<int> h{2*(mesh.conn).size(0) + 1, 2, 0};

    for (il::int_t i = 0; i < h.size(0); ++i) {

        h(i,0) = i;
        h(i,1) = i+1;

    }

    // mesh nodes {0,1,2,3..}
    il::Array<int> vertices{(mesh.conn).size(0)+1,0};
    for (il::int_t j = 0; j < vertices.size(); ++j) {

        vertices[j] = j;

    }

    //create the array of element size
    il::Array<double> EltSizes{(mesh.conn).size(0), 0.};

    for (il::int_t i = 0; i < EltSizes.size(); ++i) {

        for (il::int_t j = 0; j < (mesh.conn).size(1); ++j) {

            EltSizes[i] =  fabs(fabs(EltSizes[i]) - fabs(mesh.Coor(mesh.conn(i,j),0)));

        }

    }


    // Create the matrix of dof for just shear DD
    // Remember: 4 DOFs per element {shear_i, normal_i, shear_j, normal_j}
    il::Array2D<int> dofw{(mesh.conn).size(0), 2,0.};

    for (il::int_t k = 0; k < dofw.size(0); ++k) {

        for (il::int_t i = 0, q=0; i < dofw.size(1); ++i, q=q+2) {

            dofw(k,i) = id(k,q);

        }

    }


    il::Array<double> rho_mid{(mesh.conn).size(0), 0.};
    il::Array<double> rho_quart{2*(mesh.conn).size(0),0.};

    rho_mid = Average(rho);
    rho_quart = Quarter(rho);

    il::Array<double> d_left{(mesh.conn).size(0),0.};
    for (il::int_t k = 0; k < d_left.size(); ++k) {

        d_left[k] = d(k,0);
    }

    il::Array<double> d_right{(mesh.conn).size(0),0.};
    for (il::int_t j = 0; j < d_right.size(); ++j) {

        d_right[j] = d(j,1);
    }

    il::Array<double> Bi_left{(mesh.conn).size(0),0.} , Bi_right{(mesh.conn).size(0),0.};

    Bi_left = DDilatancy(Incr_dil,d_wd,d_left);
    Bi_right = DDilatancy(Incr_dil, d_wd, d_right);

    il::Array2D<double> Bi{(mesh.conn).size(0),2,0.};

    for (il::int_t l = 0; l < Bi.size(0); ++l) {

        Bi(l,0) = Bi_left[l];
        Bi(l,1) = Bi_right[l];

    }

    il::Array<double> Bi_mid{(mesh.conn).size(0),0.};
    Bi_mid = Average(Bi);

    il::Array<double> Bi_quart{2*(mesh.conn).size(0),0.};
    Bi_quart = Quarter(Bi);


    // Assembling the matrix

    il::Array2D<double> Vd{(mesh.conn).size(0)+1, 4*((mesh.conn).size(0)), 0.};
    il::Array2D<int> ed;
    il::Array2D<int> hi;
    il::int_t ni;
    il::int_t ej;
    il::int_t hj;
    il::Array<int> t;

    il::int_t dofwi;
    il::int_t dofwj;


    il::Array<double> Rho{2*(mesh.conn).size(0),0.}; // We need for the assembling
    il::Array<double> BI{2*(mesh.conn).size(0),0.}; // We need for the assembling

    for (il::int_t m = 0, q=0; m < rho.size(0); ++m, q=q+2) {

            Rho[q] = rho(m,0);
            Rho[q+1] = rho(m,1);

    }

    for (il::int_t m = 0, q=0; m < Bi.size(0); ++m, q=q+2) {

        BI[q] = Bi(m,0);
        BI[q+1] = Bi(m,1);

    }


    for (il::int_t i = 0; i < Vd.size(0); ++i) {

        ed = Position_2DArray(mesh.conn,((double)i));
        hi = Position_2DArray(h, ((double)(2*i)));
        ni = ed.size(0);

        for (il::int_t j = 0; j < ni; ++j) {

            ej = ed(j,0);
            hj = hi(j,0);

            dofwi = dofw(ej,ed(j,1));

            t = Auxiliary(dofw,ej);

            for (il::int_t k = 0; k < t.size(); ++k) {

                if (t[k] != dofwi) dofwj = t[k];

            }

            Vd(vertices[i],dofwi) = Vd(vertices[i],dofwi) + (EltSizes[ej]/12)*((Rho[hj]*BI[hj]) + (0.5*rho_mid[ej]*Bi_mid[ej]) + (3*rho_quart[hj]*Bi_quart[hj]));
            Vd(vertices[i],dofwj) = Vd(vertices[i],dofwj) + (EltSizes[ej]/12)*((0.5*rho_mid[ej]*Bi_mid[ej]) + (rho_quart[hj]*Bi_quart[hj]));

        }

    }

    // Taking finally only the shear related contributions...
    il::Array2D<double> VD{Vd.size(0),Vd.size(1)/2,0.};

    for (il::int_t n = 0; n < VD.size(0); ++n) {

        for (il::int_t i = 0, q=0; i < Vd.size(1); ++i, q = q+2) {

            VD(n,i) = Vd(n,q);

        }

    }

    return VD;

};