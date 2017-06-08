//
// Created by DONG LIU on 2/7/17.
//

#include "Viscosityallnodes.h"
#include <il/math.h>
#include<il/norm.h>



//all functions in this cpp needs to consider the change of the size of matrix, that's why for now we don't use the index i0,j0 yet



//the codes invovled with the fluid_parameters density should be corrected, cause the fluid pressure is not along the whole mesh


namespace hfp2d {

    il::Array<double> conductivities_newtonian_open(const il::Array<double> &rho,
                                               const il::Array<double> &vector,
                                               il::Array<double> EltSizes,
                                               Parameters_fluid &fluid_parameters,
                                               double kf, il::io_t) {

        // Inputs:
        //  - rho -> vector of fluid density values at the middle of each element
        //  (size -> Nelts)
        //  - vector -> vector of shear DD or opening DD at the middle of the elements
        //  (size -> Nelts)
        //  - EltSizes -> vector that contains the element sizes
        //  - fluid_parameters -> structure that contains all the fluid parameters
        //  - kf -> fault permeability


        il::Array<double> Res{EltSizes.size(), 0};

        for (il::int_t i = 0; i < Res.size(); ++i) {

            Res[i] = ((rho[i] * (vector[i]*vector[i]*vector[i]* kf)) / EltSizes[i]) *
                     (1 / (12 * fluid_parameters.viscosity));
        }
        return Res;
    }


    il::Array <double> element_size(Mesh mesh_total, const int &p ) {//size N*1
        il::int_t n=mesh_total.nelts();
        il::Array<double> elmt{n,0.};
        for (il::int_t i=0;i<n;i++){
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, int(i), p);
            elmt[i]=segi.size;
        }
        return elmt;
    }


    //Dofp = hfp2d::dofhandle_cg2d(2, mesh_total.nelts(), il::io);
    il::Array<double> average_open(il::Array<double> width,
                                   il::Array2D<int> col_row_i,
                                   il::Array2D<int> col_row_j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim){
        il::int_t n=col_row_j(0,0)-col_row_i(0,0)+1;
        il::Array<double> width_average{n, 0.};
        if(col_row_i(0,1)==0 and col_row_j(0,1)==1)
        {
            for (il::int_t i = 0; i < n; ++i) {
                width_average[i]=0.5*(width[id(i,1)]+width[id(i,dof_dim+1)]);
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==1){
            width_average[0]=0.5*width[1];
            for (il::int_t i = 1; i < n; ++i) {
                width_average[i]=0.5*(width[id(i,1)-dof_dim]+width[id(i,dof_dim+1)-dof_dim]);
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==0){
            width_average[0]=0.5*width[1];
            for (il::int_t i = 1; i < n-1; ++i) {
                width_average[i] = 0.5*(width[id(i,1)-dof_dim]+width[id(i,dof_dim+1)-dof_dim]);
            }
            width_average[n-1]=0.5*width[id(col_row_j(0,0)-col_row_i(0,0),1)-dof_dim];
        }

        if(col_row_i(0,1)==0 and col_row_j(0,1)==0){
            width_average[n-1]=0.5*width[id(col_row_j(0,0)-col_row_i(0,0),1)];

            for (il::int_t i = 0; i < n-1; ++i) {
                width_average[i] = 0.5*(width[id(i,1)]+width[id(i,dof_dim+1)]);
            }
        }
        return width_average;
    }



    il::Array<double> quarter_open(il::Array<double> width,
                                   il::Array2D<int> col_row_i,
                                   il::Array2D<int> col_row_j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   il::Array2D<int> col_matrix){
        // size is 2*N related element number
        // width the width vector with 4 degrees in each element
        il::int_t n=col_row_j(0,0)-col_row_i(0,0)+1;
        il::Array<double> width_quarter{2*n,0.};
        if(col_row_i(0,1)==0 and col_row_j(0,1)==1)
        {
            for (il::int_t i = 0; i < n; ++i) {
                width_quarter[col_matrix(i,0)]=
                        0.75*width[id(i,1)]+0.25*width[id(i,dof_dim+1)];
                width_quarter[col_matrix(i,1)]=
                        0.25*width[id(i,1)]+0.75*width[id(i,dof_dim+1)];
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==1){
            width_quarter[0]=0.25*width[1];
            width_quarter[1]=0.75*width[1];
            for (il::int_t i = 1; i < n; ++i) {
                width_quarter[col_matrix(i,0)]=
                        0.75*width[id(i,1)-dof_dim]+0.25*width[id(i,dof_dim+1)-dof_dim];
                width_quarter[col_matrix(i,1)]=
                        0.25*width[id(i,1)-dof_dim]+0.75*width[id(i,dof_dim+1)-dof_dim];
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==0){

            width_quarter[0]=0.25*width[1];
            width_quarter[1]=0.75*width[1];
            for (il::int_t i = 1; i < n-1; ++i) {
                width_quarter[col_matrix(i,0)] =
                        0.75*width[id(i,1)-dof_dim]+0.25*width[id(i,dof_dim+1)-dof_dim];
                width_quarter[col_matrix(i,1)] =
                        0.25*width[id(i,1)-dof_dim]+0.75*width[id(i,dof_dim+1)-dof_dim];
            }
            width_quarter[col_matrix(n-1,0)]= 0.75*width[id(n-1,1)-dof_dim];
            width_quarter[col_matrix(n-1,1)]=0.25*width[id(n-1,1)-dof_dim];
        }

        if(col_row_i(0,1)==0 and col_row_j(0,1)==0){
            for (il::int_t i = 0; i < n-1; ++i) {
                width_quarter[col_matrix(i,0)] =
                        0.75*width[id(i,1)]+0.25*width[id(i,dof_dim+1)];
                width_quarter[col_matrix(i,1)] =
                        0.25*width[id(i,1)]+0.75*width[id(i,dof_dim+1)];
            }
            width_quarter[col_matrix(n-1,0)]=0.75*width[id(n-1,1)];
            width_quarter[col_matrix(n-1,1)]=0.25*width[id(n-1,1)];
        }
        return width_quarter;
    }


    il::Array2D<double> matrix_lw(il::Array<double> widthB,
                                  il::Array2D<int> col_row_i,il::Array2D<int> col_row_j,
                                  const il::int_t &dof_dim, const int &p,
                                  il::Array<double> element_size_all,
                                  Parameters_fluid &fluid_parameters,
                                  const il::Array2D<int> &id) {

        //widthB should be the ajusted-width profile at the previous moment
        // the size is (N+1)*(N+1)
        //int dof = 2 * (p + 1);
        il::int_t n = col_row_j(0, 0) - col_row_i(0, 0) + 1;//number of related elements

        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements

//        if(col_row_i(0,1)==1){
//            n=n-1; //left node is in the right half of the element
//        };
//        if(col_row_j(0,1)==0){
//            n=n-1;
//        };

        il::Array2D<double> LL{n + 1, n + 1, 0.};//contained node number
        il::Array<double> element_size{col_row_j(0, 0) - col_row_i(0, 0) + 1, 0.};//size of related element number
        il::Array<double> w_mid=average_open(widthB,col_row_i,col_row_j,id,dof_dim);//size of related element number
        for (il::int_t m = 0; m < col_row_j(0, 0) - col_row_i(0, 0) + 1; m++) {
            element_size[m] = element_size_all[m + col_row_i(0, 0)];
        };
        il::Array2D<double> density{col_row_j(0, 0) - col_row_i(0, 0) + 1,2,0.};
        take_submatrix(density,col_row_i(0,0),col_row_j(0,0),0,1,fluid_parameters.density);//this should be size of related element number

        il::Array<double> rho_mid = average(density,il::io);
        il::Array<double> Kk = conductivities_newtonian_open(rho_mid, w_mid,
                                                             element_size,
                                                             fluid_parameters,
                                                             1., il::io);//size of related element number

        il::int_t dofj = 0;
        il::Array2D<int> Dofp = dofhandle_cg2d(2, col_row_j(0, 0) - col_row_i(0, 0) + 1, il::io);


        // Loop over the pressure nodes
        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements

//        for (il::int_t s = 0; s < n + 1; ++s) {
//
//            il::Array2D<int> ed = search(Dofp, int(s+col_row_i(0,1)),
//                                         il::io);//possible rows and coloumns
//
//            for (il::int_t sc = 0; sc < ed.size(0); ++sc) { //number of related elements
//
//                il::int_t ej = ed(sc, 0);//row number-element
//                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
//                if (ec == 0) {//dofj corresponds to the other point of the element
//                    dofj = Dofp(ej, 1)-col_row_i(0,1);
//                } else {
//                    dofj = Dofp(ej, 0)-col_row_i(0,1);
//                }
//
//                LL(s, s) = LL(s, s) - Kk[ej];
//                if(dofj>=0 && dofj<n+1){
//                    LL(s, dofj) = LL(s, dofj) + Kk[ej];
//                }
//            }
//        }

        for (il::int_t s = 0; s < n + 1; ++s) {

            il::Array2D<int> ed = search(Dofp, int(s),
                                         il::io);//possible rows and coloumns

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) { //number of related elements

                il::int_t ej = ed(sc, 0);//row number-element
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {//dofj corresponds to the other point of the element
                    dofj = Dofp(ej, 1);
                } else {
                    dofj = Dofp(ej, 0);
                }

                LL(s, s) = LL(s, s) - Kk[ej];
                if(dofj>=0 && dofj<n+1){
                    LL(s, dofj) = LL(s, dofj) + Kk[ej];
                }
            }
        }
        return LL;
    }//should times the timestep to build the whole matrix


    il::Array2D<double> matrix_vpw(Parameters_fluid &fluid_parameters,
                                   il::Array<double> element_size_all,
                                   il::Array2D<int>col_row_i,
                                   il::Array2D<int> col_row_j,
                                   il::Array<double> widthB,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   const int &p,const il::Array2D<int> &col_matrix,il::io_t) {


        // create the array of element size
        il::int_t n = col_row_j(0, 0) - col_row_i(0, 0) + 1;
        il::Array<double> element_size{n, 0.};
        for (il::int_t m = 0; m < n; m++) {
            element_size[m] = element_size_all[m + col_row_i(0, 0)];
        };

        // Create all the vectors for compressibility of fluid
        // this should be changed with the actual activated collocation points
        // Vector of compressibility of the fluid at nodal points
        il::Array<double> Cf{n+1, fluid_parameters.compressibility};
        // Vector of compressibility of the fluid at the midpoints of each element
        il::Array<double> Cfmid{n, fluid_parameters.compressibility};
        // Vector of compressibility of the fluid at +/- 1/4 of each element
        il::Array<double> Cfquart{2*n, fluid_parameters.compressibility};


        // Assembling the matrix

        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements
//        if(col_row_i(0,1)==1){
//            n=n-1;
//        };
//        if(col_row_j(0,1)==0){
//            n=n-1;
//        };

        il::Array2D<double> Vp{n+1, n+1, 0};
        il::int_t dofj = 0;
        il::Array2D<int> Dofp = dofhandle_cg2d(2, col_row_j(0, 0) - col_row_i(0, 0) + 1, il::io);

        il::Array<double> w_mid=average_open(widthB,col_row_i,col_row_j,id,dof_dim);
        il::Array<double> wquart= quarter_open(widthB,col_row_i,col_row_j, id, dof_dim,col_matrix);

        /// Assembling procedure ///
        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements


//        for (il::int_t m = 0; m < n+1; ++m) {
//
//            il::Array2D<int> ed = search(Dofp, int(m)+col_row_i(0,1), il::io);
//
//            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {//sc-one certain related element number
//
//                il::int_t ej = ed(sc, 0);//row number-element
//                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
//                if (ec == 0) {//dofj corresponds to the other point of the element
//                    dofj = Dofp(ej, 1)-col_row_i(0,1);
//                    Vp(m, m) =Vp(m, m)
//                              +(element_size[ej] / 12) *((widthB[id(ej,1)-dof_dim*col_row_i(0,1)] * Cf[m])
//                                                         + (0.5 * w_mid[ej] * Cfmid[ej])
//                                                         + (3 * wquart[col_matrix(ej,0)] * Cfquart[col_matrix(ej,0)]));
//                    if(dofj>=0 && dofj<n+1){
//                        Vp(m, dofj) = Vp(m, dofj)
//                                      + (element_size[ej] / 12) * ((0.5 * w_mid[ej] * Cfmid[ej])
//                                                                   + (wquart[col_matrix(ej,0)] * Cfquart[col_matrix(ej,0)]));
//                    }
//
//                }
//                else {
//                    dofj = Dofp(ej, 0)-col_row_i(0,1);
//                    Vp(m, m) =Vp(m, m)
//                              +(element_size[ej] / 12) *((widthB[id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)] * Cf[m])
//                                                         + (0.5 * w_mid[ej] * Cfmid[ej])
//                                                         + (3 * wquart[col_matrix(ej,1)] * Cfquart[col_matrix(ej,1)]));
//                    if(dofj>=0 && dofj<n+1){
//                        Vp(m, dofj) = Vp(m, dofj) +
//                                      (element_size[ej] / 12) * ((0.5 * w_mid[ej] * Cfmid[ej])
//                                                                 + (wquart[col_matrix(ej,1)] * Cfquart[col_matrix(ej,1)]));
//                    }
//
//                }
//            }
//        }

        for (il::int_t m = 0; m < n+1; ++m) {

            il::Array2D<int> ed = search(Dofp, int(m), il::io);

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {//sc-one certain related element number

                il::int_t ej = ed(sc, 0);//row number-element
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {//left point
                    dofj = Dofp(ej, 1);
                    if(m==0 && col_row_i(0,1)==1){
                        Vp(m, m) =Vp(m, m)
                                  +(element_size[ej] / 12) *((0.000000000001 * Cf[m])
                                                             + (0.5 * w_mid[ej] * Cfmid[ej])
                                                             + (3 * wquart[col_matrix(ej,0)] * Cfquart[col_matrix(ej,0)]));
                    }
                    else{
                        Vp(m, m) =Vp(m, m)
                                  +(element_size[ej] / 12) *((widthB[id(ej,1)-col_row_i(0,1)*dof_dim] * Cf[m])
                                                             + (0.5 * w_mid[ej] * Cfmid[ej])
                                                             + (3 * wquart[col_matrix(ej,0)] * Cfquart[col_matrix(ej,0)]));
                    }

                    if(dofj>=0 && dofj<n+1){
                        Vp(m, dofj) = Vp(m, dofj)
                                      + (element_size[ej] / 12) * ((0.5 * w_mid[ej] * Cfmid[ej])
                                                                   + (wquart[col_matrix(ej,0)] * Cfquart[col_matrix(ej,0)]));
                    }

                }
                else {//right point
                    dofj = Dofp(ej, 0);
                    if(m==n && col_row_j(0,1)==0.){
                        Vp(m, m) =Vp(m, m)
                                  +(element_size[ej] / 12) *((0.000000000001 * Cf[m])
                                                             + (0.5 * w_mid[ej] * Cfmid[ej])
                                                             + (3 * wquart[col_matrix(ej,1)] * Cfquart[col_matrix(ej,1)]));
                    }
                    else{
                        Vp(m, m) =Vp(m, m)
                                  +(element_size[ej] / 12) *((widthB[id(ej,1+dof_dim)-col_row_i(0,1)*dof_dim] * Cf[m])
                                                             + (0.5 * w_mid[ej] * Cfmid[ej])
                                                             + (3 * wquart[col_matrix(ej,1)] * Cfquart[col_matrix(ej,1)]));
                    }

                    if(dofj>=0 && dofj<n+1){
                        Vp(m, dofj) = Vp(m, dofj) +
                                      (element_size[ej] / 12) * ((0.5 * w_mid[ej] * Cfmid[ej])
                                                                 + (wquart[col_matrix(ej,1)] * Cfquart[col_matrix(ej,1)]));
                    }

                }
            }
        }

        return Vp;
    };

    //build Vw matrix
    il::Array2D<double> matrix_vw(Parameters_fluid &fluid_parameters,
                                  il::Array<double> element_size_all,
                                  il::Array2D<int>col_row_i,
                                  il::Array2D<int> col_row_j,
                                  il::Array<double> widthB,
                                  const il::Array2D<int> &id,
                                  const il::int_t &dof_dim, const int &p,
                                  const il::Array2D<int> &col_matrix,il::io_t) {


        il::int_t n = col_row_j(0, 0) - col_row_i(0, 0) + 1;
        il::Array<double> element_size{n, 0.};
        for (il::int_t m = 0; m < n; m++) {
            element_size[m] = element_size_all[m + col_row_i(0, 0)];
        };
        il::Array2D<double> density{col_row_j(0,0)-col_row_i(0,0)+1,2,0.};
        take_submatrix(density,col_row_i(0,0),col_row_j(0,0),0,1,fluid_parameters.density);//this should be size of related element number??+1??
        il::Array<double> rho_mid=average(density, il::io);;
        il::Array<double> rho_quart=quarter(fluid_parameters.density, il::io);

        // Assembling the matrix

        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements

//        if(col_row_i(0,1)==1){
//            n=n-1;
//        };
//        if(col_row_j(0,1)==0){
//            n=n-1;
//        };

        il::Array2D<double> Vw{n+1, widthB.size(), 0};

        il::Array2D<int> Dofp = dofhandle_cg2d(2, col_row_j(0, 0) - col_row_i(0, 0) + 1, il::io);

        /// Assembling procedure ///

        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements

//        for (il::int_t m = 0; m < n+1; ++m) {
//
//            il::Array2D<int> ed = search(Dofp, int(m)+col_row_i(0,1), il::io);
//
//            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {//sc-one certain related element number
//
//                il::int_t ej = ed(sc, 0);//row number-element
//                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
//                if (ec == 0) {//left point
//                    Vw(m, id(ej,1)-dof_dim*col_row_i(0,1)) =
//                            Vw(m, id(ej,1)-dof_dim*col_row_i(0,1))
//                            +(element_size[ej] / 12) *(density(ej,ec) + 0.5 * rho_mid[ej] +
//                                                        3 * rho_quart[col_matrix(ej,0)]);
//                    if(id(ej,1)-dof_dim*col_row_i(0,1)+dof_dim<widthB.size()){
//                        Vw(m, id(ej,1)-dof_dim*col_row_i(0,1)+dof_dim) =
//                                Vw(m, id(ej,1) -dof_dim*col_row_i(0,1)+dof_dim)
//                                + (element_size[ej] / 12) * (0.5 * rho_mid[ej]
//                                                             + rho_quart[col_matrix(ej,0)]);
//                    }
//
//                }
//                else {//right point
//                    Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)) =
//                            Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1))
//                            +(element_size[ej] / 12) *( density(ej,ec) + 0.5  * rho_mid[ej] +
//                                                        3 * rho_quart[col_matrix(ej,1)]);
//                    if(id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim>=0){
//                        Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim) =
//                                Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim) +
//                                (element_size[ej] / 12) * (0.5  * rho_mid[ej] + rho_quart[col_matrix(ej,1)]);
//                    }
//                }
//            }
//        }

        for (il::int_t m = 0; m < n+1; ++m) {

            il::Array2D<int> ed = search(Dofp, int(m), il::io);

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {//sc-one certain related element number

                il::int_t ej = ed(sc, 0);//row number-element
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {//left point
                    if(id(ej,1)-dof_dim*col_row_i(0,1)<widthB.size() && id(ej,1)-dof_dim*col_row_i(0,1)>=0){
                        Vw(m, id(ej,1)-dof_dim*col_row_i(0,1)) =
                                Vw(m, id(ej,1)-dof_dim*col_row_i(0,1))
                                +(element_size[ej] / 12) *(density(ej,ec) + 0.5 * rho_mid[ej] +
                                                           3 * rho_quart[col_matrix(ej,0)]);
                    }
                    if(id(ej,1)+dof_dim-dof_dim*col_row_i(0,1)<widthB.size() && id(ej,1)+dof_dim-dof_dim*col_row_i(0,1)>=0){
                        Vw(m, id(ej,1)+dof_dim-dof_dim*col_row_i(0,1)) =
                                Vw(m, id(ej,1) +dof_dim-dof_dim*col_row_i(0,1))
                                + (element_size[ej] / 12) * (0.5 * rho_mid[ej]
                                                             + rho_quart[col_matrix(ej,0)]);
                    }

                }
                else {//right point
                    if(id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)<widthB.size() && id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)>=0){
                        Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)) =
                                Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1))
                                +(element_size[ej] / 12) *( density(ej,ec) + 0.5  * rho_mid[ej] +
                                                            3 * rho_quart[col_matrix(ej,1)]);
                    }

                    if(id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim>=0 && id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim<widthB.size()){
                        Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim) =
                                Vw(m, id(ej,1+dof_dim)-dof_dim*col_row_i(0,1)-dof_dim) +
                                (element_size[ej] / 12) * (0.5  * rho_mid[ej] + rho_quart[col_matrix(ej,1)]);
                    }
                }
            }
        }
        return Vw;
    };





    il::int_t find_source(il::int_t i0, il::int_t j0) {
        //we consider i0 and j0 represent the element numbers and i0 and j0 begin from 0.1.2.3.4...
        //and there's a requirement for i0 and j0, j0-i0=2, or j0-i0=1. j0-i0=2 would be better for the paired Mesh
        //returns the node number of injection point

        il::int_t input;

        if (j0 - i0 == 1) {
            input= j0;
        } else {
            input = (j0 - i0+1) / 2+i0;
        }
        return input;
    }



// the matrix for now is built in a not-fully-implicit way, but we can build it in a fully-implicit way

    void construct_matrix_visco(il::Array2D<double> &kmatC,
                              il::Array2D<double> vwc, il::Array2D<double> vp,
                                il::Array2D<double> ll,
                              il::Array<double> cohf, const Material &material,
                              il::Array2D<double> matrix_edge_to_col,
                                il::Array<double> m_source,
                              const Initial_condition &initial_condition,
                              il::Status &status,
                              il::Array<double> pressure_f, il::Array<double> widthB,
                              const int &p, const il::int_t &dof_dim,il::Array<double> delta_p_itera, il::io_t,
                              il::Array<double> &width, il::Array<double> &pressure, il::Array<double> &volume_vary,il::Array<double> & elastic_vary) {
        il::int_t n = kmatC.size(0);
        il::int_t n2= vwc.size(0);
        IL_EXPECT_FAST(n2 == vp.size(0));
        IL_EXPECT_FAST(n2 == ll.size(0));
        IL_EXPECT_FAST(n2 == matrix_edge_to_col.size(1));
        IL_EXPECT_FAST(n == matrix_edge_to_col.size(0));
        IL_EXPECT_FAST(pressure_f.size() == matrix_edge_to_col.size(1));
        il::Array2D<double> newmatrix{n + n2, n + n2, 0.};
        il::Array<double> newvector{n + n2, 0.};


        il::Array<double> fK = il::dot(kmatC, widthB);
        il::Array<double> fP=il::dot(matrix_edge_to_col,pressure_f);
        il::Array<double> conf{n,initial_condition.sigma0};
        il::blas(-1.,conf,1.,il::io,fP);
        il::Array<double> lp_b=il::dot(ll,pressure_f);
        for (int c = 0; dof_dim * c < n; c++) {
            fP[2 * c] = 0;
        }

        for (int m = 0; m < n; ++m) {
            newvector[m]=cohf[m] - fP[m] - fK[m];
            for (int s = 0; s < n; ++s){
                newmatrix(s, m) = kmatC(s, m);
            }
            for (int s1=0; s1<n2;++s1){
                newmatrix(s1+n,m)=vwc(s1,m);
            }
        };
        for (int m1=0;m1<n2;m1++){
            newvector[n+m1]=lp_b[m1]*initial_condition.timestep+m_source[m1]*initial_condition.timestep;
            for (int s2=0;s2<n;++s2){
                newmatrix(s2,m1+n)=matrix_edge_to_col(s2,m1);
            }
            for (int s3=0;s3<n2;s3++){
                newmatrix(s3+n,m1+n)=vp(s3,m1)-ll(s3,m1)*initial_condition.timestep;
            }
        }

//        //add this part to conisder the effect of lag
//        il::Array<double>pressure_transfer=delta_p_itera;
//        il::blas(1.0,pressure_f,1.0,il::io,pressure_transfer);
//        //il::Array<double>pressure_middle=il::dot(matrix_edge_to_col,pressure_transfer);
//        for(int pm=0;pm<pressure_transfer.size();pm++){
//            if(pressure_transfer[pm]<0){
//                for (int pn = 0; pn <n+n2 ; ++pn) {
//                    newmatrix(pn,n+pm)=0.;
//                    newmatrix(n+pm,pn)=0.;
//                };
//                newmatrix(n+pm,n+pm)=1.;
//                newvector[pm]=0.;
//            };
//        };

        il::Array<double> widthinter{n, 0.};
        il::Array<double> pressureinter{n2,0.};
        il::Array<double> dd = il::linear_solve(newmatrix, newvector, il::io,
                                                status);
        il::Array<double> volume_residual=il::dot(newmatrix,dd);
        il::blas(-1.,newvector,1.,il::io,volume_residual);
        il::Array<double> volume_point_residual{n2,0.};
        il::Array<double> elastic_point_residual{n,0.};

        for (int t = 0; t < n; ++t) {
            widthinter[t] = dd[t];
            elastic_point_residual[t]=volume_residual[t];
        };
        for (int t1=0;t1<n2;++t1){
            pressureinter[t1]=dd[n+t1];
            volume_point_residual[t1]=volume_residual[n+t1];
        }
        width = widthinter;
        pressure = pressureinter;
        volume_vary=volume_point_residual;
        elastic_vary=elastic_point_residual;

        status.abort_on_error();
    }

    void plasticity_loop_visco(const Material &material,
                               const Initial_condition &initial_condition,
                               il::int_t c_i, il::int_t c_j,
                               const il::Array2D<double> &kmat,
                               const il::Array2D<int> &id,const int &p,
                               il::Array<double> widthB, il::Array<double> pressure_f,
                               il::Array<double> width_history,
                               const il::int_t &dof_dim,
                               const il::Array2D<int>&col_matrix,
                               Parameters_fluid fluid_parameters,
                               il::Array<double> element_size_all,
                               const il::Array2D<double> &matrix_edge_to_col_all,
                               il::int_t input,il::io_t,
                               il::Array<double> &delta_width,
                               il::Array <double> &pressure_change,
                               il::Array<double> &coht,int &mm,il::Array<double> &volume_vary,il::Array<double> &elastic_vary) {
        il::Array2D<int> col_row_i=search(col_matrix,int(c_i),il::io);
        il::Array2D<int> col_row_j=search(col_matrix,int(c_j),il::io);

        il::Array2D<double> kmatnew{dof_dim * p*(col_row_j(0,0)-col_row_i(0,0) + 1), dof_dim * p*(col_row_j(0,0)-col_row_i(0,0) + 1), 0.};
        take_submatrix(kmatnew, id(col_row_i(0,0),1),
                       id(col_row_j(0,0),1+dof_dim),
                       id(col_row_i(0,0),1),
                       id(col_row_j(0,0),1+dof_dim),
                       kmat);

        il::int_t n = col_row_j(0, 0) - col_row_i(0, 0) + 1;
        il::int_t left_node=col_row_i(0,0);
        il::int_t right_node=col_row_j(0,0)+1;


        // This part to consider when only one collocation point
        // in the element is activated, the pressure is only considered
        // till the inner activated node of the crack
        // However, pressure nodes should involve all the nodes in the activated elements

//        if(col_row_i(0,1)==1){
//            n=n-1;
//            left_node=col_row_i(0,0)+1;
//        }
//        if(col_row_j(0,1)==0){
//            n=n-1;
//            right_node=col_row_j(0,0);
//        };

        il::Array2D<double> vwc{n+1,dof_dim * (c_j - c_i + 1), 0.};
        vwc=matrix_vw(fluid_parameters, element_size_all, col_row_i,col_row_j,
                      widthB, id, dof_dim, p,col_matrix,il::io);
        il::Array2D<double> ll{n+1,n+1,0.};
        il::Array2D<double> vp{n+1,n+1,0.};
        il::Array2D<double> matrix_edge_to_col{id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))+1,right_node-left_node+1,0.};
        take_submatrix(matrix_edge_to_col,id(col_row_i(0,0),dof_dim*col_row_i(0,1)),
                id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1),
                       int(left_node),int(right_node),matrix_edge_to_col_all);
        il::Array<double> m_source{n+1,0.};
        m_source[input-left_node]=initial_condition.Q0*fluid_parameters.density(input,0);

        il::Array<double> delta_w_ini{widthB.size(), 0.};

        il::Array<double> coh{widthB.size(), 0.};

        int k = 0;
        double error_w = 1.;
        il::Array<il::Status> statusk{200,};
        double error_p=1.;
        //double error_r=1;

        il::Array<double> delta_w_bitera=delta_w_ini;
        il::Array<double> delta_p_ini{n+1,0.};
        il::Array<double> delta_p_itera=delta_p_ini;

        il::Array<double> volume_vary_inter{n+1,0.};
        il::Array<double> elastic_vary_inter{widthB.size(),0.};


        while (k < 200 && (error_w > 0.00001 or error_p>0.00001)) {
            il::Array<double> width_inm = widthB;
            il::blas(1.0, delta_w_ini, 1.0, il::io, width_inm);
            //Dugdale cohesive model
                       coh = cohesive_force_col(material, width_inm, col_row_i,
                                     col_row_j, width_history, p,dof_dim,id);

            //Linear Softening model
//            coh = cohesive_force_linear(material, width_inm, col_row_i,
//            col_row_j, width_history, p,dof_dim,id);

            ll=matrix_lw(width_inm, col_row_i,col_row_j,dof_dim, p,
                         element_size_all, fluid_parameters, id);
                   vp=matrix_vpw(fluid_parameters,element_size_all,col_row_i,col_row_j,
                          width_inm, id, dof_dim, p,col_matrix,il::io);


            construct_matrix_visco(kmatnew, vwc, vp,ll,coh, material,
                                   matrix_edge_to_col,m_source,
                                 initial_condition, statusk[k],
                                 pressure_f, widthB, p, dof_dim,delta_p_itera,il::io,
                                 delta_width, pressure_change,volume_vary_inter,elastic_vary_inter);//iteration pressure to eliminate r&c
            il::Array<double> intermediate = delta_width;
            il::blas(-1.0, delta_w_bitera, 1.0, il::io, intermediate);
            double error_inter = il::norm(intermediate, il::Norm::L2);
            error_w = error_inter / il::norm(delta_width, il::Norm::L2);
            il::blas(0.85, delta_width, 0.15, il::io, delta_w_ini);//iteration condition
            delta_w_bitera=delta_width;

            il::Array<double> inter_pressure=pressure_change;
            il::blas(-1.0,delta_p_ini,1.0,il::io,inter_pressure);
            double error_inter_p=il::norm(inter_pressure,il::Norm::L2);
            error_p=error_inter_p/il::norm(pressure_change,il::Norm::L2);
//            error_p=0.0;//this line tries to make the pressure iteration criteria unavalible

//            delta_p_itera=delta_p_ini;
//            il::blas(0.85, pressure_change, 0.15, il::io, delta_p_itera);//iteration pressure to eliminate rows and columns






            delta_p_ini=pressure_change;


//            error_r=fabs(il::dot(vwc,delta_width)
// +pressure_change*material.U-initial_condition.Q0*initial_condition.timestep);

            ++k;

        }
        mm=k;
        //from here you can consider to plot the current coh vector

        for (il::int_t r = 0; r < widthB.size(); r++) {
            coht[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+ r] = coh[r];
        }

        volume_vary =volume_vary_inter;
        elastic_vary=elastic_vary_inter;

    }


    void
    propagation_loop_visco(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                         const Material &material,
                         const Initial_condition &initial_condition,
                         il::int_t i0, il::int_t j0, int nstep,
                         il::Status &status, Parameters_fluid fluid_parameters, il::io_t,
                         il::Array2C<double> &widthlist,
                         il::Array2D<double> &plist, il::Array<double> &l_coh,
                         il::Array<double> &l, il::Array2C<double> &coh_list,
                         il::Array<int> &mvalue,int &break_time,
                         il::Array2C<double> &stress_list, il::Array<double> &energy_g,
                           il::Array2D<double> &volume_vary_list,il::Array2D<double> &elastic_vary_list) {

        int n = mesh_total.nelts();
        IL_EXPECT_FAST(n == id.size(0));
        il::Array<double> element_size_all=element_size(mesh_total,p);

        il::Array2D<int> col_matrix=collocation_matrix(n,p);
        il::int_t dof_dim=id.size(1)/2;
        il::int_t dof = dof_dim * (p + 1);
        il::int_t ttn=id.size(0)*id.size(1);
        il::Array2D<double> kmat{ttn,ttn,0.};
        kmat=basic_assembly(mesh_total, id, p, material.Ep);
        il::Array2D<int> Dofp = dofhandle_cg2d(int(dof_dim),n,il::io);
        il::Array2D<double> matrix_edge_to_col_all = from_edge_to_col_cg(int(dof_dim),
        id,Dofp, il::io);

        il::Array<double> width_history{ttn, 0.};

        //i0 and j0 are the beginning element number in C++, starting from 0,1,2.....
        for (il::int_t h = i0; h < j0 + 1; ++h) {
            width_history[id(h,dof_dim*0+1)] = 1.;
            width_history[id(h,dof_dim*1+1)] = 1.;
        };



        il::Array<double> widthB = initialwidth_col(kmat, id, p,
                                                    initial_condition, i0, j0,
                                                    status,dof_dim);
        il::Array<double> width;
        il::Array2C<double> widthlistinter{nstep, ttn,0.};
        il::Array2C<double> stress_profile{nstep,ttn,0.};
        il::Array2D<double> plistinter{nstep + 1, ttn,0.};
        il::Array<double> pressure;
        il::Array<double> pressure_change;
        il::Array<double> pressure_b{j0-i0+1+1,initial_condition.pini};
        //plistinter[0] = initial_condition.pini;

        il::Array<int> m_value{nstep,0};
        int mk=0;

        int s = 0;
        int m = 0;
        il::int_t i = i0;
        il::int_t j = j0;
        il::int_t input=find_source(i0,j0);
        il::int_t left_node=i0;
        //il::int_t right_node=j0+1;

        double time = 0;
        bool break_value=false;
        bool break_value_2=false;
        il::Array<double> crack_length{nstep+1, 0.};
        il::Array<double> length_coh{nstep, 0.};
        il::Array<double> energy_j{nstep, 0.};

        il::Array2C<double> coh_tt{nstep, ttn, 0.};
        il::Array<double> coht{ttn, 0.};

        il::Array2D<double> volume_vary_matrix{nstep,n+1,0.};//to output the residual of the volume control equation
        il::Array2D<double> elastic_vary_matrix{nstep,ttn,0.};

        il::int_t c_i =col_matrix(i0,0);
        il::int_t c_j =col_matrix(j0,1);
        double l0=0.;
        for(il::int_t ei=i0;ei<j0+1;ei++){
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, int(ei), p);
            l0+=segi.size;
        }
        crack_length[0]=l0;

        break_time=nstep;


        while (s < nstep) {

            while (m < 10000) {//the stop limit value for m should be big enough
                if (break_value) {
                    break_value_2=true;
                    break;
                };
                il::Array<double> delta_width;
                il::Array<double> volume_vary;
                il::Array<double> elastic_vary;

                plasticity_loop_visco(material, initial_condition,c_i, c_j,
                                      kmat, id, p, widthB, pressure_b,
                                      width_history,dof_dim,col_matrix,
                                      fluid_parameters, element_size_all,
                                      matrix_edge_to_col_all,input,il::io,
                                      delta_width, pressure_change, coht,mk,volume_vary,elastic_vary);

                m_value[s]=mk;

                il::Array<il::int_t> index = stress_criteria_col(kmat, id, material,
                                                                 initial_condition,
                                                                 widthB, delta_width,
                                                                 i, j, c_i,c_j,
                                                                 dof_dim,col_matrix,p);

                if (c_i == index[0] and c_j == index[1]) {
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    pressure=pressure_b;
                    il::blas(1.0,pressure_change,1.0,il::io,pressure);
                    break;
                }
                if (index[0] <= 0 or index[1] >= dof_dim*n - 1) {//not sure the 2 here is (p+1) or dof_dim
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    pressure=pressure_b;
                    il::blas(1.0,pressure_change,1.0,il::io,pressure);

                    break_time=s+1;//because s starts from 0
                    break_value=true;
                    break;
                }
                il::Array<double> widthBnew{dof_dim * (index[1] - index[0] + 1),
                                            0.};
                for (int wbn = 0; wbn < widthB.size(); wbn += 1) {
                    widthBnew[dof_dim * (c_i - index[0]) + wbn] = widthB[wbn];
                }
                il::Array<double> delta_width_new{
                        dof_dim * (index[1] - index[0] + 1), 0.};
                for(int ela_i=0;ela_i<elastic_vary_matrix.size(1);ela_i++){
                    elastic_vary_matrix(s,ela_i)=0.;
                }


                for (int dwn = 0; dwn < widthB.size(); dwn += 1) {
                    delta_width_new[dof_dim * (c_i - index[0]) +
                                    dwn] = delta_width[dwn];
                    elastic_vary_matrix(s,dof_dim*c_i+dwn)=elastic_vary[dwn];
                }
                widthB.resize(dof_dim * (index[1] - index[0] + 1));
                widthB = widthBnew;
                width = widthBnew;

                il::Array2D<int> col_row_i=search(col_matrix,int(index[0]),il::io);
                il::Array2D<int> col_row_j=search(col_matrix,int(index[1]),il::io);
                left_node=col_row_i(0,0);
                //right_node=col_row_j(0,1)+1;
                il::int_t pn=col_row_j(0,0)-col_row_i(0,0)+1+1;


                // This part to consider when only one collocation point
                // in the element is activated, the pressure is only considered
                // till the inner activated node of the crack
                // However, pressure nodes should involve all the nodes in the activated elements

//                if(col_row_i(0,1)==1){
//                    pn=pn-1;
//                    left_node=col_row_i(0,0)+1;
//                };
//                if(col_row_j(0,1)==0){
//                    pn=pn-1;
//                    //right_node=col_row_j(0,0);
//                };

                il::Array2D<int> col_row_i_b=search(col_matrix,int(c_i),il::io);
                il::Array2D<int> col_row_j_b=search(col_matrix,int(c_j),il::io);
                il::int_t left_node_b=col_row_i_b(0,0);
                //il::int_t right_node_b=col_row_j_b(0,1)+1;
                il::int_t pn_b=col_row_j_b(0,0)-col_row_i_b(0,0)+1+1;

                // This part to consider when only one collocation point
                // in the element is activated, the pressure is only considered
                // till the inner activated node of the crack
                // However, pressure nodes should involve all the nodes in the activated elements

//                if(col_row_i_b(0,1)==1){
//                    pn_b=pn_b-1;
//                    left_node_b=col_row_i_b(0,0)+1;
//                };
//                if(col_row_j_b(0,1)==0){
//                    pn_b=pn_b-1;
//                    //right_node_b=col_row_j_b(0,0);
//                };
                il::Array<double> pressure_b_new{pn,0.};
                il::Array<double> pressure_change_new{pn,0.};

                for(int v_clear=0;v_clear<n+1;v_clear++){
                    volume_vary_matrix(s,v_clear)=0.;
                };

                for(int pb_i=0;pb_i<pn_b;pb_i++){
                    pressure_b_new[pb_i-left_node+left_node_b]=pressure_b[pb_i];
                    pressure_change_new[pb_i-left_node+left_node_b]=pressure_change[pb_i];
                    volume_vary_matrix(s,pb_i+left_node_b)=volume_vary[pb_i];
                }
                pressure_b.resize(pn);
                pressure_b = pressure_b_new;
                pressure = pressure_b_new;


                c_i = index[0];
                c_j = index[1];
                i=index[2];
                j=index[3];


                il::blas(1., delta_width_new, 1.0, il::io, width);
                il::blas(1.,pressure_change_new,1.0,il::io,pressure);
                ++m;
            }
            if(break_value_2){
                break;
            };



            il::Array<double> width_large{dof * n, 0.};
            for (int r = 0; r < width.size(); ++r) {
                width_large[dof_dim * c_i + r] = width[r];
            }//needs to be carefully checked
            for (int sc = 0; sc < dof * n; ++sc) {
                widthlistinter(s, sc) = width_large[sc];
            }
            il::Array<double> pressure_large{n+1,0.};
            for(int sp=0;sp<pressure.size();sp++){
                plistinter(s,left_node+sp)=pressure[sp];
                pressure_large[left_node+sp]=pressure[sp];
            }


            cc_length_col(mesh_total, material, width_large, i, j, c_i,c_j,
                          width_history, p,id,col_matrix,dof_dim, il::io,
                          length_coh[s], crack_length[s+1],energy_j[s]);

            width_history = write_history_col(width, width_history, c_i,c_j,
                                              col_matrix, material,p,dof_dim,id);
            //Here we use width_large, could be simpler for the coding

            for (int cohi = 0; cohi < dof * n; cohi++) {
                coh_tt(s, cohi) = coht[cohi];
            }

            il::Array<double> pressure_current=il::dot(matrix_edge_to_col_all,pressure_large);
            il::Array<double> stress_current=il::dot(kmat,width_large);
            for (int stre=0;stre<dof*n;stre++){
                stress_profile(s,stre)=stress_current[stre];//+pressure_current[stre];
            }

            widthB = width;
            pressure_b=pressure;
            time += initial_condition.timestep;
            ++s;

        }
        widthlist = widthlistinter;
        plist = plistinter;
        l = crack_length;
        l_coh = length_coh;
        coh_list = coh_tt;
        energy_g=energy_j;
        mvalue=m_value;//to output the iteration times in plasticity_loop to get the solution
        stress_list=stress_profile;
        volume_vary_list=volume_vary_matrix;
        elastic_vary_list=elastic_vary_matrix;
        status.abort_on_error();
    }

}