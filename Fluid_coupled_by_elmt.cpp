//
// Created by DONG LIU on 7/25/17.
//

// Taking the fluid viscosity into account, this code
// can consider the fluid compressibility and a non-uniform fluid density
// We try to do the propagation regime of element by element

// Different CZMs(linear-softening and three exponential cohesive models)
// can be considered if changing part of the function
// cohesive_force_linear_elemt

//We have the timestep adaptive method which can control the number of elements
//propagated per time step
//calculate the velocity each time and determines the next timestep
//to get better convergence


// we can change the permeability of the fracture tip with this code

#include "Fluid_coupled_by_elmt.h"
#include <il/math.h>
#include<il/norm.h>


namespace hfp2d {

    il::Array2D<double> width_ec(const il::int_t &dof_dim, const il::Array2D<int> &id){
        il::int_t ttn=id.size(0)*id.size(1);
        il::Array2D<double> width_edge_to_col{ttn,ttn,0.};
        for(il::int_t m=0;m<ttn;m+=id.size(1)){
            // including the shear displacement arrangement, pay attention this can't be
            // simply replaced by the EC matrix, because this is 4Nx4N while EC
            // matrix is 4N*(N+1)
            width_edge_to_col(m,m)=sqrt(2.)/2.;
            width_edge_to_col(m,m+dof_dim)=1.-sqrt(2.)/2.;
            width_edge_to_col(m+1,m+1)=sqrt(2.)/2.;
            width_edge_to_col(m+1,m+1+dof_dim)=1.-sqrt(2.)/2.;
            width_edge_to_col(m+dof_dim,m)=1.-sqrt(2.)/2.;
            width_edge_to_col(m+dof_dim,m+dof_dim)=sqrt(2.)/2.;
            width_edge_to_col(m+1+dof_dim,m+1)=1.-sqrt(2.)/2.;
            width_edge_to_col(m+1+dof_dim,m+1+dof_dim)=sqrt(2.)/2.;
        }
        return width_edge_to_col;
    };

    void cc_length_elemt(Mesh mesh_total, const Material &material,
                       il::Array<double> width_large,
                       il::int_t i, il::int_t j,
                       il::Array<double> width_history, const int &p,
                       const il::Array2D<int> &id,
                       const il::Array2D<int> &col_matrix,
                       const il::int_t &dof_dim,il::io_t,
                       double &length_coh, double &crack_length) {

        crack_length = 0.;
        length_coh = 0.;

        //calculation of crack length - sum of the activaed element length
        for (il::int_t s = 0; s < j - i + 1; ++s) {
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, int(s + i), p);
            crack_length += segi.size;
        }


        //calculation of cohesive length
        for(il::int_t s1=i;s1<j+1;++s1){//s1 is the element number, starting from 0
            SegmentCharacteristic segi1 = get_segment_DD_characteristic(
                    mesh_total, int(s1), p);
            if ((width_large[id(s1,1)]>=0
                 and width_large[id(s1,1)] < material.wc
                 and width_history[id(s1,1)] < material.wc)
                or (width_large[id(s1,dof_dim+1)]>=0
                    and width_large[id(s1,dof_dim+1)] < material.wc
                    and width_history[id(s1,dof_dim+1)] < material.wc)) {
                // if we don't apply cohesive force for the negative opening part,
                // we don't count this part as cohesive length
                if (width_large[id(s1,1)] < material.wc
                    and width_large[id(s1,dof_dim+1)] < material.wc) {
                    length_coh += segi1.size;
                }
                else {
                    double small;
                    if (width_large[id(s1,1)] > width_large[id(s1,dof_dim+1)]) {
                        small = width_large[id(s1,dof_dim+1)];
                    }
                    else {
                        small = width_large[id(s1,1)];
                    };
                    length_coh += segi1.size * (material.wc - small) /
                                  fabs(width_large[id(s1,dof_dim+1)] - width_large[id(s1,1)]);
                }
            };
        };
    }

    il::Array<double> write_history_elemt(il::Array<double> width,
                                          il::Array<double> width_history,
                                          il::int_t i, il::int_t j,
                                          const Material &material,
                                          const int &p, const il::int_t &dof_dim,
                                          const il::Array2D<int> &id,
                                          il::Array2D<double> width_e_to_c) {

        il::int_t n = width.size();
        il::int_t n_t = width_history.size();
        il::Array<double> width_history_new{n_t, 0.};
        il::Array<double> width_new{n, 0.};

        il::Array2D<double>width_edge_col{id(j,1+dof_dim)-id(i,0)+1,id(j,1+dof_dim)-id(i,0)+1,0.};
        take_submatrix(width_edge_col,id(i,0),id(j,1+dof_dim),id(i,0),id(j,1+dof_dim),width_e_to_c);
        width_new = il::dot(width_edge_col,width);
        //changes the 7th April, m initial value from 0 to 1 a
        // nd the varying step from 1 to 2,
        // which makes the width_history only relates to the opening not the slip.
        for (il::int_t m = 0; m < n; m += 1) {
            width_history_new[id(i,0) + m] = width_history[id(i,0)+ m];
            if (width_new[m] > width_history[id(i,0) + m]) {
                width_history_new[id(i,0)+ m] = width_new[m];
            }
        }
        return width_history_new;
    }




    il::Array<double> cohesive_force_linear_elemt(const Material &material,
                                                  il::Array<double> width,
                                                  il::int_t i, il::int_t j,
                                                  il::Array<double> width_history,
                                                  const int &p,
                                                  const il::int_t &dof_dim,
                                                  const il::Array2D<int> &id,
                                                  il::Array2D<double> width_e_to_c) {
        il::Array<double> f{width.size(), 0.};
        //we switch the dislocation from edge point to collocation point

        il::Array2D<double>width_edge_col{id(j,1+dof_dim)-id(i,0)+1,id(j,1+dof_dim)-id(i,0)+1,0.};
        take_submatrix(width_edge_col,id(i,0),id(j,1+dof_dim),id(i,0),id(j,1+dof_dim),width_e_to_c);
        il::Array<double> width_etoc=il::dot(width_edge_col,width);

        for (il::int_t s = 1; s < width.size(); s += dof_dim) {
            // the linear-softening model and the partly-exponenetial model
            // considers the opening history

//            if (width_etoc[s]>0 and width_etoc[s] < material.wc and
//                width_history[id(i,0)+s] == 0.) {
//                //expression for linear softening model
//                f[s] = material.sigma_t*(material.wc-width_etoc[s])/material.wc;
//                //expression for exponential model
//                //f[s] = material.sigma_t*6*width_etoc[s]/material.wc*
//                // exp(1-width_etoc[s]*6/material.wc);
//            };
//
//            if (width_etoc[s]>0 and width_etoc[s] < material.wc and
//                width_history[id(i,0)+s] > 0. and
//                width_history[id(i,0)+s] < material.wc) {
//                if (width_etoc[s] < width_history[id(i,0)+s]) {
//                    //expression for linear softening model
//                    f[s] = width_etoc[s] / width_history[id(i,0)+s] *
//                           material.sigma_t*(material.wc-width_history[id(i,0)+s])/material.wc;
//                    //experssion for exponential model
//                    //f[s]=width_etoc[s] / width_history[id(i,0)+s] *
//                    //     material.sigma_t*6*width_history[id(i,0)+s]/material.wc
//                    // *exp(1-6*width_history[id(i,0)+s]/material.wc);
//                }
//                else {
//                    //expression for linear softening model
//                    f[s] = material.sigma_t*(material.wc-width_etoc[s])/material.wc;
//                    //expression for exponential model
//                    //f[s] = material.sigma_t*6*width_etoc[s]/material.wc*
//                    // exp(1-width_etoc[s]*6/material.wc);
//                }
//            };
            // exponential cohesive force over all the elements
            // in this case, we don't take the opening history into account which
            // means that there's no unloading during the propagation process
            if(width_etoc[s]>0  and
               width_history[id(i,0)+s] >= 0.){
                // exponential with one increasing part G=sigmaT*wc*exp(1.0) it's not working with the increasing branch
                //f[s]=material.sigma_t* width_etoc[s]/material.wc *exp(1.0-width_etoc[s]/material.wc);
                // exponential with decreasing part G=sigmaT*wc*0.886227
                f[s]=material.sigma_t* exp(-(width_etoc[s]/material.wc)*(width_etoc[s]/material.wc));
                // exponential with decreasing not square part G=sigmaT*wc
                //f[s]=material.sigma_t*exp(-width_etoc[s]/material.wc);
            };
        };
        return f;
    };//return the local cohesive force, and the force in total




    il::Array<il::int_t> stress_criteria_elmt(const il::Array2D<double> &kmat,
                                             const il::Array2D<int> &id,
                                             const Material &material,
                                             const Initial_condition &initial_condition,
                                             il::Array<double> widthB,
                                             il::Array<double> delta_w,
                                             il::int_t i, il::int_t j,
                                             const il::int_t &dof_dim,
                                             const int &p) {
        //here we need to get the next possible failed elements

        //i, j are elements numbers starting from 0,1,2,3,4...
        il::int_t ttn=(p+1)*dof_dim*id.size(0);
        il::int_t n = j-i+1;
        il::Array2D<double> kt{ttn, (p+1)*dof_dim*n, 0.};

        take_submatrix(kt, 0, int(ttn - 1), id(i,0), id(j,1+dof_dim), kmat);
        il::blas(1.0, delta_w,il::io, widthB);
        il::Array<double> syy = il::dot(kt, widthB);
        il::int_t min_i=i;
        il::int_t max_j=j;

        for (il::int_t r = 0; r < i; r++) {
            if (syy[id(r,1)] >=
                material.sigma_t  and i >= 0) {
                //here personally thinking we should not + initial_condition.sigma0,
                //but not sure
                if(r<=min_i){
                    min_i=r;
                }
            }
        }

        min_i=long(0.85*min_i+0.15*i+0.5);


        for (il::int_t r2 =id.size(0)-1; r2 > j; r2--) {

            if (syy[id(r2,dof_dim+1)] >= material.sigma_t   and
                j <= id.size(0)-1) {
                //here not sure the failure criteria should be + sigma0
                //after some analysis, should not add sigma0
                //the same case for the upper code

                if(max_j<=r2){
                    max_j=r2;
                }
            }
        }
        //we suppose the far point stress is smaller than the near point stress

        max_j=long(0.85*max_j+0.15*j);// before we have +0.5


        if (max_j >= id.size(0)-1) {
            max_j= id.size(0)-1;
        };
        if (min_i < 1) {
            min_i=0;
        };

        il::Array<il::int_t> index{2, 0};
        index[0]=min_i;
        index[1]=max_j;
        return index;
    };





    il::Array<double> conductivities_newtonian_open(const il::Array<double> &rho,
                                                    const il::Array<double> &vector,
                                                    il::Array<double> EltSizes,
                                                    Parameters_fluid &fluid_parameters,
                                                    double kf, Material material, il::io_t) {

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

            if(vector[i]>= material.wc){
                Res[i] = ((rho[i] * (vector[i]*vector[i]*vector[i]))
                          / EltSizes[i]) * (1 / (12 * fluid_parameters.viscosity));
            }
            else{
                // //here we change the permeability of the fracture tip,
                // // one needs to activate this part of the code
                //Res[i]=((rho[i] * kf * (vector[i]*vector[i]*material.wc))
                //        / EltSizes[i]) * (1 / (12 * fluid_parameters.viscosity));
                //if we don't change the permeability around the tip
                Res[i]=((rho[i] * (vector[i]*vector[i]*vector[i]))
                 / EltSizes[i]) * (1 / (12 * fluid_parameters.viscosity));
            }
        }
        return Res;
    }


    il::Array <double> element_size(Mesh mesh_total, const int &p ) {
        //size N*1
        //returns the size of each element in the total mesh
        il::int_t n=mesh_total.nelts();
        il::Array<double> elmt{n,0.};
        for (il::int_t i=0;i<n;i++){
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, int(i), p);
            elmt[i]=segi.size;
        }
        return elmt;
    }


    il::Array<double> average_open(il::Array<double> width,
                                   il::int_t i,
                                   il::int_t j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim){
        il::int_t n=j-i+1;
        il::Array<double> width_average{n, 0.};
            for (il::int_t m = 0; m < n; ++m) {
                width_average[m]=0.5*(width[id(m,1)]+width[id(m,dof_dim+1)]);
            }
        return width_average;
    }



    il::Array<double> quarter_open(il::Array<double> width,
                                   il::int_t i,
                                   il::int_t j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   il::Array2D<int> col_matrix){
        // size is 2*N related element number
        // width the width vector with 4 degrees in each element
        il::int_t n=j-i+1;
        il::Array<double> width_quarter{2*n,0.};

        //the crack covers the whole elements at the tip
            for (il::int_t m = 0; m < n; ++m) {
                width_quarter[col_matrix(m,0)]=
                        0.75*width[id(m,1)]+0.25*width[id(m,dof_dim+1)];
                width_quarter[col_matrix(m,1)]=
                        0.25*width[id(m,1)]+0.75*width[id(m,dof_dim+1)];
            }

        return width_quarter;
    }


    il::Array2D<double> matrix_lw(il::Array<double> widthB,
                                  il::int_t i,
                                  il::int_t j,
                                  const il::int_t &dof_dim, const int &p,
                                  il::Array<double> element_size_all,
                                  Parameters_fluid &fluid_parameters,
                                  const il::Array2D<int> &id, Material material) {

        //widthB should be the adjusted-width profile at the previous moment
        // the size is (N+1)*(N+1)
        // N here means the number of the related elements
        il::int_t n = j - i + 1;
        //number of related elements

        il::Array2D<double> LL{n + 1, n + 1, 0.};
        //contained node number
        il::Array<double> element_size{j - i + 1, 0.};
        //size of related element number
        il::Array<double> w_mid=
                average_open(widthB,i,j,id,dof_dim);
        //size of related element number
        for (il::int_t m = 0; m < j - i + 1; m++) {
            element_size[m] = element_size_all[m + i];
        };
        il::Array2D<double> density{j - i + 1,2,0.};
        take_submatrix(density,int(i), int(j),0,1,fluid_parameters.density);
        //this should be size of related element number

        il::Array<double> rho_mid = average(density,il::io);
        il::Array<double> Kk =
                conductivities_newtonian_open(rho_mid, w_mid, element_size,
                                              fluid_parameters, 0., material, il::io);//this is where we should change the permeability coefficient kf
        //size of related element number

        il::int_t dofj = 0;
        il::Array2D<int> Dofp = dofhandle_cg2d(2, j - i + 1, il::io);

        // Loop over the pressure nodes
        for (il::int_t s = 0; s < n + 1; ++s) {

            il::Array2D<int> ed = search(Dofp, int(s), il::io);
            //related element and the position of the pressure node
            // possible rows and columns

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {
                //number of related elements
                il::int_t ej = ed(sc, 0);//row number-element number
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {//pressure node is the left node of the element
                    //dofj corresponds to the other point of the element
                    dofj = Dofp(ej, 1);
                } else {//pressure node is the right node of the element
                    dofj = Dofp(ej, 0);
                }
                LL(s, s) = LL(s, s) - Kk[ej];
                if(dofj>=0 && dofj<n+1){
                    LL(s, dofj) = LL(s, dofj) + Kk[ej];
                }
            }
        }

        return LL;
    }//should times the time step to build the whole matrix


    il::Array2D<double> matrix_vpw(Parameters_fluid &fluid_parameters,
                                   il::Array<double> element_size_all,
                                   il::int_t i,
                                   il::int_t j,
                                   il::Array<double> widthB,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   const int &p,
                                   const il::Array2D<int> &col_matrix,
                                   il::io_t) {
        //(N+1)*(N+1)
        //N is the number of the related elements

        // create the array of element size
        il::int_t n = j - i + 1;
        il::Array<double> element_size{n, 0.};
        for (il::int_t m = 0; m < n; m++) {
            element_size[m] = element_size_all[m + i];
        };

        // Create all the vectors for compressibility of fluid
        // this should be changed with the actual activated collocation points
        // Vector of compressibility of the fluid at nodal points
        il::Array<double> Cf{n+1, fluid_parameters.compressibility};
        // Vector of compressibility of the fluid at the midpoints of each element
        il::Array<double> Cfmid{n, fluid_parameters.compressibility};
        // Vector of compressibility of the fluid at +/- 1/4 of each element
        il::Array<double> Cfquart{2*n, fluid_parameters.compressibility};



        il::Array2D<double> Vp{n+1, n+1, 0};
        il::int_t dofj = 0;
        il::Array2D<int> Dofp = dofhandle_cg2d(2, j - i + 1, il::io);

        il::Array<double> w_mid= average_open(widthB,i,j,id,dof_dim);
        il::Array<double> wquart= quarter_open(widthB,i,j, id, dof_dim,col_matrix);

        /// Assembling procedure ///
        //Loop over the pressure node
        for (il::int_t m = 0; m < n+1; ++m) {

            il::Array2D<int> ed = search(Dofp, int(m), il::io);
            //find the related element and the position of the pressure node

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {
                /// sc-one certain related element number

                il::int_t ej = ed(sc, 0);//row number-element
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {
                    // pressure node situates at the left node of the element
                    //dofj corresponds to the other point of the element
                    dofj = Dofp(ej, 1);
                    Vp(m, m) = Vp(m, m) +(element_size[ej] / 12) *
                                         ((widthB[id(ej,1)] * Cf[m])
                                          + (0.5 * w_mid[ej] * Cfmid[ej])
                                          + (3 * wquart[col_matrix(ej,0)] *
                                             Cfquart[col_matrix(ej,0)]));
                    if(dofj>=0 && dofj<n+1){
                        Vp(m, dofj) = Vp(m, dofj) + (element_size[ej] / 12) *
                                                    ((0.5 * w_mid[ej] * Cfmid[ej])
                                                     + (wquart[col_matrix(ej,0)] *
                                                        Cfquart[col_matrix(ej,0)]));
                    }

                }
                else {// pressure node situates at the right node of the element
                    //dofj corresponds to the other point of the element
                    dofj = Dofp(ej, 0);
                    Vp(m, m) = Vp(m, m) +(element_size[ej] / 12) *
                                         ((widthB[id(ej,1+dof_dim)] * Cf[m])
                                          + (0.5 * w_mid[ej] * Cfmid[ej])
                                          + (3 * wquart[col_matrix(ej,1)] *
                                             Cfquart[col_matrix(ej,1)]));
                    if(dofj>=0 && dofj<n+1){
                        Vp(m, dofj) = Vp(m, dofj) + (element_size[ej] / 12) *
                                                    ((0.5 * w_mid[ej] * Cfmid[ej])
                                                     + (wquart[col_matrix(ej,1)] *
                                                        Cfquart[col_matrix(ej,1)]));
                    }

                }
            }
        }

        return Vp;
    };

    //build Vw matrix
    il::Array2D<double> matrix_vw(Parameters_fluid &fluid_parameters,
                                  il::Array<double> element_size_all,
                                  il::int_t i,
                                  il::int_t j,
                                  il::Array<double> widthB,
                                  const il::Array2D<int> &id,
                                  const il::int_t &dof_dim, const int &p,
                                  const il::Array2D<int> &col_matrix,il::io_t) {


        il::int_t n = j - i + 1;
        il::Array<double> element_size{n, 0.};
        for (il::int_t m = 0; m < n; m++) {
            element_size[m] = element_size_all[m + i];
        };
        il::Array2D<double> density{j-i+1,2,0.};
        take_submatrix(density,int(i),int(j),0,1,fluid_parameters.density);
        //this should be size of all related element number
        il::Array<double> rho_mid=average(density, il::io);;
        il::Array<double> rho_quart=quarter(fluid_parameters.density, il::io);


        il::Array2D<double> Vw{n+1, widthB.size(), 0};

        il::Array2D<int> Dofp = dofhandle_cg2d(2, j - i + 1, il::io);

        /// Assembling procedure ///
        for (il::int_t m = 0; m < n+1; ++m) {

            il::Array2D<int> ed = search(Dofp, int(m), il::io);

            for (il::int_t sc = 0; sc < ed.size(0); ++sc) {
                //sc-one certain related element number

                il::int_t ej = ed(sc, 0);//row number-element
                il::int_t ec = ed(sc, 1);//coloumn number-left or right node
                if (ec == 0) {//left point
                    Vw(m, id(ej,1)) =
                            Vw(m, id(ej,1))
                            +(element_size[ej] / 12) *(density(ej,ec) + 0.5 * rho_mid[ej] +
                                                       3 * rho_quart[col_matrix(ej,0)]);
                    if(id(ej,1)+dof_dim<widthB.size()){
                        Vw(m, id(ej,1)+dof_dim) =
                                Vw(m, id(ej,1) +dof_dim)
                                + (element_size[ej] / 12) * (0.5 * rho_mid[ej]
                                                             + rho_quart[col_matrix(ej,0)]);
                    }

                }
                else {//right point
                    Vw(m, id(ej,1+dof_dim)) =
                            Vw(m, id(ej,1+dof_dim))
                            +(element_size[ej] / 12) *
                             ( density(ej,ec) + 0.5  * rho_mid[ej] +
                               3 * rho_quart[col_matrix(ej,1)]);
                    if(id(ej,1+dof_dim)-dof_dim>=0){
                        Vw(m, id(ej,1+dof_dim)-dof_dim) =
                                Vw(m, id(ej,1+dof_dim) -dof_dim)
                                + (element_size[ej] / 12) * (0.5  * rho_mid[ej]
                                                             + rho_quart[col_matrix(ej,1)]);
                    }
                }
            }
        }
        return Vw;
    };





    il::int_t find_source(il::int_t i0, il::int_t j0) {
        //we consider i0 and j0 represent the element numbers and i0
        // and j0 begin from 0.1.2.3.4...
        //and there's a requirement for i0 and j0, j0-i0=2, or j0-i0=1. j0-i0=2
        // would be better for the paired Mesh
        //returns the node number of injection point

        il::int_t input;

        if (j0 - i0 == 1) {
            input= j0;
        } else {
            input = (j0 - i0+1) / 2+i0;
        }
        return input;
    }



// the matrix for now is built in a not-fully-implicit way,
// but we can build it in a fully-implicit way

    void construct_matrix_visco(il::Array2D<double> &kmatC,
                                il::Array2D<double> vwc, il::Array2D<double> vp,
                                il::Array2D<double> ll,
                                il::Array<double> cohf, const Material &material,
                                il::Array2D<double> matrix_edge_to_col,
                                il::Array<double> m_source,
                                const Initial_condition &initial_condition,
                                il::Status &status,
                                il::Array<double> pressure_f,
                                il::Array<double> widthB,
                                const int &p, const il::int_t &dof_dim,
                                double time_inter,il::io_t,
                                il::Array<double> &width,
                                il::Array<double> &pressure,
                                il::Array<double> &volume_vary,
                                il::Array<double> & elastic_vary) {
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
        // we embed the effect of the confining stress into fP
        il::Array<double> lp_b=il::dot(ll,pressure_f);
        for (int c = 0; dof_dim * c < n; c++) {
            fP[2 * c] = 0;
        }

        for (int m = 0; m < n; ++m) {
            newvector[m]=cohf[m] - fP[m] - fK[m];
            for (int s = 0; s < n; ++s){
                newmatrix(m, s) = kmatC(m, s);
            }
            for (int s1=0; s1<n2;++s1){
                newmatrix(s1+n,m)=vwc(s1,m);
            }
        };
        for (int m1=0;m1<n2;m1++){
            newvector[n+m1]= (lp_b[m1]+m_source[m1])*time_inter;
            for (int s2=0;s2<n;++s2){
                newmatrix(s2,m1+n)=matrix_edge_to_col(s2,m1);
            }
            for (int s3=0;s3<n2;s3++){
                newmatrix(s3+n,m1+n)=
                        vp(s3,m1)-ll(s3,m1)*time_inter;
            }
        };


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
                               il::int_t i, il::int_t j,
                               const il::Array2D<double> &kmat,
                               const il::Array2D<int> &id,const int &p,
                               il::Array<double> widthB,
                               il::Array<double> pressure_f,
                               il::Array<double> width_history,
                               const il::int_t &dof_dim,
                               const il::Array2D<int>&col_matrix,
                               Parameters_fluid fluid_parameters,
                               il::Array<double> element_size_all,
                               const il::Array2D<double> &matrix_edge_to_col_all,
                               il::int_t input,double time_current,
                               il::Array2D<double> width_e_to_c,
                               il::io_t,
                               il::Array<double> &delta_width,
                               il::Array <double> &pressure_change,
                               il::Array<double> &coht,int &mm,
                               il::Array<double> &volume_vary,
                               il::Array<double> &elastic_vary,
                               il::Array<double> &error_w_list) {

        il::Array2D<double> kmatnew{dof_dim *(p+1)* (j - i + 1),
                                    dof_dim *(p+1)*(j - i + 1), 0.};
        take_submatrix(kmatnew, id(i,0), id(j,dof_dim+dof_dim-1),
                       id(i,0), id(j,dof_dim+dof_dim-1), kmat);

        il::int_t n = j - i + 1;

        il::Array2D<double> vwc{n+1,dof_dim *(p+1)*(j - i + 1), 0.};
        vwc=matrix_vw(fluid_parameters, element_size_all, i,j,
                      widthB, id, dof_dim, p,col_matrix,il::io);
        il::Array2D<double> ll{n+1,n+1,0.};
        il::Array2D<double> vp{n+1,n+1,0.};
        il::Array2D<double> matrix_edge_to_col{id(j, dof_dim+dof_dim-1)
                                               -id(i, 0)+1, j+1-i+1,0.};
        // size 4N(all related element number)*(N+1)(only concerned pressure node)
        take_submatrix(matrix_edge_to_col,id(i,0), id(j,dof_dim+dof_dim-1),
                       int(i),int(j+1),matrix_edge_to_col_all);




        il::Array<double> m_source{n+1,0.};
        m_source[input-i]=
                initial_condition.Q0*fluid_parameters.density(input,0);

        il::Array<double> delta_w_ini{widthB.size(), 0.};

        il::Array<double> coh{widthB.size(), 0.};

        int k = 0;
        double error_w = 1.;
        il::Array<il::Status> statusk{200,};
        double error_p=1.;
        //double error_r=1;
        il::Array<double> error_list{200,0.};

        il::Array<double> delta_w_bitera=delta_w_ini;
        il::Array<double> delta_p_ini{n+1,0.};
        il::Array<double> delta_p_itera=delta_p_ini;

        il::Array<double> volume_vary_inter{n+1,0.};
        il::Array<double> elastic_vary_inter{widthB.size(),0.};


        //double alpha=0.1;

        while (k < 200 && (error_w > 0.00001 or error_p>0.00001)) {
            il::Array<double> width_inm = widthB;
            il::blas(1.0, delta_w_bitera, 1.0, il::io, width_inm);

            //Linear Softening model
            coh = cohesive_force_linear_elemt(material, width_inm, i,
                                        j, width_history, p,dof_dim,id,width_e_to_c);

            ll=matrix_lw(width_inm, i,j,dof_dim, p,
                         element_size_all, fluid_parameters, id, material);
            vp=matrix_vpw(fluid_parameters,element_size_all,i,j,
                          width_inm, id, dof_dim, p,col_matrix,il::io);


            construct_matrix_visco(kmatnew, vwc, vp,ll,coh, material,
                                   matrix_edge_to_col,m_source,
                                   initial_condition, statusk[k],
                                   pressure_f, widthB, p, dof_dim,
                                   time_current,
                                   il::io, delta_width, pressure_change,
                                   volume_vary_inter,elastic_vary_inter);


            //to force the increment of w who contributes a negative
            // opening equal to zero
            il::Array<double> width_check=widthB;
            il::blas(1.0, delta_width, 1.0, il::io, width_check);
            //delta_w_ini could also be replaced by delta_width
            for(int w_check=0;w_check<widthB.size();w_check++){
                if(width_check[w_check]<0.){
                    delta_width[w_check]=-1.*widthB[w_check];
                    //to force the total opening is zero
                }
            }
            //relaxation of the increment of the dislocation
            il::blas(0.15, delta_w_bitera, 0.85, il::io, delta_width);

            //calculation of the relative errors
            il::Array<double> intermediate = delta_width;
            il::blas(-1.0, delta_w_bitera, 1.0, il::io, intermediate);
            double error_inter = il::norm(intermediate, il::Norm::L2);
            error_w = error_inter / il::norm(delta_width, il::Norm::L2);

            //iteration condition
            delta_w_bitera=delta_width;


            il::Array<double> inter_pressure=pressure_change;
            il::blas(-1.0,delta_p_ini,1.0,il::io,inter_pressure);
            double error_inter_p=il::norm(inter_pressure,il::Norm::L2);
            error_p=error_inter_p/il::norm(pressure_change,il::Norm::L2);
            //error_p=0.0;
            //this line tries to make the pressure iteration criteria unavalible
            //if activated

            //update of the relative error for increment of the pressure
            delta_p_ini=pressure_change;

            //check the iteration scheme is suitable or not
            error_list[k]=error_w;//error_w;error_p;
            ++k;

        };
        mm=k;

        //plot the current coh vector
        for (il::int_t r = 0; r < widthB.size(); r++) {
            coht[id(i,0)+ r] = coh[r];
        }

        volume_vary =volume_vary_inter;
        elastic_vary=elastic_vary_inter;
        error_w_list=error_list;

    }


    void
    propagation_loop_visco(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                           const Material &material,
                           const Initial_condition &initial_condition,
                           il::int_t i0, il::int_t j0, int nstep,
                           il::Status &status, Parameters_fluid fluid_parameters,
                           il::io_t,
                           il::Array2C<double> &widthlist,
                           il::Array2D<double> &plist, il::Array<double> &l_coh,
                           il::Array<double> &l, il::Array2C<double> &coh_list,
                           il::Array<int> &mvalue,int &break_time,
                           il::Array2C<double> &stress_list,
                           il::Array2D<double> &volume_vary_list,
                           il::Array2D<double> &elastic_vary_list,
                           il::Array2D<double> &error_matrix,
                           il::Array<double> &timelist) {

        //i0 and j0 are the beginning element number in C++, starting from 0,1,2.....

        int n = mesh_total.nelts();
        IL_EXPECT_FAST(n == id.size(0));
        il::int_t dof_dim=id.size(1)/(p+1);

        //element size for all
        il::Array<double> element_size_all=element_size(mesh_total,p);
        //preparation for w_quarter numbering
        il::Array2D<int> col_matrix=collocation_matrix(n,p);

        //bluid the elastic matrix K
        il::int_t ttn=id.size(0)*id.size(1);
        il::Array2D<double> kmat{ttn,ttn,0.};
        kmat=basic_assembly(mesh_total, id, p, material.Ep);

        //Preparation for bluiding the mass-conservation matrix
        il::Array2D<int> Dofp = dofhandle_cg2d(int(dof_dim),n,il::io);

        //bluid the Edge to Collocation point matrix
        il::Array2D<double> matrix_edge_to_col_all =
                from_edge_to_col_cg(int(dof_dim), id,Dofp, il::io);

        //maximum of the opening
        il::Array<double> width_history{ttn, 0.};
        //give one to the initial opening to avoid the cohesive force
        for (il::int_t h = i0; h < j0 + 1; ++h) {
            width_history[id(h,dof_dim*0+1)] = 1.;
            width_history[id(h,dof_dim*1+1)] = 1.;
        };

        //matrix needed to calculate the dislocation at the collocation point
        il::Array2D<double> width_e_to_c=width_ec(dof_dim,id);


        //inital crack opening
        il::Array<double> widthB =
                initialwidth_col(kmat, id, p, initial_condition, i0, j0,
                                 status,dof_dim);
        //vector of opening for the failed elements
        il::Array<double> width;
        il::Array2C<double> widthlistinter{nstep, ttn,0.};
        //stress plot
        il::Array2C<double> stress_profile{nstep,ttn,0.};
        //vector of time's increment for each time step
        il::Array<double> time_inter{nstep,0.};
        time_inter[0]=initial_condition.timestep;
        //vector of the acculated time and initialization
        il::Array<double> timelist_inter{nstep+1,0.};

        //array of pressure along the whole mesh for each time step
        il::Array2D<double> plistinter{nstep + 1, n+1,0.};
        il::Array<double> pressure;
        il::Array<double> pressure_change;
        //initialization of the pressure vector of the previous time step
        il::Array<double> pressure_b{j0-i0+1+1,initial_condition.pini};

        //initialization of the iteration times during plasticity
        //for each time step
        il::Array<int> m_value{nstep,0};
        int mk=0;

        //Initialization of the current number of timestep
        int s = 0;
        //Control of iteration times of number of plasticity and stress-criteria
        int m = 0;
        //the initial left element number of the crack
        il::int_t i = i0;
        //the initial right element number of the crack
        il::int_t j = j0;
        //the location of the source(injection point)
        il::int_t input=find_source(i0,j0);

        il::int_t left_node=i0;

        //Initialization of the crack length
        il::Array<double> crack_length{nstep+1, 0.};
        double l0=0.;
        for(il::int_t ei=i0;ei<j0+1;ei++){
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, int(ei), p);
            l0+=segi.size;
        }
        crack_length[0]=l0;

        //Initializaiton of the coehsive length
        il::Array<double> length_coh{nstep, 0.};

        //Initialization of cohesive forces along all the mesh
        //for each time step
        il::Array2C<double> coh_tt{nstep, ttn, 0.};
        il::Array<double> coht{ttn, 0.};

        //to output the residual of the volume control equation
        il::Array2D<double> volume_vary_matrix{nstep,n+1,0.};
        //output the residul of the elastic equations
        il::Array2D<double> elastic_vary_matrix{nstep,ttn,0.};
        // evolution of error_w of each iteration for every time step
        il::Array2D<double> error_lists{nstep,200,0.};


        //Control of the total time taken for the crack to reach the mesh end
        break_time=nstep;
        bool break_value=false;
        bool break_value_2=false;


        while (s < nstep) {

            //to check the iteration scheme is suitable or not
            il::Array<double> error_w_list;

            while (m < 10000) {
                //the stop limit value for m should be big enough
                if (break_value) {
                    break_value_2=true;
                    // to stop the calculation for the next time step
                    break;
                };
                il::Array<double> delta_width;
                il::Array<double> volume_vary;
                il::Array<double> elastic_vary;

                plasticity_loop_visco(material, initial_condition,i, j,
                                      kmat, id, p, widthB, pressure_b,
                                      width_history,dof_dim,col_matrix,
                                      fluid_parameters, element_size_all,
                                      matrix_edge_to_col_all,input,time_inter[s],
                                      width_e_to_c,il::io,
                                      delta_width, pressure_change, coht,mk,
                                      volume_vary,elastic_vary,error_w_list);

                m_value[s]=mk;

                il::Array<il::int_t> index =
                        stress_criteria_elmt(kmat, id, material, initial_condition,
                                            widthB, delta_width, i, j,
                                            dof_dim,p);

                //the crack doesn't propagate, goes to the next time step

                if (i == index[0] and j == index[1]) {
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    pressure=pressure_b;
                    il::blas(1.0,pressure_change,1.0,il::io,pressure);
                    timelist_inter[s+1]=timelist_inter[s]+time_inter[s];
                    break;
                }
                //the crack reaches the end of the mesh
                if (index[0] <= 0 or index[1] >= n- 1) {

                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    pressure=pressure_b;
                    il::blas(1.0,pressure_change,1.0,il::io,pressure);
                    timelist_inter[s+1]=timelist_inter[s]+time_inter[s];

                    break_time=s+1;//because s starts from 0
                    break_value=true;
                    // to stop the plasticity and the front position loop
                    break;
                }
                //the crack propagates
                //rearrange the width vector
                il::Array<double> widthBnew{
                        dof_dim *(p+1)*(index[1] - index[0] + 1), 0.};
                for (int wbn = 0; wbn < widthB.size(); wbn += 1) {
                    widthBnew[dof_dim *(p+1)*(i - index[0]) + wbn] = widthB[wbn];
                }

                //rearrange the vector of width increment and elastic residual
                il::Array<double> delta_width_new{
                        dof_dim *(p+1)* (index[1] - index[0] + 1), 0.};
                for(int ela_i=0;ela_i<elastic_vary_matrix.size(1);ela_i++){
                    elastic_vary_matrix(s,ela_i)=0.;
                }
                for (int dwn = 0; dwn < widthB.size(); dwn += 1) {
                    delta_width_new[dof_dim * (p+1)*(i - index[0]) + dwn]
                            = delta_width[dwn];
                    elastic_vary_matrix(s,dof_dim*(p+1)*i+dwn)=elastic_vary[dwn];
                }

                //Updating the vector of width of previous time step
                widthB.resize(dof_dim *(p+1)* (index[1] - index[0] + 1));
                widthB = widthBnew;
                width = widthBnew;

                //current left tip, the node starting from 0,1,2,3...
                left_node=index[0];

                //the current number of pressure nodes
                il::int_t pn=index[1]-index[0]+1+1;


                //left tip of previous time step, the node starting from 0,1,2...
                il::int_t left_node_b=i;
                //number of pressure nodes of previous time step
                il::int_t pn_b=j-i+1+1;

                //Updating the pressure vector and the volume residual
                il::Array<double> pressure_b_new{pn,0.};
                il::Array<double> pressure_change_new{pn,0.};

                for(int v_clear=0;v_clear<n+1;v_clear++){
                    volume_vary_matrix(s,v_clear)=0.;
                };

                for(int pb_i=0;pb_i<pn_b;pb_i++){
                    pressure_b_new[pb_i-left_node+left_node_b]=pressure_b[pb_i];
                    pressure_change_new[pb_i-left_node+left_node_b]=
                            pressure_change[pb_i];
                    volume_vary_matrix(s,pb_i+left_node_b)=volume_vary[pb_i];
                }
                pressure_b.resize(pn);
                pressure_b = pressure_b_new;
                pressure = pressure_b_new;

                //Updating the tip position or the failed element number
                i = index[0];
                j = index[1];

                //Updating the width and the pressure at the end of the time step
                il::blas(1., delta_width_new, 1.0, il::io, width);
                il::blas(1.,pressure_change_new,1.0,il::io,pressure);
                ++m;
            }
            if(break_value_2){
                break;
            };



            //Arrangement for the width plot over the whole mesh
            il::Array<double> width_large{dof_dim*(p+1)*n, 0.};
            for (int r = 0; r < width.size(); ++r) {
                width_large[dof_dim *(p+1)*i + r] = width[r];
            }//needs to be carefully checked

            //Arrangement for the width plot for all time steps
            for (int sc = 0; sc < dof_dim*(p+1)* n; ++sc) {
                widthlistinter(s, sc) = width_large[sc];
            }
            il::Array<double> pressure_large{n+1,0.};

            // Arrangement for the pressure over the whole Mesh
            for(int sp=0;sp<pressure.size();sp++){
                plistinter(s,left_node+sp)=pressure[sp];
                pressure_large[left_node+sp]=pressure[sp];
            }


            // Calculation of the crack length and the cohesive length
            cc_length_elemt(mesh_total, material, width_large, i, j,
                          width_history, p,id,col_matrix,dof_dim, il::io,
                          length_coh[s], crack_length[s+1]);

            // Inner control of the elements failed(propagation) per time step
            // this is based on the constant propagation velocity assumption
            double psi=40;//before there's no *20
            if(crack_length[s+1]==crack_length[s]){
                time_inter[s+1]=2.*time_inter[s];
            }
            else{
                time_inter[s+1]= psi*element_size_all[i]
                                 /(crack_length[s+1]-crack_length[s])*time_inter[s];
            }


            // Recording the maximum opening
            width_history =
                    write_history_elemt(width, width_history,
                                        i,j, material,p,dof_dim,id,width_e_to_c);
            //Here we use width_large, could be simpler for the coding

            // Arrangement for the cohesive force vector for different time steps
            for (int cohi = 0; cohi < dof_dim*(p+1)* n; cohi++) {
                coh_tt(s, cohi) = coht[cohi];
            }


            // Calculation of the current stress plot caused by the dislocation
            il::Array<double> pressure_current=
                    il::dot(matrix_edge_to_col_all,pressure_large);
            il::Array<double> stress_current=il::dot(kmat,width_large);
            for (int stre=0;stre<dof_dim*(p+1)*n;stre++){
                stress_profile(s,stre)=
                        stress_current[stre]+pressure_current[stre]
                        -initial_condition.sigma0;
            }

            // Error list arrangement for all timesteps, 200 is the
            // maximum iteration number in plastisity loop
            for(int iter=0;iter<200;iter++){
                error_lists(s,iter)=error_w_list[iter];
            }

            widthB = width;
            pressure_b=pressure;
            //time += time_inter[s];
            timelist_inter[s+1]=timelist_inter[s]+time_inter[s];
            ++s;

        }
        widthlist = widthlistinter;
        plist = plistinter;
        l = crack_length;
        l_coh = length_coh;
        coh_list = coh_tt;
        timelist=timelist_inter;
        mvalue=m_value;
        //to output the iteration times in plasticity_loop to get the solution
        stress_list=stress_profile;
        volume_vary_list=volume_vary_matrix;
        elastic_vary_list=elastic_vary_matrix;
        error_matrix=error_lists;
        status.abort_on_error();
    }

}