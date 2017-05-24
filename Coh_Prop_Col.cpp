//
// Created by DONG LIU on 2/23/17.
//

#include "Coh_Prop_Col.h"
#include<il/norm.h>



namespace hfp2d {

    void
    material_condition_col(const double& wc1, const double& sigma_t1,
                           const double& Ep1, const double& U1,
                           const double& pini1, const double& Q01,
                           const double& timestep1, const double& sigma01,
                           il::io_t,Material& material,
                           Initial_condition& initial_condition) {
        material.Ep = Ep1;
        material.sigma_t = sigma_t1;
        material.U = U1;
        material.wc = wc1;
        initial_condition.pini = pini1;
        initial_condition.Q0 = Q01;
        initial_condition.sigma0 = sigma01;
        initial_condition.timestep = timestep1;
    }

    il::Array<double> get_vwc_vc_col(Mesh mesh_total,const int &p,
                                     const il::int_t &dof_dim) {

        int n = mesh_total.nelts();
        il::int_t dof = dof_dim * (p + 1);
        il::Array<double> vwcinter{dof * n, 0.};
        for (int i = 0; i < n; ++i) {
            SegmentCharacteristic segi = get_segment_DD_characteristic(mesh_total,i, p);
            vwcinter[dof * i + 1] = 0.5 * segi.size;
            vwcinter[dof * i + dof_dim+1] = 0.5 * segi.size;
        };
        return vwcinter;
    }


    il::Array<il::int_t> stress_criteria_col(const il::Array2D<double> &kmat,
                                             const il::Array2D<int> &id,
                                             const Material &material,
                                             const Initial_condition &initial_condition,
                                             il::Array<double> widthB,
                                             il::Array<double> delta_w,
                                             il::int_t i, il::int_t j,
                                             il::int_t c_i, il::int_t c_j,
                                             const il::int_t &dof_dim,
                                             const il::Array2D<int> &col_matrix,
                                             const int &p) {
        //here we need to get the collocation points of the next-possible
        // failed elements not those of the activated elements

        il::int_t ttn=(p+1)*dof_dim*id.size(0);
        il::int_t col_n = c_j-c_i+1;
        il::Array2D<double> kt{ttn, dof_dim*col_n, 0.};
        il::Array2D<int> col_row_i=search(col_matrix,int(c_i),il::io);
        il::Array2D<int> col_row_j=search(col_matrix,int(c_j),il::io);


        take_submatrix(kt, 0, int(ttn - 1), id(col_row_i(0,0),dof_dim*col_row_i(0,1)),
                       id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1), kmat);
        il::blas(1.0, delta_w,il::io, widthB);
        il::Array<double> syy = il::dot(kt, widthB);
        il::int_t min_c_i=c_i;
        il::int_t max_c_j=c_j;

        for (il::int_t r = 0; r < c_i; r++) {
            il::Array2D<int> col_row_r=search(col_matrix,int(r),il::io);
            if (syy[id(col_row_r(0,0),dof_dim*col_row_r(0,1)+1)] >=
                material.sigma_t  and c_i >= 0) {
                //here the failure criteria should add the sigma0 to become the total insitu-stress,
                // not sure. + initial_condition.sigma0
                //
                if(r<=min_c_i){
                    min_c_i=r;
                }
            }
        }

        c_i=long(0.85*min_c_i+0.15*c_i+0.5);
        il::Array2D<int> col_row_min_ci=search(col_matrix,int(c_i),il::io);
        i=col_row_min_ci(0,0);


        for (il::int_t r2 =ttn/(p+1)-1; r2 > c_j; r2--) {
            il::Array2D<int> col_row_r2=search(col_matrix,int(r2),il::io);
            if (syy[id(col_row_r2(0,0),dof_dim*col_row_r2(0,1)+1)] >=
                material.sigma_t   and
                c_j <= ttn/(p+1)-1) {
                //here not sure the failure criteria should be + sigma0
                //after some analysis, should not add sigma0
                //the same case for the upper code

                if(max_c_j<=r2){
                    max_c_j=r2;
                }
            }
        }
        //with this, we can accelerate the propagation,
        // but we suppose the far point stress is smaller than the near point stress

        c_j=long(0.85*max_c_j+0.15*c_j+0.5);
        il::Array2D<int> col_row_max_cj=search(col_matrix,int(c_j),il::io);
        j=col_row_max_cj(0,0);


        if (c_j >= ttn/(p+1)-1) {
            c_j =ttn/(p+1)-1;
            j= id.size(0)-1;
        };
        if (c_i < 1) {
            c_i = 0;
            i=0;
        };
        il::Array<il::int_t> index{4, 0};
        index[0] = c_i;
        index[1] = c_j;
        index[2]=i;
        index[3]=j;
        return index;

    };



    il::Array<double> width_edge_col_col(il::Array<double> width,
                                         il::Array2D<int> col_row_i,
                                         il::Array2D<int> col_row_j,
                                         const il::int_t &dof_dim,
                                         const il::Array2D<int> &id) {
        il::Array<double> width_new{width.size(), 0.};
        if(col_row_i(0,1)==0 and col_row_j(0,1)==1)
        {
            for (il::int_t i = 0; i < col_row_j(0,0)-col_row_i(0,0)+1; ++i) {
                width_new[id(i,1)] = (0.5 + sqrt(2) / 4.) * width[id(i,1)] +
                                         (0.5 - sqrt(2) / 4.) * width[id(i,dof_dim+1)];
                width_new[id(i,dof_dim+1)] = (0.5 - sqrt(2) / 4.) * width[id(i,1)] +
                                         (0.5 + sqrt(2) / 4.) * width[id(i,dof_dim+1)];
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==1){
            width_new[1]=(0.5 + sqrt(2) / 4.) * width[1];
            for (il::int_t i = 1; i < col_row_j(0,0)-col_row_i(0,0)+1; ++i) {
                width_new[id(i,1)-dof_dim] = (0.5 + sqrt(2) / 4.) * width[id(i,1)-dof_dim] +
                                     (0.5 - sqrt(2) / 4.) * width[id(i,dof_dim+1)-dof_dim];
                width_new[id(i,dof_dim+1)-dof_dim] = (0.5 - sqrt(2) / 4.) * width[id(i,1)-dof_dim] +
                                                 (0.5 + sqrt(2) / 4.) * width[id(i,dof_dim+1)-dof_dim];
            }
        }

        if(col_row_i(0,1)==1 and col_row_j(0,1)==0){
            width_new[1]=(0.5 + sqrt(2) / 4.) * width[1];
            for (il::int_t i = 1; i < col_row_j(0,0)-col_row_i(0,0); ++i) {
                width_new[id(i,1)-dof_dim] = (0.5 + sqrt(2) / 4.) * width[id(i,1)-dof_dim] +
                                     (0.5 - sqrt(2) / 4.) * width[id(i,dof_dim+1)-dof_dim];
                width_new[id(i,dof_dim+1)-dof_dim] = (0.5 - sqrt(2) / 4.) * width[id(i,1)-dof_dim] +
                                             (0.5 + sqrt(2) / 4.) * width[id(i,dof_dim+1)-dof_dim];
            }

            width_new[id(col_row_j(0,0),dof_dim*col_row_j(0,1)+1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))]=
                    (0.5 + sqrt(2) / 4.) * width[id(col_row_j(0,0),dof_dim*col_row_j(0,1)+1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))];
        }

        if(col_row_i(0,1)==0 and col_row_j(0,1)==0){
            width_new[id(col_row_j(0,0),dof_dim*col_row_j(0,1)+1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))]=
                    (0.5 + sqrt(2) / 4.) * width[id(col_row_j(0,0),dof_dim*col_row_j(0,1)+1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))];
            for (il::int_t i = 0; i < col_row_j(0,0)-col_row_i(0,0); ++i) {
                width_new[id(i,1)] = (0.5 + sqrt(2) / 4.) * width[id(i,1)] +
                                     (0.5 - sqrt(2) / 4.) * width[id(i,dof_dim+1)];
                width_new[id(i,dof_dim+1)] = (0.5 - sqrt(2) / 4.) * width[id(i,1)] +
                                             (0.5 + sqrt(2) / 4.) * width[id(i,dof_dim+1)];
            }
        }
        return width_new;
    }


    il::Array<double>
    write_history_col(il::Array<double> width, il::Array<double> width_history,
                      il::int_t c_i, il::int_t c_j,il::Array2D<int> col_matrix,
                      const Material &material, const int &p,
                      const il::int_t &dof_dim,const il::Array2D<int> &id) {

        il::int_t n = width.size();
        il::int_t n_t = width_history.size();
        il::Array<double> width_history_new{n_t, 0.};
        il::Array<double> width_new{n, 0.};
        il::Array2D<int> col_row_i=search(col_matrix,int(c_i),il::io);
        il::Array2D<int> col_row_j=search(col_matrix,int(c_j),il::io);
        width_new = width_edge_col_col(width,col_row_i, col_row_j,dof_dim,id);
        //changes the 7th April, m initial value from 0 to 1 a
        // nd the varying step from 1 to 2,
        // which makes the width_history only relates to the opening not the slip.
        for (il::int_t m = 0; m < n; m += 1) {
            width_history_new[id(col_row_i(0,0),dof_dim*col_row_i(0,1)) + m] =
                    width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+ m];
            if (width_new[m] > width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1)) + m]) {

                width_history_new[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+ m] = width_new[m];
            }
        }
        return width_history_new;
    }


    il::Array<double> cohesive_force_col(const Material &material,
                                         il::Array<double> width,
                                         il::Array2D<int> col_row_i,
                                         il::Array2D<int> col_row_j,
                                         il::Array<double> width_history,
                                         const int &p, const il::int_t &dof_dim,
                                         const il::Array2D<int> &id) {
        il::Array<double> f{width.size(), 0.};//material.sigma_t*0.00001 one way to avoid the jump, seems not working

        for (il::int_t s = 1; s < width.size(); s += dof_dim) {


            il::Array<double> width_etoc =
                    width_edge_col_col(width,col_row_i,col_row_j,dof_dim,id);
            //width_etoc[s]>0 and not decided yet
            if (width_etoc[s]>0 and width_etoc[s] < material.wc and
                width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+s] == 0.) {
                f[s] = material.sigma_t;
            };

            if (width_etoc[s]>0 and width_etoc[s] < material.wc and
                width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+s] > 0. and
                width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+s] <
                material.wc) {
                if (width_etoc[s] < width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+s]) {
                    f[s] = width_etoc[s] / width_history[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+s] *
                           material.sigma_t;
                }
                else {
                    f[s] = material.sigma_t;
                }
            };
            //we add this part the 11th April
        };
        return f;
    };//return the local cohesive force, and the force in total

    void cc_length_col(Mesh mesh_total, const Material &material,
                       il::Array<double> width_large,
                       il::int_t i, il::int_t j,il::int_t c_i,il::int_t c_j,
                       il::Array<double> width_history, const int &p,
                       const il::Array2D<int> &id,
                       const il::Array2D<int> &col_matrix,
                       const il::int_t &dof_dim,il::io_t,
                       double &length_coh, double &crack_length ,double &energy_j) {

        crack_length = 0.;
        length_coh = 0.;
        energy_j=0.;
        //il::Array<double> width_etoc=width_edge_col(width,p);

        for (il::int_t s = 0; s < c_j - c_i + 1; ++s) {

            il::Array2D<int> col_row_s = search(col_matrix, int(s),il::io);
            SegmentCharacteristic segi = get_segment_DD_characteristic(
                    mesh_total, col_row_s(0, 0) + i, p);
            crack_length += 0.5 * segi.size;
        }



        for(il::int_t s1=i;s1<j+1;++s1){//s1 is the element number, starting from 0
            SegmentCharacteristic segi1 = get_segment_DD_characteristic(
                    mesh_total, s1, p);
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
                    energy_j+=(width_large[id(s1,dof_dim+1)]+
                            width_large[id(s1,1)])*material.sigma_t*0.5*segi1.size;//
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
                    energy_j+=(small+material.wc)*(material.wc-small)/
                            fabs(width_large[id(s1,dof_dim+1)] - width_large[id(s1,1)])
                              *material.sigma_t*0.5*segi1.size;//
                }
            };
        };
    }//now true to any shape of crack



    il::Array<double> initialwidth_col(const il::Array2D<double> &kmat,
                                       const il::Array2D<int> &id, int p,
                                       const Initial_condition &initial_condition,
                                       il::int_t i, il::int_t j,
                                       il::Status &status,
                                       const il::int_t &dof_dim) {
        il::int_t dof = dof_dim * (p + 1);
        il::Array2D<double> kini{dof * (j - i + 1), dof * (j - i + 1), 0.};
        take_submatrix(kini, id(i,0), id(j,dof_dim*1+1),
                       id(i,0),id(j,dof_dim*1+1), kmat);
        il::Array<double> F{dof * (j - i + 1),
                            -initial_condition.pini + initial_condition.sigma0};
        for (il::int_t m = 0; m < dof*(j- i+1); m += dof_dim) {
            F[m] = 0;
        };
        il::Array<double> delta_w_ini = il::linear_solve(kini, F, il::io, status);
        return delta_w_ini;

    };

    void construct_matrix_col(il::Array2D<double> &kmatC,
                              il::Array<double> vwc,
                              il::Array<double> cohf, const Material &material,
                              il::Array<double> unit,
                              const Initial_condition &initial_condition,
                              il::Status &status,
                              double pressure_f, il::Array<double> widthB,
                              const int &p, const il::int_t &dof_dim, il::io_t,
                              il::Array<double> &width, double &pressure) {
        il::int_t n = kmatC.size(0);
        il::Array2D<double> newmatrix{n + 1, n + 1, 0};
        il::Array<double> newvector{n + 1, 0};
        for (int m = 0; m < n; ++m) {
            for (int s = 0; s < n; ++s)
                newmatrix(s, m) = kmatC(s, m);
        };
        newmatrix(n, n) = material.U;
        il::Array<double> fK = il::dot(kmatC, widthB);
        il::Array<double> fP{n, pressure_f - initial_condition.sigma0};
        for (int c = 0; dof_dim * c < n; c++) {
            fP[2 * c] = 0;
        }

        for (int q = 0; q < n; ++q) {

            newmatrix(q, n) = unit[q];
            newmatrix(n, q) = vwc[q];
            newvector[q] = cohf[q] - fP[q] - fK[q];
        };
        newvector[n] = initial_condition.Q0 * initial_condition.timestep;
        il::Array<double> widthinter{n, 0.};
        il::Array<double> dd = il::linear_solve(newmatrix, newvector, il::io,
                                                status);
        for (int t = 0; t < n; ++t) {
            widthinter[t] = dd[t];
        };
        width = widthinter;
        pressure = dd[n];
        status.abort_on_error();
    }


    void plasticity_loop_col(const Material &material,
                             const Initial_condition &initial_condition,
                             il::int_t c_i, il::int_t c_j,
                             const il::Array2D<double> &kmat,
                             const il::Array<double> &vwc0,
                             const il::Array2D<int> &id,const int &p,
                             il::Array<double> widthB, double pressure_f,
                             il::Array<double> width_history,
                             const il::int_t &dof_dim,
                             const il::Array2D<int>&col_matrix, il::io_t,
                             il::Array<double> &delta_width,
                             double &pressure_change,
                             il::Array<double> &coht,int &mm) {
        il::Array2D<int> col_row_i=search(col_matrix,int(c_i),il::io);
        il::Array2D<int> col_row_j=search(col_matrix,int(c_j),il::io);

        il::Array2D<double> kmatnew{dof_dim * (c_j - c_i + 1), dof_dim * (c_j - c_i + 1), 0.};
        take_submatrix(kmatnew, id(col_row_i(0,0),dof_dim*col_row_i(0,1)),
                       id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1),
                       id(col_row_i(0,0),dof_dim*col_row_i(0,1)),
                       id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1),
                       kmat);

        il::Array<double> vwc{dof_dim * (c_j - c_i + 1), 0.};
        //il::Array<double> vwc0= get_vwc_vc_col(mesh_total, p,dof_dim);
        for (il::int_t m =0; m < id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))+1; ++m) {
            vwc[m] = vwc0[m+id(col_row_i(0,0),dof_dim*col_row_i(0,1))];
        }
        il::Array<double> unitm{dof_dim * (c_j - c_i + 1), 1.};
        for (il::int_t s = 0; s < dof_dim * (c_j - c_i + 1); s += dof_dim) {
            unitm[s] = 0.;
        }

        il::Array<double> delta_w_ini{widthB.size(), 0.};

        il::Array<double> coh{widthB.size(), 0.};

        int k = 0;
        double error_w = 1.;
        il::Array<il::Status> statusk{100,};
        double error_p=1.;
        double error_r=1;
        double inter_pressure=0.;


        il::Array<double> delta_w_bitera=delta_w_ini;


        while (k < 100 && (error_w > 0.00001 or error_p>0.00001 or error_r>0.000000001)) {
            il::Array<double> width_inm = widthB;
            il::blas(1.0, delta_w_ini, 1.0, il::io, width_inm);
            coh = cohesive_force_col(material, width_inm, col_row_i,
                                     col_row_j, width_history, p,dof_dim,id);

            construct_matrix_col(kmatnew, vwc, coh, material, unitm,
                                 initial_condition, statusk[k],
                                 pressure_f, widthB, p, dof_dim,il::io,
                                 delta_width, pressure_change);
            il::Array<double> intermediate = delta_width;
            il::blas(-1.0, delta_w_bitera, 1.0, il::io, intermediate);
            double error_inter = il::norm(intermediate, il::Norm::L2);
            error_w = error_inter / il::norm(delta_width, il::Norm::L2);

            il::blas(0.85, delta_width, 0.15, il::io, delta_w_ini);
            delta_w_bitera=delta_width;

            double error_inter_p= fabs(pressure_change-inter_pressure);
            error_p=error_inter_p/fabs(pressure_change);
            inter_pressure=pressure_change;

            error_r=fabs(il::dot(vwc,delta_width)+pressure_change*material.U-initial_condition.Q0*initial_condition.timestep);




            ++k;


        }
        mm=k;
        //from here you can consider to plot the current coh vector

        for (il::int_t r = 0; r < widthB.size(); r++) {
            coht[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+ r] = coh[r];
        }

    }

    il::Array2D<int> collocation_matrix(int &nelts, const int &p){
        il::Array2D<int> col{nelts, p+1, 0};
        int j;

        for (int i = 0; i < nelts; ++i) {
            j = i * (p + 1);
            for (int k = 0; k < (p + 1); ++k) {
                col(i, k) = j + k;
            }
        }
        return col;
    }

    void
    propagation_loop_col(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                         const Material &material,
                         const Initial_condition &initial_condition,
                         il::int_t i0, il::int_t j0, int nstep,
                         il::Status &status, il::io_t,
                         il::Array2C<double> &widthlist,
                         il::Array<double> &plist, il::Array<double> &l_coh,
                         il::Array<double> &l, il::Array2C<double> &coh_list,
                         il::Array<int> &mvalue,int &break_time,
                         il::Array2C<double> &stress_list, il::Array<double> &energy_g) {

        int n = mesh_total.nelts();
        IL_EXPECT_FAST(n == id.size(0));
        il::Array2D<int> col_matrix=collocation_matrix(n,p);
        il::int_t dof_dim=id.size(1)/2;
        il::int_t dof = dof_dim * (p + 1);
        il::int_t ttn=id.size(0)*id.size(1);
        il::Array2D<double> kmat{ttn,ttn,0.};
        kmat=basic_assembly(mesh_total, id, p, material.Ep);
        il::Array<double> width_history{ttn, 0.};

        //i0 and j0 are the beginning element number in C++, starting from 0,1,2.....
        for (il::int_t h = i0; h < j0 + 1; ++h) {
            width_history[id(h,dof_dim*0+1)] = 1.;
            width_history[id(h,dof_dim*1+1)] = 1.;
        };

        il::Array<double> vwc0= get_vwc_vc_col(mesh_total, p,dof_dim);


        il::Array<double> widthB = initialwidth_col(kmat, id, p,
                                                    initial_condition, i0, j0,
                                                    status,dof_dim);
        il::Array<double> width;
        il::Array2C<double> widthlistinter{nstep, ttn,0.};
        il::Array2C<double> stress_profile{nstep,ttn,0.};
        il::Array<double> plistinter{nstep + 1, 0.};
        il::Array<int> m_value{nstep,0};
        int mk=0;
        plistinter[0] = initial_condition.pini;
        int s = 0;
        int m = 0;
        il::int_t i = i0;
        il::int_t j = j0;
        double pressure_change = 0.;//not sure we should initialize here or not
        double pressure = initial_condition.pini;
        double time = 0;
        bool break_value=false;
        bool break_value_2=false;
        il::Array<double> crack_length{nstep+1, 0.};
        il::Array<double> length_coh{nstep, 0.};
        il::Array<double> energy_j{nstep, 0.};

        il::Array2C<double> coh_tt{nstep, ttn, 0.};
        il::Array<double> coht{ttn, 0.};

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
                plasticity_loop_col(material, initial_condition,c_i, c_j,
                                    kmat, vwc0, id, p, widthB, pressure,
                                    width_history,dof_dim,col_matrix,il::io,
                                    delta_width, pressure_change, coht,mk);
                m_value[s]=mk;

                il::Array<il::int_t> index = stress_criteria_col(kmat, id, material,
                                                                 initial_condition,
                                                                 widthB, delta_width,
                                                                 i, j, c_i,c_j,
                                                                 dof_dim,col_matrix,p);

                if (c_i == index[0] and c_j == index[1]) {
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    break;
                }
                if (index[0] <= 0 or index[1] >= dof_dim*n - 1) {//not sure the 2 here is (p+1) or dof_dim
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);

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
                for (int dwn = 0; dwn < widthB.size(); dwn += 1) {
                    delta_width_new[dof_dim * (c_i - index[0]) +
                                    dwn] = delta_width[dwn];
                }
                widthB.resize(dof_dim * (index[1] - index[0] + 1));
                widthB = widthBnew;
                width = widthBnew;


                c_i = index[0];
                c_j = index[1];
                i=index[2];
                j=index[3];


                il::blas(1., delta_width_new, 1.0, il::io, width);
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


            cc_length_col(mesh_total, material, width_large, i, j, c_i,c_j,
                          width_history, p,id,col_matrix,dof_dim, il::io,
                          length_coh[s], crack_length[s+1],energy_j[s]);

            width_history = write_history_col(width, width_history, c_i,c_j,
                                              col_matrix, material,p,dof_dim,id);
            //Here we use width_large, could be simpler for the coding

            for (int cohi = 0; cohi < dof * n; cohi++) {
                coh_tt(s, cohi) = coht[cohi];
            }

            il::Array<double> stress_current=il::dot(kmat,width_large);
            for (int stre=0;stre<dof*n;stre++){
                stress_profile(s,stre)=stress_current[stre];
            }


            pressure += pressure_change;
            plistinter[s + 1] = pressure;//pressure output
            widthB = width;
            time += initial_condition.timestep;
            ++s;

        }
        widthlist = widthlistinter;
        plist = plistinter;
        l = crack_length;
        l_coh = length_coh;
        coh_list = coh_tt;
        energy_g=energy_j;
        mvalue=m_value;
        stress_list=stress_profile;
        //to output the iteration times in plasticity_loop to get the solution
        status.abort_on_error();
    }

    void energy_output(il::Array2C<double> widthlist,il::Array<double> plist,
                       il::Array<double> l_c,il::Array<double> l_coh,
                       Material material,Mesh mesh_total,il::Array2D<int> &id,
                       const int &p,const int &dof_dim, Initial_condition initial_condition,
                       il::io_t,
                       il::Array<double> &energy_ff,
                       il::Array<double> &energy_coh,
                       il::Array<double> &energy_j_int){
        il::int_t n=widthlist.size(0);
        il::int_t nc=widthlist.size(1);
        //il::Array<double> energy_p {n,0.};
        il::Array<double> energy_f {n,0.};
        il::Array<double> energy_coh_k{n,0.};
        //il::Array<double> energy_g{n,0.};
        il::Array<double> energy_j_integral{n+1,0.};



        for(int s=0;s<n;s++){
            energy_f[s]=3.1415926*l_c[s+1]*plist[s+1]*plist[s+1]/material.Ep;
            energy_coh_k[s]=material.sigma_t*material.sigma_t*8/3.1415926/material.Ep*l_coh[s];

            energy_j_integral[s+1]=0.5*initial_condition.Q0*plist[s+1]*initial_condition.timestep;



//            if(l_c[s+1]-l_c[s]>l_coh[s]){
//                if(s==0){
//                    energy_g[s]=energy[s]+material.sigma_t*material.wc*(l_c[s+1]-l_c[s]-l_coh[s]);
//                }
//                else{
//                    energy_g[s]=energy[s]+material.sigma_t*material.wc*(l_c[s+1]-l_c[s]+l_coh[s-1]-l_coh[s])-energy[s-1];
//                }
//            }
//            else{
//                if(s==0){
//                    energy_g[s]=energy[s];
//                }
//                else{
//                    energy_g[s]=energy[s]-energy[s-1]+material.sigma_t*material.wc*(l_c[s+1]-l_c[s]+l_coh[s-1]-l_coh[s]);
//                }
//
//            }
            for(int sc=1;sc<nc;sc+=dof_dim){
                il::Array2D<int> ne=search(id,sc,il::io);
                SegmentCharacteristic segi = get_segment_DD_characteristic(
                        mesh_total, ne(0,0), p);

                if(s==0){
                    //energy_p[s]+=1./2.*widthlist(s,sc)*plist[s+1]*0.5*segi.size;
                    energy_j_integral[s+1]-=0.5*(plist[s+1]*widthlist(s,sc))*0.5*segi.size;
                    if(widthlist(s,sc)<material.wc && widthlist(s,sc)>0.){
                        energy_j_integral[s+1]-=0.5*material.sigma_t*widthlist(s,sc)*0.5*segi.size;
                    }

                }
                else{
                    //energy_p[s]+=1./2.*(2*plist[s+1]*widthlist(s,sc)-plist[s+1]*widthlist(s-1,sc)-plist[s]*widthlist(s,sc))*0.5*segi.size;
                    energy_j_integral[s+1]-=0.5*(plist[s+1]-plist[s])*widthlist(s,sc)*0.5*segi.size;
                    if(widthlist(s,sc)<material.wc && widthlist(s,sc)>0.){
                        energy_j_integral[s+1]-=0.5*material.sigma_t*(widthlist(s,sc)-widthlist(s-1,sc))*0.5*segi.size;
                    }
//                    if(widthlist(s-1,sc)>0.){
//                        if(l_c[s+1]-l_c[s]<=l_coh[s] && widthlist(s-1,sc)<material.wc){
//                            if(widthlist(s,sc)<material.wc){
//                                energy_g[s]-=material.sigma_t*widthlist(s-1,sc)*0.5*segi.size;
//                            }
//                            else{
//                                energy_g[s]+=material.sigma_t*(material.wc-widthlist(s-1,sc))*0.5*segi.size;
//                            }
//                        }
//                    }
                }
            }
            double energy_middle=energy_j_integral[s+1];
            //this is the energy ratio //energy_j_integral[s+1]=energy_middle/(l_c[s+1]-l_c[s])/material.sigma_t/material.wc;
            energy_j_integral[s+1]=energy_middle-(l_c[s+1]-l_c[s])*material.sigma_t*material.wc;
            energy_j_integral[s+1]=energy_j_integral[s+1]/initial_condition.Q0/initial_condition.timestep/material.sigma_t;

            //to calculate the slope, corresponding to the energy release rate
//            if(s==0){
//                energy_f[s]=energy_p[s];
//            }
//            else{
//                energy_f[s]=energy_p[s]-energy_p[s-1];
//            }

        }
        energy_ff=energy_f;
        //energy_pp=energy_p;
        energy_coh=energy_coh_k;
        //energy_gg=energy_g;
        energy_j_int=energy_j_integral;
    }
    il::Array<double> volume_output(il::Array2C<double> &widthlist, const int &dof_dim){
        il::int_t n=widthlist.size(0);
        il::int_t nc=widthlist.size(1);
        il::Array<double> volume_change{n,0.};
        for(int s=0;s<n;s++){
            for(int sc=1;sc<nc;sc+=dof_dim){
                if(s==0){
                    volume_change[s]+=widthlist(s,sc);
                }
                else{
                    volume_change[s]+=widthlist(s,sc)-widthlist(s-1,sc);
                }
            }
        }
        return volume_change;
    }

    il::Array2C<double> deal_with_stress(il::Array2C<double> stresslist,
                                         il::Array2C<double> cohlist,
                                         il::Array<double> plist,
                                         Initial_condition initial_condition){
        il::Array2C<double> stresslist_new=stresslist;
        il::int_t n=stresslist.size(0);
        il::int_t nc=stresslist.size(1);
        for(il::int_t s=0;s<n;s++){
            for(il::int_t sc=0;sc<nc;sc++){
                if(cohlist(s,sc)!=0){
                    stresslist_new(s,sc)=stresslist(s,sc)+plist[s+1]-initial_condition.sigma0;
                    //because the stress=K*w+pf-confining pressure, and this is only true for the stress inside te crack
                }
            }
        }
        return stresslist_new;
    }



    void get_xlist_col(il::Array<double> &xlist, Mesh mesh_total) {
        int n = mesh_total.nelts();
        for (int i = 0; i < n; ++i) {
            xlist[2 * i] = mesh_total.coor()(i, 0);
            xlist[2 * i + 1] = mesh_total.coor()(i + 1, 0);
        }
    };

}
//}