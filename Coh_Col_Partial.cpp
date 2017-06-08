//
// Created by DONG LIU on 2/23/17.
//

#include "Coh_Col_Partial.h"

#include<il/norm.h>



namespace hfp2d {

    il::Array2D<double>partial_matrix(il::Array<double> coh_f){
        il::int_t n=coh_f.size();
        il::Array2D<double> partial_m{n,n,0.};
        for(il::int_t m=0;m<n;m++){
            partial_m(m,m)=1.0;
            if(coh_f[m]!=0.){
                partial_m(m,m)=0.;
            }
        }
        return partial_m;
    }



    void construct_matrix_col_partial(il::Array2D<double> &kmatC,
                              il::Array2D<double> vwc,
                              il::Array<double> cohf, const Material &material,
                              il::Array<double> unit,
                              const Initial_condition &initial_condition,
                              il::Status &status,
                              double pressure_f, il::Array<double> widthB,
                              const int &p, const il::int_t &dof_dim,
                              il::Array2D<double> partial_m, il::io_t,
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
        il::Array<double> fP{n, pressure_f};

        for (int c = 0; dof_dim * c < n; c++) {
            fP[2 * c] = initial_condition.sigma0;//we minus here and add after to make the shear disappear
        };

        il::Array2D<double> vwc_p=il::dot(vwc,partial_m);
        il::Array<double> unit_p=il::dot(partial_m,unit);
        il::Array<double> fP_p=il::dot(partial_m,fP);




        for (int q = 0; q < n; ++q) {

            newmatrix(q, n) = unit_p[q];
            newmatrix(n,q)=vwc(0,q);
//            newmatrix(n, q) = vwc_p(0,q);
            newvector[q] = cohf[q] - fP_p[q] - fK[q]+initial_condition.sigma0;
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


    void plasticity_loop_col_partial(const Material &material,
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

        il::Array2D<double> vwc{1,dof_dim * (c_j - c_i + 1), 0.};
        //il::Array<double> vwc0= get_vwc_vc_col(mesh_total, p,dof_dim);
        for (il::int_t m =0; m < id(col_row_j(0,0),dof_dim*col_row_j(0,1)+dof_dim-1)-id(col_row_i(0,0),dof_dim*col_row_i(0,1))+1; ++m) {
            vwc(0,m) = vwc0[m +id(col_row_i(0,0),dof_dim*col_row_i(0,1))];
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


        while (k < 100 && (error_w > 0.00001 or error_p>0.00001 or error_r>0.000000001)) {//the number of k should be correlated to the statusk
            il::Array<double> width_inm = widthB;
            il::blas(1.0, delta_w_ini, 1.0, il::io, width_inm);
            coh = cohesive_force_col(material, width_inm, col_row_i,
                                     col_row_j, width_history, p,dof_dim,id);

            il::Array2D<double> partial_m=partial_matrix(coh);//new

            construct_matrix_col_partial(kmatnew, vwc, coh, material, unitm,
                                 initial_condition, statusk[k],
                                 pressure_f, widthB, p, dof_dim,partial_m,il::io,
                                 delta_width, pressure_change);//construct_matrix needs change also
            il::Array<double> intermediate = delta_width;
            il::blas(-1.0, delta_w_bitera, 1.0, il::io, intermediate);
            double error_inter = il::norm(intermediate, il::Norm::L2);
            error_w = error_inter / il::norm(delta_width, il::Norm::L2);

            il::blas(0.85, delta_width, 0.15, il::io, delta_w_ini);
            delta_w_bitera=delta_width;

            double error_inter_p= fabs(pressure_change-inter_pressure);
            error_p=error_inter_p/fabs(pressure_change);
            inter_pressure=pressure_change;


            error_r=fabs(il::dot(vwc,delta_width)[0]+pressure_change*material.U-initial_condition.Q0*initial_condition.timestep);




            ++k;


        }
        mm=k;
        //from here you can consider to plot the current coh vector

        for (il::int_t r = 0; r < widthB.size(); r++) {
            coht[id(col_row_i(0,0),dof_dim*col_row_i(0,1))+ r] = coh[r];
        }

    }


    void
    propagation_loop_col_partial(Mesh mesh_total, il::Array2D<int> &id, const int &p,
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
                plasticity_loop_col_partial(material, initial_condition,c_i, c_j,
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


    void energy_output_partial(il::Array2C<double> widthlist,il::Array<double> plist,
                       il::Array<double> l_c,il::Array<double> l_coh, il::Array2C<double> cohlist,
                       Material material,Mesh mesh_total,il::Array2D<int> &id,
                       const int &p,const int &dof_dim, Initial_condition initial_condition,
                       il::io_t,
                       il::Array<double> &energy_ff,
                       il::Array<double> &energy_coh,
                       il::Array<double> &energy_j_int){
        il::int_t n=widthlist.size(0);
        il::int_t nc=widthlist.size(1);
        il::Array<double> energy_f {n,0.};
        il::Array<double> energy_coh_k{n,0.};
        il::Array<double> energy_j_integral{n+1,0.};

        for(int s=0;s<n;s++){
            energy_f[s]=3.1415926*l_c[s+1]*plist[s+1]*plist[s+1]/material.Ep;
            energy_coh_k[s]=material.sigma_t*material.sigma_t*8/3.1415926/material.Ep*l_coh[s];
            energy_j_integral[s+1]=initial_condition.Q0*plist[s+1]*initial_condition.timestep;

            for(int sc=1;sc<nc;sc+=dof_dim){
                il::Array2D<int> ne=search(id,sc,il::io);
                SegmentCharacteristic segi = get_segment_DD_characteristic(
                        mesh_total, ne(0,0), p);
                if(cohlist(s,sc)==0){
                    if(s==0){
                        energy_j_integral[s+1]-=0.5*(plist[s+1]*widthlist(s,sc))*0.5*segi.size;
                    }
                    else{
                        energy_j_integral[s+1]-=0.5*(plist[s+1]*widthlist(s,sc)-plist[s]*widthlist(s-1,sc))*0.5*segi.size;
                    }
                }
            }

            double energy_middle=energy_j_integral[s+1];

            energy_j_integral[s+1]=energy_middle-(l_c[s+1]-l_c[s])*material.sigma_t*material.wc;

        }
        energy_ff=energy_f;
        energy_coh=energy_coh_k;
        energy_j_int=energy_j_integral;
    }


}
//}