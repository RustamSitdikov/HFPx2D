//
// Created by DONG LIU on 2/23/17.
//

#include "Coh_Prop_new.h"
#include<il/norm.h>
#include <il/math.h>


//void take_submatrixint(il::Array2D<int>& sub, int i0, int i1, int j0, int j1,
//                       const il::Array2D<int>& A){
//
//    for(int i = i0; i <= i1;++i) {
//        for (int j=j0; j<= j1;++j){
//            sub(i-i0,j-j0)=A(i,j);
//        }
//
//    }
//}

//Mesh take_submesh(int n1, int n2, Mesh mesh_total){//we consider n1 and n2
//// represents for the element number
//    Mesh meshC;
//    il::Array2D<int> connC{n2-n1+1,2,0};
//    take_submatrixint(connC,0,n2-n1,0,1,mesh_total.conn);//here also, the connectivity begins from 0.1.2.3...
//    il::Array2D<double> CoorC{n2-n1+2,2,0.};
//    take_submatrix(CoorC,n1-1,n2,0,1,mesh_total.Coor);
//    meshC.set_values(CoorC,connC);
//    return meshC;
//};

//void material_condition(Material &material, Initial_condition &initial_condition,
//                        double wc1, double sigma_t1,double Ep1,double U1,
//                        double pini1,double Q01,double timestep1,double sigma01){
//    material.Ep=Ep1;
//    material.sigma_t=sigma_t1;
//    material.U=U1;
//    material.wc=wc1;
//    initial_condition.pini=pini1;
//    initial_condition.Q0=Q01;
//    initial_condition.sigma0=sigma01;
//    initial_condition.timestep=timestep1;
//}

//2*(p+1)

//element number i begins from 1.2.3.4...

//namespace hfp_opening{

namespace hfp2d {

    void
    material_condition(const double& wc1, const double& sigma_t1, const double& Ep1, const double& U1,
                       const double& pini1, const double& Q01, const double& timestep1,
                       const double& sigma01,il::io_t,Material& material, Initial_condition& initial_condition) {
        material.Ep = Ep1;
        material.sigma_t = sigma_t1;
        material.U = U1;
        material.wc = wc1;
        initial_condition.pini = pini1;
        initial_condition.Q0 = Q01;
        initial_condition.sigma0 = sigma01;
        initial_condition.timestep = timestep1;
    }

    void get_vwc_vc(il::Array<double> &vwc, Mesh mesh_total, int p) {

        int n = mesh_total.nelts();
        int dof = 2 * (p + 1);
        il::Array<double> vwcinter{dof * n, 0.};
        for (int i = 0; i < n; ++i) {
            //il::Array2D<double> xi{2, 2, 0.};
            //take_submatrix(xi, i, i + 1, 0, 1, mesh_total.coor());


            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SegmentCharacteristic segi = get_segment_DD_characteristic(mesh_total,i, p);
            vwcinter[dof * i] = 0.;
            vwcinter[dof * i + 1] = 0.5 * segi.size;
            vwcinter[dof * i + 2] = 0.;
            vwcinter[dof * i + 3] = 0.5 * segi.size;
        };
        vwc = vwcinter;
    }


    il::Array<il::int_t> stress_criteria_new(const il::Array2D<double> &kmat,
                                             const il::Array2D<int> &id,
                                             const Material &material,
                                             const Initial_condition &initial_condition,
                                             il::Array<double> widthB,
                                             il::Array<double> delta_w,
                                             il::int_t i, il::int_t j,
                                             int p) {//here we need to get the collocation points of the next-possible failed elements not those of the activated elements
        int dof = 2 * (p + 1);
        il::Array2D<double> kt{dof * id.size(0), dof * (j - i + 1), 0.};
        take_submatrix(kt, 0, dof * id.size(0) - 1, dof * i - dof,
                       dof * j - dof + 3, kmat);
        il::blas(1.0, delta_w, 1.0, il::io, widthB);
        il::Array<double> syy = il::dot(kt, widthB);
        double cr = 0;
        double cl = 0;

        for (il::int_t r = 1; r < i; r++) {
            if (syy[dof * r - dof + 1] >=
                material.sigma_t + initial_condition.sigma0 &&
                syy[dof * r - dof + 3] >=
                material.sigma_t + initial_condition.sigma0 and i >= 1) {
                i = r;//the
            }
        }
        for (il::int_t r2 = id.size(0); r2 > j; r2--) {
            if (syy[dof * r2 - dof + 1] >=
                material.sigma_t + initial_condition.sigma0 &&
                syy[dof * r2 - dof + 3] >=
                material.sigma_t + initial_condition.sigma0 and
                j <= id.size(0)) {
                j = r2;
            }
        }  //with this, we can accelerate the propagation, but we suppose the far point stress is smaller than the near point stress

        //below, it's the propagation one by one strategy.


//    int dof=2*(p+1);
//    il::Array2D<double> kl{dof,dof*(j-i+1),0.};
//    il::Array2D<double> kr{dof,dof*(j-i+1),0.};
//    take_submatrix(kl,dof*i-dof-dof,dof*i-dof-1,dof*i-dof,dof*j+3-dof,kmat);//needs to be verified
//    take_submatrix(kr,dof*j,dof*j+3,dof*i-dof,dof*j+3-dof,kmat);
//    il::blas(1.0,delta_w,1.0,il:: io,widthB);
//
//
//    double syyl1=il::dot(kl,widthB)[1];
//    double syyl2=il::dot(kl,widthB)[3];
//    double syyr1=il::dot(kr,widthB)[1];
//    double syyr2=il::dot(kr,widthB)[3];
//
//    //we change this the 13th April
//    if((syyl1>= (material.sigma_t+initial_condition.sigma0) and
//        syyl2>=(material.sigma_t+initial_condition.sigma0)) and i>=1){
//        --i;
//    };
//    if((syyr1>= (material.sigma_t+initial_condition.sigma0) and
//        syyr2>=(material.sigma_t+initial_condition.sigma0)) and j<=id.size(0)){
//        //needs to verify this part
//        ++j;
//    };
//




//    if((syyl1>= (material.sigma_t+initial_condition.sigma0) or
//        syyl2>=(material.sigma_t+initial_condition.sigma0)) and i>=1){
//        --i;
//    };
//    if((syyr1>= (material.sigma_t+initial_condition.sigma0) or
//        syyr2>=(material.sigma_t+initial_condition.sigma0)) and j<=id.size(0)){
//        //needs to verify this part
//        ++j;
//    };
        if (j >= id.size(0)) {
            j = id.size(0);
        };
        if (i < 1) {
            i = 1;
        };
        il::Array<il::int_t> index{2, 0};
        index[0] = i;
        index[1] = j;
        return index;

    };


//return index, the first is the new left element number i, the second is the
//new right element number j of the crack

//il::Array<double> former_pressure(double pressure_f,int i,int j,int i0,int j0){
//    il::Array<double> fP{4*(j-i+1),0.};
//    if(i<i0){
//        for(int m=1;m<4*(i0-i);m+=2){
//            fP[m]=-pressure_f;
//        }
//    }
//    if(j>j0){
//        for(int s=4*(j0-i+1)+1;s<4*(j-i+1);s+=2){
//            fP[s]=-pressure_f;
//        }
//    }
//    return fP;
//}
//
//il::Array<double> former_k(il::Array2D<double> kmatC,int i,int j,
//                           int i0,int j0,il::Array<double> widthB){
//    il::Array<double> fK{4*(j-i+1),0.};
//    il::Array<double> vectem=il::dot(kmatC,widthB);
//    if(i<i0){
//        for (int k = 0; k <4*(i0-i) ; ++k) {
//            fK[k]=-vectem[k];
//        }
//    }
//    if(j>j0){
//        for (int s = 4*(j0-i+1); s <4*(j-i+1) ; ++s) {
//            fK[s]=-vectem[s];
//        }
//    }
//    return fK;
//}


//void construct_matrix(il::Array<double> &width, double &pressure,
//                      il::Array2D<double> kmatC,il::Array<double> vwc,
//                      il::Array<double> cohf,Material material,
//                      il::Array<double> unit,Initial_condition initial_condition,
//                      il::Status &status,int i,int j,int i0,int j0,
//                      double pressure_f,il::Array<double> widthB){
//    int n =int(kmatC.size(0));
//    il::Array2D<double> newmatrix{n+1,n+1,0};
//    il::Array<double> newvector{n+1,0};
//    for(int m=0;m<n;++m){
//        for(int s=0;s<n;++s)
//            newmatrix(m,s)=kmatC(m,s);
//    };
//    newmatrix(n,n)=material.U;
//    il::Array<double> fP=former_pressure(pressure_f,i,j,i0,j0);
//    il::Array <double> fK=former_k(kmatC,i,j,i0,j0,widthB);//added the 6th Februray
//
//
//    for(int q=0;q<n;++q){
//        newmatrix(q,n)=unit[q];
//        newmatrix(n,q)=vwc[q];
//        newvector[q]=cohf[q]+fP[q]+fK[q];
//    };
//    newvector[n]=initial_condition.Q0 * initial_condition.timestep;
//    il::Array<double> widthinter{n,0.};
//    il::Array<double> dd = il::linear_solve(newmatrix,newvector, il::io, status);
//    for(int t=0;t<n;++t){
//        widthinter[t]=dd[t];
//    };
//    width=widthinter;
//    pressure=dd[n];
//    status.abort_on_error();//every time we execuate the constructmatrix,
//    // we will get the status_
//    // status.abort_on_error();//we can add this line to do the same thing,
//    // so actually you don't have to change the "Status.h"
//}
//
//void get_vwc_vc(il::Array<double> &vwc,Mesh mesh_total,int p){
//
//    int n=mesh_total.nelts();
//    il::Array<double> vwcinter{4*n,0.};
//    for (int i=0;i<n;++i){
//        il::Array2D<double> xi{2,2,0.};
//        take_submatrix(xi,i,i+1,0,1,mesh_total.Coor);
//        SegmentCharacteristic segi=get_segment_DD_characteristic(xi,p);
//        vwcinter[4*i]=0.;
//        vwcinter[4*i+1]=0.5*segi.size;
//        vwcinter[4*i+2]=0.;
//        vwcinter[4*i+3]=0.5*segi.size;
//    };
//    vwc=vwcinter;
//}

    il::Array<double> width_edge_col(il::Array<double> width, int p) {
        il::int_t n = width.size();
        il::int_t dof = 2 * (p + 1);
        il::Array<double> width_new{n, 0.};
        for (il::int_t i = 0; i < n / dof; ++i) {
            width_new[dof * i + 1] = (0.5 + sqrt(2) / 4.) * width[dof * i + 1] +
                                     (0.5 - sqrt(2) / 4.) * width[dof * i + 3];
            width_new[dof * i + 3] = (0.5 - sqrt(2) / 4.) * width[dof * i + 1] +
                                     (0.5 + sqrt(2) / 4.) * width[dof * i + 3];
        }
        return width_new;
    }


    il::Array<double>
    write_history(il::Array<double> width, il::Array<double> width_history,
                  il::int_t i, const Material &material, int p) {

        il::int_t n = width.size();
        il::int_t n_t = width_history.size();
        int dof = 2 * (p + 1);
        il::Array<double> width_history_new{n_t, 0.};
        il::Array<double> width_new{n, 0.};
        width_new = width_edge_col(width, p);
        //changes the 7th April, m initial value from 0 to 1 and the varying step from 1 to 2, which makes the width_history only relates to the opening not the slip.
        for (il::int_t m = 1; m < n; m += 2) {
            width_history_new[dof * i - dof + m] = width_history[dof * i - dof +
                                                                 m];
            if (width_new[m] > width_history[dof * i - dof + m]) {
                //if(width_new[m]<material.wc && width_history[dof*i-dof+m]<1.){
                width_history_new[dof * i - dof + m] = width_new[m];
            }
        }
        return width_history_new;
    }


    il::Array<double> cohesive_force_new(const Material &material,
                                         il::Array<double> width, il::int_t i,
                                         il::int_t j,
                                         il::Array<double> width_history,
                                         int p) {
        il::Array<double> f{width.size(), 0.};
        int dof = 2 * (p + 1);
        for (il::int_t s = 1; s < width.size(); s += 2) {
            il::Array<double> width_etoc = width_edge_col(width, p);
            if (width_etoc[s] < material.wc and
                width_history[dof * i - dof + s] ==
                0.) {//and width_history[dof*i-dof+s]==0.
                // to keep the unchanged version//if(0.<width_etoc[s] and width_etoc[s]<material.wc and width_history[dof*i-dof+s]<1.)
                f[s] = material.sigma_t;
            };

            if (width_etoc[s] > 0 and width_etoc[s] < material.wc and
                width_history[dof * i - dof + s] > 0. and
                width_history[dof * i - dof + s] <
                material.wc) {//and width_history[dof*i-dof+s]<material.wc
                if (width_etoc[s] < width_history[dof * i - dof + s]) {
                    f[s] = width_etoc[s] / width_history[dof * i - dof + s] *
                           material.sigma_t;
                } else {
                    f[s] = material.sigma_t;
                }
            }; //we add this part the 11th April

        };
        return f;
    };//return the local cohesive force, and the force in total

    void cc_length(Mesh mesh_total,
                   const Material &material, il::Array<double> width,
                   il::int_t i, il::int_t j,
                   il::Array<double> width_history, int p, il::io_t,
                   double &length_coh, double &crack_length) {
        int n = mesh_total.nelts();
        int dof = 2 * (p + 1);
        crack_length = 0.;
        length_coh = 0.;
        //il::Array<double> width_etoc=width_edge_col(width,p);
        double m = 0.;
        for (il::int_t s = 0; s < j - i + 1; ++s) {
            //il::Array2D<double> xi{2, 2, 0.};
            //take_submatrix(xi, s + i - 1, s + i, 0, 1, mesh_total.coor());

            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!should be correct
            SegmentCharacteristic segi = get_segment_DD_characteristic(mesh_total,s+i, p);
            crack_length += segi.size;
            if ((width[dof * s + 1] < material.wc and
                 width_history[dof * i - dof + dof * s + 1] < material.wc)
                or (width[dof * s + 3] < material.wc and
                    width_history[dof * i - dof + dof * s + 3] <
                    material.wc)) {// if we don't apply cohesive force for the negative opening part, we don't count this part as cohesive length
                if (width[dof * s + 1] < material.wc and
                    width[dof * s + 3] < material.wc) {
                    length_coh += segi.size;
                } else {
                    double small = 0.;
                    if (width[dof * s + 1] > width[dof * s + 3]) {
                        small = width[dof * s + 3];
                    } else {
                        small = width[dof * s + 1];
                    };
                    length_coh += segi.size * (material.wc - small) /
                                  fabs(width[dof * s + 3] - width[dof * s + 1]);
                }
            };
        };
    }//now true to any shape of crack



    il::Array<double> initialwidth_new(const il::Array2D<double> &kmat,
                                       const il::Array2D<int> &id, int p,
                                       const Material &material,
                                       const Initial_condition &initial_condition,
                                       il::int_t i, il::int_t j,
                                       il::Status &status) {
        int dof = 2 * (p + 1);
        il::Array2D<double> kini{dof * (j - i + 1), dof * (j - i + 1), 0.};
        take_submatrix(kini, dof * i - dof, dof * j + 3 - dof, dof * i - dof,
                       dof * j + 3 - dof, kmat);
        il::Array<double> F{dof * j - dof * i + dof,
                            -initial_condition.pini + initial_condition.sigma0};
        for (il::int_t m = 0; m < dof * j - dof * i + dof; m += p + 1) {
            F[m] = 0;
        };
        il::Array<double> delta_w_ini = il::linear_solve(kini, F, il::io,
                                                         status);
        return delta_w_ini;

    };

    void construct_matrix_new(il::Array2D<double> &kmatC,
                              il::Array<double> vwc,//here not sure whether there should be a reference for kmatC
                              il::Array<double> cohf, const Material &material,
                              il::Array<double> unit,
                              const Initial_condition &initial_condition,
                              il::Status &status, il::int_t i, il::int_t j,
                              il::int_t i0, il::int_t j0,
                              double pressure_f, il::Array<double> widthB,
                              int p, il::io_t, il::Array<double> &width,
                              double &pressure) {
        il::int_t n = kmatC.size(0);
        il::Array2D<double> newmatrix{n + 1, n + 1, 0};
        il::Array<double> newvector{n + 1, 0};
        for (int m = 0; m < n; ++m) {
            for (int s = 0; s < n; ++s)
                newmatrix(s, m) = kmatC(s, m);
        };
        newmatrix(n, n) = material.U;
        il::Array<double> fK = il::dot(kmatC, widthB);
        il::Array<double> fP{n, pressure_f + initial_condition.sigma0};
        for (int c = 0; 2 * c < n; c++) {
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


    void plasticity_loop_new(const Material &material,
                             const Initial_condition &initial_condition,
                             il::int_t i, il::int_t j, il::int_t i0,
                             il::int_t j0,
                             const il::Array2D<double> &kmat, Mesh mesh_total,
                             const il::Array2D<int> &id, int p,
                             il::Array<double> widthB, double pressure_f,
                             il::Array<double> width_history, il::io_t,
                             il::Array<double> &delta_width,
                             double &pressure_change,
                             il::Array<double> &coht) {// i and j stands for the element number starting from 1
        int dof = 2 * (p + 1);
        il::Array2D<double> kmatnew{dof * (j - i + 1), dof * (j - i + 1), 0.};
        take_submatrix(kmatnew, dof * i - dof, dof * j + 3 - dof, dof * i - dof,
                       dof * j + 3 - dof,
                       kmat);//take_submatrix() function, the line number index is the input;
        il::Array<double> vwc0;
        il::Array<double> vwc{dof * (j - i + 1), 0.};
        get_vwc_vc(vwc0, mesh_total, p);
        for (il::int_t m = dof * i; m < dof * j + dof; ++m) {
            vwc[m - dof * i] = vwc0[m];
        }
        il::Array<double> unitm{dof * (j - i + 1), 1.};
        for (il::int_t s = 0; s < dof * (j - i + 1); s += 2) {
            unitm[s] = 0.;
        }
        il::Array<double> delta_w_ini{widthB.size(), 0.};
        //il::Array<double> coh_ini{widthB.size(),0.};
        //!!!!!!!!!!!!!!!!!!!!!!!!!!! If every time we have a zero coheisive force, we could have a normal shape of the crack!!!!!!!!!!!!!!Needs to find out the reason
        //coh_ini=cohesive_force_new(material,widthB,i,j,i0,j0,old_history,p);//here we have a problem, this cohesive force should use the width-history of last time
////when the fracture propagates, the size changes, the zeros added in vector widthB won't give a positive cohesive force, because with these zeros we always have width_edoc[s] is positive!!
//    if(i<i_inter){
//        for(il::int_t h=i;h<i_inter;h++){
//            coh_ini[dof*(h-i)+1]=0.;
//            coh_ini[dof*(h-i)+3]=0.;
//        }
//    }
//    if(j>j_inter){
//        for(il::int_t h2=j_inter;h2<j;h2++){
//            coh_ini[dof*(h2-i+1)+1]=0.;
//            coh_ini[dof*(h2-i+1)+3]=0.;
//        }
//    }
        il::Array<double> coh{widthB.size(), 0.};

        int k = 0;
        double error_w = 1.;
        il::Array<il::Status> statusk{100,};
        double positivechecker = 0.;

        while (k < 100 && error_w > 0.00001) {
            il::Array<double> width_inm = widthB;
            il::blas(1.0, delta_w_ini, 1.0, il::io, width_inm);
            coh = cohesive_force_new(material, width_inm, i, j, width_history,
                                     p);
            il::Array<double> coh_inm = coh;
            //il::blas(-1.0,coh_ini,1.0,il::io,coh_inm);

            construct_matrix_new(kmatnew, vwc,
                                 coh_inm, material, unitm, initial_condition,
                                 statusk[k],
                                 i, j, i0, j0, pressure_f, widthB, p, il::io,
                                 delta_width, pressure_change);
            il::Array<double> intermediate = delta_width;
            il::blas(-1.0, delta_w_ini, 1.0, il::io, intermediate);
            double error_inter = il::norm(intermediate, il::Norm::L2);
            error_w = error_inter / il::norm(delta_width, il::Norm::L2);
            il::blas(0.85, delta_width, 0.15, il::io, delta_w_ini);
            ++k;


        }
        //from here you can consider to plot the current coh vector

        for (il::int_t r = 0; r < widthB.size(); r++) {
            coht[dof * i - dof + r] = coh[r];
        }

    }



//void plasticity_loop_new(const Material &material, const Initial_condition &initial_condition,
//                     il::int_t i,il::int_t j, il::int_t i0, il::int_t j0,
//                     const il::Array2D<double> &kmat,Mesh mesh_total, const il::Array2D<int> &id,int p,
//                     il::Array<double> widthB, double pressure_f, il::Array<double> width_history,il::Array<double> old_history,il::int_t i_inter,il::int_t j_inter,il::io_t,il::Array<double> &delta_width,double &pressure_change,il::Array<double> &coht){// i and j stands for the element number starting from 1
//    int dof=2*(p+1);
//    il::Array2D<double> kmatnew{dof*(j-i+1),dof*(j-i+1),0.};
//    take_submatrix(kmatnew,dof*i-dof,dof*j+3-dof,dof*i-dof,dof*j+3-dof,kmat);//take_submatrix() function, the line number index is the input;
//    il::Array<double> vwc0;
//    il::Array<double> vwc{dof*(j-i+1),0.};
//    get_vwc_vc(vwc0,mesh_total,p);
//    for (il::int_t m=dof*i; m<dof*j+dof;++m){
//        vwc[m-dof*i]=vwc0[m];
//    }
//    il::Array<double> unitm{dof*(j-i+1),1.};
//    for(il::int_t s=0;s<dof*(j-i+1);s+=2){
//        unitm[s]=0.;
//    }
//    il::Array<double> delta_w_ini{widthB.size(),0.};
//    il::Array<double> coh_ini{widthB.size(),0.};
//
//    //!!!!!!!!!!!!!!!!!!!!!!!!!!! If every time we have a zero coheisive force, we could have a normal shape of the crack!!!!!!!!!!!!!!Needs to find out the reason
//    coh_ini=cohesive_force_new(material,widthB,i,j,i0,j0,old_history,p);//here we have a problem, this cohesive force should use the width-history of last time
//////when the fracture propagates, the size changes, the zeros added in vector widthB won't give a positive cohesive force, because with these zeros we always have width_edoc[s] is positive!!
//    if(i<i_inter){
//        for(il::int_t h=i;h<i_inter;h++){
//            coh_ini[dof*(h-i)+1]=0.;
//            coh_ini[dof*(h-i)+3]=0.;
//        }
//    }
//    if(j>j_inter){
//        for(il::int_t h2=j_inter;h2<j;h2++){
//            coh_ini[dof*(h2-i+1)+1]=0.;
//            coh_ini[dof*(h2-i+1)+3]=0.;
//        }
//    }
//    il::Array<double> coh{widthB.size(),0.};
//
//    int k=0;
//    double error_w=1.;
//    il::Array<il::Status> statusk{100,};
//    double positivechecker=0.;
//
//    while(k<100  && error_w> 0.00001){
//       // while((k<60 or positivechecker==1.) && error_w> 0.00001){
//        il::Array<double> width_inm=widthB;
//        il::blas(1.0,delta_w_ini,1.0,il:: io,width_inm);
//        coh=cohesive_force_new(material,width_inm,i,j,i0,j0,width_history,p);
//        il::Array<double> coh_inm=coh;
//        il::blas(-1.0,coh_ini,1.0,il::io,coh_inm);
//
//        construct_matrix(kmatnew,vwc,
//                         coh_inm,material,unitm,initial_condition,statusk[k],
//                         i,j,i0,j0,pressure_f,widthB,p,il::io,delta_width,pressure_change);
//        il::Array<double> intermediate=delta_width;
//        il::blas(-1.0,delta_w_ini,1.0,il::io,intermediate);
//        double error_inter=il::norm(intermediate,il::Norm::L2);
//        error_w=error_inter/il::norm(delta_width,il::Norm::L2);
//        il::blas(0.85,delta_width,0.15,il::io,delta_w_ini);
//        ++k;
//
////        il::Array<double> width_inf=widthB;
////        il::blas(1.0,delta_w_ini,1.0,il:: io,width_inf);
////        positivechecker=0.;
////        for(int check=0;check<width_inf.size();check+=1){
////            if(width_inf[check]<0){
////                positivechecker=1.;
////                error_w=1.;
////            }
////        }//we add this part to make sure we can always have a positive opening
////
////
////        if(k>598){
////            positivechecker=0.;
////        }
//
//    }
//    //from here you can consider to plot the current coh vector
//
//    for(il::int_t r=0;r<widthB.size();r++){
//        coht[dof*i-dof+r]=coh[r];
//    }
//
//}


    void
    propagation_loop_new(Mesh mesh_total, il::Array2D<int> &id, int p,
                         const Material &material,
                         const Initial_condition &initial_condition,
                         il::int_t i0, il::int_t j0, int nstep,
                         il::Status &status, il::io_t,
                         il::Array2C<double> &widthlist,
                         il::Array<double> &plist,
                         il::Array<double> &l_coh, il::Array<double> &l,
                         il::Array2C<double> &coh_list) {
        int n = mesh_total.nelts();
        int dof = 2 * (p + 1);
        //double l;//crack length
        il::Array2D<double> kmat{id.size(0) * id.size(1),
                                 id.size(0) * id.size(1), 0.};
        kmat=basic_assembly(mesh_total, id, p, material.Ep);
        il::Array<double> width_history{dof * n, 0.};
        for (il::int_t h = i0; h < j0 + 1; ++h) {
            width_history[dof * h - 3] = 1.;
            width_history[dof * h - 1] = 1.;
        };

        il::Array<double> old_history = width_history;

        il::Array<double> widthB = initialwidth_new(kmat, id, p, material,
                                                    initial_condition, i0, j0,
                                                    status);
        il::Array<double> width;
        il::Array2C<double> widthlistinter{nstep, dof * n,
                                           0.};//needs to be modified here
        il::Array<double> plistinter{nstep + 1, 0.};
        plistinter[0] = initial_condition.pini;
        int s = 0;
        int m = 0;
        il::int_t i = i0;
        il::int_t j = j0;
        double pressure_change = 0.;//not sure we should initialize here or not
        double pressure = initial_condition.pini;
        double time = 0;
        il::Array<double> crack_length{nstep, 0.};
        il::Array<double> length_coh{nstep, 0.};

        il::Array2C<double> coh_tt{nstep, dof * n, 0.};
        il::Array<double> coht{dof * n, 0.};
        il::int_t i_inter = i0;
        il::int_t j_inter = j0;


        while (s < nstep) {

            while (m < 10000) {//the stop limit value for m should be big enough
                if (i == 1 or j == mesh_total.nelts()) {
                    break;
                };
                il::Array<double> delta_width;
                plasticity_loop_new(material, initial_condition,
                                    i, j, i0, j0, kmat, mesh_total, id, p,
                                    widthB, pressure, width_history, il::io,
                                    delta_width, pressure_change, coht);

                il::Array<il::int_t> index = stress_criteria_new(kmat, id,
                                                                 material,
                                                                 initial_condition,
                                                                 widthB,
                                                                 delta_width, i,
                                                                 j, p);

                if (i == index[0] and j == index[1]) {
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    break;
                }
                if (index[0] <= 0 or index[1] >= n - 1) {
                    width = widthB;
                    il::blas(1.0, delta_width, 1.0, il::io, width);
                    break;
                }
                il::Array<double> widthBnew{dof * (index[1] - index[0] + 1),
                                            0.};
                for (int wbn = 1; wbn < widthB.size(); wbn += 1) {
                    widthBnew[dof * (i - index[0]) + wbn] = widthB[wbn];
                }
                il::Array<double> delta_width_new{
                        dof * (index[1] - index[0] + 1), 0.};
                for (int dwn = 0; dwn < widthB.size(); dwn += 1) {
                    delta_width_new[dof * (i - index[0]) +
                                    dwn] = delta_width[dwn];
                }
                widthB.resize(dof * (index[1] - index[0] + 1));
                widthB = widthBnew;
                width = widthBnew;


                i = index[0];
                j = index[1];


                il::blas(1., delta_width_new, 1.0, il::io, width);
                ++m;
            }

            i_inter = i;
            j_inter = j;

            il::Array<double> width_large{dof * n, 0.};
            for (int r = 0; r < width.size(); ++r) {
                width_large[dof * i + r] = width[r];
            }//needs to be carefully checked
            for (int sc = 0; sc < dof * n; ++sc) {
                widthlistinter(s, sc) = width_large[sc];
            }
            cc_length(mesh_total, material, width, i, j, width_history, p,
                      il::io, length_coh[s], crack_length[s]);

            old_history = width_history;

            width_history = write_history(width, width_history, i, material, p);

            for (int cohi = 0; cohi < dof * n; cohi++) {
                coh_tt(s, cohi) = coht[cohi];
            }


            pressure += pressure_change;
            plistinter[s + 1] = pressure;//pressure output
            widthB = width;
            time += initial_condition.timestep;//time output, or not,
            // time output can be easily deduced by nstep*timestep
            ++s;

        }
        widthlist = widthlistinter;
        plist = plistinter;
        l = crack_length;
        l_coh = length_coh;
        coh_list = coh_tt;
        status.abort_on_error();
    }


    void get_xlist(il::Array<double> &xlist, Mesh mesh_total) {
        int n = mesh_total.nelts();
        for (int i = 0; i < n; ++i) {
            xlist[2 * i] = mesh_total.coor()(i, 0);
            xlist[2 * i + 1] = mesh_total.coor()(i + 1, 0);
        }
    };
}
//}