//
// Created by DONG LIU on 2/23/17.
//

#include "Coh_Prop_new.h"
#include<il/linear_algebra/dense/blas/norm.h>

//void take_submatrixint(il::Array2D<int>& sub, int i0, int i1, int j0, int j1,
//                       const il::Array2D<int>& A){
//    IL_ASSERT((i1 - i0+1) == sub.size(0));
//    IL_ASSERT((j1-j0+1) == sub.size(1));
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



//element number i begins from 1.2.3.4...
void stress_criteria_new(il::Array<int> &index,il::Array2D<double> kmat,il::Array2D<int> id,
                     Material material,Initial_condition initial_condition,
                     il::Array<double> widthB,
                     il::Array<double> delta_w, int i, int j,int p){//here we need to get the collocation points of the next-possible failed elements not those of the activated elements

    il::Array2D<double> kl{4,4*(j-i+1),0.};
    il::Array2D<double> kr{4,4*(j-i+1),0.};
    take_submatrix(kl,4*i-4-4,4*i-4-1,4*i-4,4*j+3-4,kmat);//needs to be verified
    take_submatrix(kr,4*j,4*j+3,4*i-4,4*j+3-4,kmat);
    il::blas(1.0,delta_w,1.0,il:: io,widthB);


    double syyl1=il::dot(kl,widthB)[1];
    double syyl2=il::dot(kl,widthB)[3];
    double syyr1=il::dot(kr,widthB)[1];
    double syyr2=il::dot(kr,widthB)[3];

    if((syyl1>= (material.sigma_t+initial_condition.sigma0) or
        syyl2>=(material.sigma_t+initial_condition.sigma0)) and i>=1){
        --i;
    };
    if((syyr1>= (material.sigma_t+initial_condition.sigma0) or
        syyr2>=(material.sigma_t+initial_condition.sigma0)) and j<=id.size(0)){
        //needs to verify this part
        ++j;
    };
    if(j>=id.size(0)){
        j=int(id.size(0));
    };
    if(i<1){
        i=1;
    };
    index[0]=i;
    index[1]=j;

};//return index, the first is the new left element number i, the second is the
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

il::Array<double> width_edge_col(il::Array<double> width){
    int n=int(width.size());
    il::Array<double> width_new {n,0.};
    for(int i=0;i<n/4;++i){
        width_new[4*i+1]=(0.5+sqrt(2)/4.)*width[4*i+1]+(0.5-sqrt(2)/4.)*width[4*i+3];
        width_new[4*i+3]=(0.5-sqrt(2)/4.)*width[4*i+1]+(0.5+sqrt(2)/4.)*width[4*i+3];
    }
    return width_new;
}


il::Array<double> write_history(il::Array<double> width,il::Array<double> width_history,int i,Material material){
    int n=int(width.size());
    int n_t=int(width_history.size());
    il::Array<double> width_history_new{n_t,0.};
    il::Array<double> width_new{n,0.};
    width_new=width_edge_col(width);
    for(int m=0;m<n;++m){
        width_history_new[4*i-4+m]=width_history[4*i-4+m];
        if(width_new[m]<material.wc && width_history[4*i-4+m]<1.){
            width_history_new[4*i-4+m]=1.;
        }
    }
    return width_history_new;
}


il::Array<double> cohesive_force_new(Material material,
                                 il::Array<double> width, int i, int j, int i0, int j0,il::Array<double> width_history){
    il::Array<double> f{width.size(), 0.};
    for(int s=1; s<width.size();s+=2){
        il::Array<double> width_etoc=width_edge_col(width);
        if(0.<width_etoc[s] and width_etoc[s]<material.wc and width_history[4*i-4+s]<1.){
            f[s]=material.sigma_t;
        };
        };
        return f;
    };//return the local cohesive force, and the force in total

void cc_length(double &length_coh,double &crack_length,Mesh mesh_total,
               Material material,il::Array<double> width,int i,int j,
               il::Array<double> width_history,int p){
    int n=mesh_total.nelts();
    il::Array2D<double> xi{2,2,0.};
    il::Array2D<double> xj{2,2,0.};
    take_submatrix(xi,i-1,i,0,1,mesh_total.Coor);
    take_submatrix(xj,j-1,j,0,1,mesh_total.Coor);
    crack_length =xj(1,0)-xi(0,0);//only true to a straight linear crack
    SegmentCharacteristic segi=get_segment_DD_characteristic(xi,p);
    double m=0.;
    for(int s=1; s<width.size();s+=2){
        il::Array<double> width_etoc=width_edge_col(width);
        if(0.<width_etoc[s] and width_etoc[s]<material.wc and width_history[4*i-4+s]<1.){
            ++m;
        };
    };
    length_coh =m/2*segi.size;//only true to uniform mesh

}



void initialwidth_new(il::Array<double> &delta_w_ini, il::Array2D<double> kmat,
                  il::Array2D<int> id,int p,
                  Material material,Initial_condition initial_condition,
                  int i,int j,il::Status &status){
    il::Array2D<double> kini{4*(j-i+1),4*(j-i+1),0.};
    take_submatrix(kini,4*i-4,4*j+3-4,4*i-4,4*j+3-4,kmat);
    il::Array<double> F{4*j-4*i+4,
                        -initial_condition.pini+initial_condition.sigma0};
    for(int m=0;m<4*j-4*i+4;m+=2){
        F[m]=0;
    };
    delta_w_ini=il::linear_solve(kini,F,il::io,status);

};


void plasticity_loop_new(il::Array<double> &delta_width,double &pressure_change,
                     Material material, Initial_condition initial_condition,
                     int i,int j, int i0, int j0,
                     il::Array2D<double> kmat,Mesh mesh_total, il::Array2D<int> id,int p,
                     il::Array<double> widthB, double pressure_f, il::Array<double> width_history){// i and j stands for the element number starting from 1
    il::Array2D<double> kmatnew{4*(j-i+1),4*(j-i+1),0.};
    take_submatrix(kmatnew,4*i-4,4*j+3-4,4*i-4,4*j+3-4,kmat);//take_submatrix() function, the line number index is the input;
    il::Array<double> vwc0;
    il::Array<double> vwc{4*(j-i+1),0.};
    get_vwc_vc(vwc0,mesh_total,p);
    for (int m=4*i; m<4*j+4;++m){
        vwc[m-4*i]=vwc0[m];
    }
    il::Array<double> unitm{4*(j-i+1),1.};
    for(int s=0;s<4*(j-i+1);s+=2){
        unitm[s]=0.;
    }
    il::Array<double> delta_w_ini{widthB.size(),0.};
    il::Array<double> coh_ini{widthB.size(),0.};
    coh_ini=cohesive_force_new(material,widthB,i,j,i0,j0,width_history);

    int k=0;
    double error_w=1.;
    il::Array<il::Status> statusk{60,};

    while(k<60 && error_w>pow(10.,-5)){

        il::Array<double> width_inm=widthB;
        il::blas(1.0,delta_w_ini,1.0,il:: io,width_inm);
        il::Array<double> coh=cohesive_force_new(material,width_inm,i,j,i0,j0,width_history);
        il::Array<double> coh_inm=coh;
        il::blas(-1.0,coh_ini,1.0,il::io,coh_inm);

        construct_matrix(delta_width,pressure_change,kmatnew,vwc,
                         coh_inm,material,unitm,initial_condition,statusk[k],
                         i,j,i0,j0,pressure_f,widthB);
        il::Array<double> intermediate=delta_width;
        il::blas(-1.0,delta_w_ini,1.0,il::io,intermediate);
        double error_inter=il::norm(intermediate,il::Norm::L2);
        error_w=error_inter/il::norm(delta_width,il::Norm::L2);
        il::blas(0.85,delta_width,0.15,il::io,delta_w_ini);
        ++k;

    }
}


void propagation_loop_new(il::Array2D<double> &widthlist,il::Array<double> &plist,
                          il::Array<double> &l_coh,il::Array<double> &l,
                      Mesh mesh_total,il::Array2D<int> id,int p,
                      Material material, Initial_condition initial_condition,
                      int i0, int j0, int nstep,il::Status status){
    int n=mesh_total.nelts();
    //double l;//crack length
    il::Array2D<double> kmat{id.size(0)*id.size(1),id.size(0)*id.size(1),0.};
    BasicAssembly(kmat,mesh_total,id,p,material.Ep);
    il::Array<double> width_history{4*n,0.};
    for(int h=i0;h<j0+1;++h){
        width_history[4*h-3]=1.;
        width_history[4*h-1]=1.;
    };

    il::Array<double> widthB;
    initialwidth_new(widthB,kmat,id,p,material,initial_condition,i0,j0,status);
    il::Array<double> width;
    il::Array2D<double> widthlistinter{nstep,4*n,0.};//needs to be modified here
    il::Array<double> plistinter{nstep+1,0.};
    plistinter[0]=initial_condition.pini;
    int s=0;
    int m=0;
    int i=i0;
    int j=j0;
    double pressure_change;
    double pressure=initial_condition.pini;
    double time=0;
    il::Array<double> crack_length{nstep,0.};
    il::Array<double> length_coh{nstep,0.};


    while(s<nstep){

        while(m<10000){
            if(i==1 or j==mesh_total.nelts()){
                break;
            };
            il::Array<double> delta_width;
            plasticity_loop_new(delta_width,pressure_change,material,initial_condition,
                            i,j,i0,j0,kmat,mesh_total,id,p,widthB,pressure,width_history);

            il::Array<int> index{2,0};
            stress_criteria_new(index,kmat,id,material,initial_condition,
                            widthB,delta_width,i,j,p);

            if(i==index[0] and j==index[1]){
                width=widthB;
                il::blas(1.0,delta_width,1.0,il::io,width);
                break;
            }
            if(index[0]<=0 or index[1]>=n-1){
                width=widthB;
                il::blas(1.0,delta_width,1.0,il::io,width);
                break;
            }
            il::Array<double> widthBnew{4*(index[1]-index[0]+1),0.};
            for(int wbn=1;wbn<widthB.size();wbn+=1){
                widthBnew[4*(i-index[0])+wbn]=widthB[wbn];
            }
            il::Array<double> delta_width_new{4*(index[1]-index[0]+1),0.};
            for(int dwn=0;dwn<widthB.size();dwn+=1){
                delta_width_new[4*(i-index[0])+dwn]=delta_width[dwn];
            }
            widthB.resize(4*(index[1]-index[0]+1));
            widthB=widthBnew;
            width=widthBnew;
            i=index[0];
            j=index[1];


            il::blas(1.,delta_width_new,1.0,il::io,width);
            ++m;
        }

        il::Array<double> width_large{4*n,0.};
        for(int r=0;r<width.size();++r){
            width_large[4*i+r]=width[r];
        }//needs to be carefully checked
        for(int sc=0;sc<4*n;++sc){
            widthlistinter(s,sc)=width_large[sc];
        }
        cc_length(length_coh[s],crack_length[s],mesh_total,material,width,i,j,width_history,p);

        width_history=write_history(width,width_history,i,material);



        pressure+=pressure_change;
        plistinter[s+1]=pressure;//pressure output
        widthB=width;
        time+=initial_condition.timestep;//time output, or not,
        // time output can be easily deduced by nstep*timestep
        ++s;

    }
    widthlist=widthlistinter;
    plist=plistinter;
    l=crack_length;
    l_coh=length_coh;
    status.abort_on_error();
}

