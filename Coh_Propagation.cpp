//
// Created by DONG LIU on 1/30/17.
//

#include "Coh_Propagation.h"
//#include <il/linear_algebra/dense/blas/blas.h>
// if you add this, there will be an error in building the project
#include <il/linear_algebra.h>
#include<il/linear_algebra/dense/blas/norm.h>
//#include <il/linear_algebra/dense/factorization/LU.h>
// if you add this, there will be an error in building the project
#include "math.h"
#include <il/core/Status.h>



void take_submatrixint(il::Array2D<int>& sub, int i0, int i1, int j0, int j1,
                    const il::Array2D<int>& A){
    IL_ASSERT((i1 - i0+1) == sub.size(0));
    IL_ASSERT((j1-j0+1) == sub.size(1));

    for(int i = i0; i <= i1;++i) {
        for (int j=j0; j<= j1;++j){
            sub(i-i0,j-j0)=A(i,j);
        }

    }
}

Mesh take_submesh(int n1, int n2, Mesh mesh_total){//we consider n1 and n2
// represents for the element number
    Mesh meshC;
    il::Array2D<int> connC{n2-n1+1,2,0};
    take_submatrixint(connC,0,n2-n1,0,1,mesh_total.conn);//here also, the connectivity begins from 0.1.2.3...
    il::Array2D<double> CoorC{n2-n1+2,2,0.};
    take_submatrix(CoorC,n1-1,n2,0,1,mesh_total.Coor);
    meshC.set_values(CoorC,connC);
    return meshC;
};

void material_condition(Material &material, Initial_condition &initial_condition,
                        double wc1, double sigma_t1,double Ep1,double U1,
                        double pini1,double Q01,double timestep1,double sigma01){
    material.Ep=Ep1;
    material.sigma_t=sigma_t1;
    material.U=U1;
    material.wc=wc1;
    initial_condition.pini=pini1;
    initial_condition.Q0=Q01;
    initial_condition.sigma0=sigma01;
    initial_condition.timestep=timestep1;
}

void stress_criteria(il::Array<int> &index,Mesh mesh_total,il::Array2D<int> id,
                     Material material,Initial_condition initial_condition,
                     il::Array<double> widthB,
                     il::Array<double> delta_w, int i, int j,int p){//here we need to get the collocation points of the next-possible failed elements not those of the activated elements
    Mesh meshC =take_submesh(i,j,mesh_total);// i eme line and j eme line
    Mesh mesh_large=take_submesh(i-1,j+1,mesh_total);
    il::blas(1.0,delta_w,1.0,il:: io,widthB);
    //inline void blas(double alpha, const il::Array<double>& x, double beta,
    //  il::io_t, il::Array<double>& y)
    // x, y are vectors
    // y <- alpha x + beta y
    il::Array2D<double> colocationlseg{2,2,0.},colocationrseg{2,2,0.};
    take_submatrix(colocationlseg,0,1,0,1,mesh_large.Coor);//the first line
    take_submatrix(colocationrseg,mesh_large.nelts()-1,mesh_large.nelts(),0,1,mesh_large.Coor);
    //the last line

    SegmentCharacteristic segl=get_segment_DD_characteristic(colocationlseg,p);
    SegmentCharacteristic segr=get_segment_DD_characteristic(colocationrseg,p);

    il::Array2D<double> colocationl{2,2,0.},colocationr{2,2,0.};
    colocationl=segl.CollocationPoints;
    colocationr=segr.CollocationPoints;

    //ajust the id matrix to the current mesh
    il::Array2D<int> idC{j-i+1,4,0};
    take_submatrixint(idC,0,j-i,0,3,id);//here you need to check,because the number of degree of freedom always begins at 0.1.2.3...

    double syyl1=stressoutput(meshC,idC,p,material.Ep,colocationl(0,0),
                              colocationl(0,1),widthB)[1];
    double syyl2=stressoutput(meshC,idC,p,material.Ep,colocationl(1,0),
                              colocationl(1,1),widthB)[1];
    double syyr1=stressoutput(meshC,idC,p,material.Ep,colocationr(0,0),
                              colocationr(0,1),widthB)[1];
    double syyr2=stressoutput(meshC,idC,p,material.Ep,colocationr(1,0),
                              colocationr(1,1),widthB)[1];

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

il::Array<double> former_pressure(double pressure_f,int i,int j,int i0,int j0){
    il::Array<double> fP{4*(j-i+1),0.};
    if(i<i0){
        for(int m=1;m<4*(i0-i);m+=2){
            fP[m]=-pressure_f;
        }
    }
    if(j>j0){
        for(int s=4*(j0-i+1)+1;s<4*(j-i+1);s+=2){
            fP[s]=-pressure_f;
        }
    }
    return fP;
}

il::Array<double> former_k(il::Array2D<double> kmatC,int i,int j,
                           int i0,int j0,il::Array<double> widthB){
    il::Array<double> fK{4*(j-i+1),0.};
    il::Array<double> vectem=il::dot(kmatC,widthB);
    if(i<i0){
        for (int k = 0; k <4*(i0-i) ; ++k) {
            fK[k]=-vectem[k];
        }
    }
    if(j>j0){
        for (int s = 4*(j0-i+1); s <4*(j-i+1) ; ++s) {
            fK[s]=-vectem[s];
        }
    }
    return fK;
}


void construct_matrix(il::Array<double> &width, double &pressure,
                      il::Array2D<double> kmatC,il::Array<double> vwc,
                      il::Array<double> cohf,Material material,
                      il::Array<double> unit,Initial_condition initial_condition,
                      il::Status &status,int i,int j,int i0,int j0,
                      double pressure_f,il::Array<double> widthB){
    int n =int(kmatC.size(0));
    il::Array2D<double> newmatrix{n+1,n+1,0};
    il::Array<double> newvector{n+1,0};
    for(int m=0;m<n;++m){
        for(int s=0;s<n;++s)
        newmatrix(m,s)=kmatC(m,s);
    };
    newmatrix(n,n)=material.U;
    il::Array<double> fP=former_pressure(pressure_f,i,j,i0,j0);
    il::Array <double> fK=former_k(kmatC,i,j,i0,j0,widthB);//added the 6th Februray


    for(int q=0;q<n;++q){
        newmatrix(q,n)=unit[q];
        newmatrix(n,q)=vwc[q];
        newvector[q]=cohf[q]+fP[q]+fK[q];
    };
    newvector[n]=initial_condition.Q0 * initial_condition.timestep;
    il::Array<double> widthinter{n,0.};
    il::Array<double> dd = il::linear_solve(newmatrix,newvector, il::io, status);
    for(int t=0;t<n;++t){
        widthinter[t]=dd[t];
    };
    width=widthinter;
    pressure=dd[n];
    status.abort_on_error();//every time we execuate the constructmatrix,
    // we will get the status_
    // status.abort_on_error();//we can add this line to do the same thing,
    // so actually you don't have to change the "Status.h"
}

void get_vwc_vc(il::Array<double> &vwc,Mesh mesh_total,int p){

    int n=mesh_total.nelts();
    il::Array<double> vwcinter{4*n,0.};
    for (int i=0;i<n;++i){
        il::Array2D<double> xi{2,2,0.};
        take_submatrix(xi,i,i+1,0,1,mesh_total.Coor);
        SegmentCharacteristic segi=get_segment_DD_characteristic(xi,p);
        vwcinter[4*i]=0.;
        vwcinter[4*i+1]=0.5*segi.size;
        vwcinter[4*i+2]=0.;
        vwcinter[4*i+3]=0.5*segi.size;
    };
    vwc=vwcinter;
}

il::Array<double> cohesive_force(Material material,
                    il::Array<double> width, int i, int j, int i0, int j0){
    il::Array<double> f{width.size(), 0.};
    if(i==i0 && j==j0){
        return f;
    }
    else{
        for(int s=1; s<width.size();s+=2){
            if(0.<width[s] and width[s]<material.wc){
                f[s]=material.sigma_t;
            };
        };
        return f;
    }


};//return the local cohesive force, and the force in total

void initialwidth(il::Array<double> &delta_w_ini, Mesh mesh_total,
                  il::Array2D<int> id,int p,
                  Material material,Initial_condition initial_condition,
                  int i,int j,il::Status &status){
    il::Array2D<double> kmat{4*mesh_total.nelts(), 4*mesh_total.nelts(), 0.};
    BasicAssembly(kmat,mesh_total,id,p,material.Ep);
    il::Array2D<double> kini{4*(j-i+1),4*(j-i+1),0.};
    take_submatrix(kini,4*i,4*j+3,4*i,4*j+3,kmat);
    il::Array<double> F{4*j-4*i+4,
                        -initial_condition.pini+initial_condition.sigma0};
    for(int m=0;m<4*j-4*i+4;m+=2){
        F[m]=0;
    };
    delta_w_ini=il::linear_solve(kini,F,il::io,status);

};


void plasticity_loop(il::Array<double> &delta_width,double &pressure_change,
                     Material material, Initial_condition initial_condition,
                     int i,int j, int i0, int j0,
                     Mesh mesh_total, il::Array2D<int> id,int p,
                     il::Array<double> widthB, double pressure_f){// i and j stands for the element number starting from 1
    il::Array2D<double> kmat{id.size(0)*id.size(1),id.size(0)*id.size(1),0.};
    il::Array2D<double> kmatnew{4*(j-i+1),4*(j-i+1),0.};
    BasicAssembly(kmat,mesh_total,id,p,material.Ep);
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
    if(i<i0 or j>j0){
        coh_ini=cohesive_force(material,widthB,i,j,i0,j0);
    }//added the 6th February to make sure no cohesive force exists in the already opening area
    //initialwidth(delta_w_ini,mesh_total,id,p,material,initial_condition,i,j);
    int k=0;
    double error_w=1.;
    il::Array<il::Status> statusk{60,};

    while(k<60 && error_w>pow(10.,-5)){

        il::Array<double> width_inm=widthB;
        il::blas(1.0,delta_w_ini,1.0,il:: io,width_inm);
        il::Array<double> coh=cohesive_force(material,width_inm,i,j,i0,j0);
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
//        if(error_w <= pow(10.,-5)){
//            for(int m=0;m<k;++m){
//                statusk[m].set_false();
//            }
//       };
    }
}


void propagation_loop(il::Array2D<double> &widthlist,il::Array<double> &plist,
                      Mesh mesh_total,il::Array2D<int> id,int p,
                      Material material, Initial_condition initial_condition,
                      int i0, int j0, int nstep,il::Status status){
    int n=mesh_total.nelts();
    //double l;//crack length

    il::Array<double> widthB;
    initialwidth(widthB,mesh_total,id,p,material,initial_condition,i0,j0,status);
    il::Array<double> width;
    il::Array2D<double> widthlistinter{nstep,4*n,0.};//needs to be modified here
    il::Array<double> plistinter{nstep,0.};
    int s=0;
    int m=0;
    int i=i0;
    int j=j0;
    double pressure_change;
    double pressure=initial_condition.pini;
    double time=0;
    while(s<nstep){

        while(m<160){
            if(i==1 or j==mesh_total.nelts()){
                break;
            };
            il::Array<double> delta_width;
            plasticity_loop(delta_width,pressure_change,material,initial_condition,
                            i,j,i0,j0,mesh_total,id,p,widthB,pressure);

            il::Array<int> index{2,0};
            stress_criteria(index,mesh_total,id,material,initial_condition,
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
        pressure+=pressure_change;
        plistinter[s]=pressure;//pressure output
        widthB=width;
        time+=initial_condition.timestep;//time output, or not,
        // time output can be easily deduced by nstep*timestep
        ++s;

    }
    widthlist=widthlistinter;
    plist=plistinter;
    status.abort_on_error();
}



void get_xlist(il::Array<double> &xlist,Mesh mesh_total){
    int n=mesh_total.nelts();
    for(int i=0;i<n;++i){
        xlist[2*i]=mesh_total.Coor(i,0);
        xlist[2*i+1]=mesh_total.Coor(i+1,0);
    }
};
