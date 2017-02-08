//
// Created by DONG LIU on 2/7/17.
//

#include "Viscosity.h"
#include <cmath>


//all functions in this cpp needs to consider the change of the size of matrix, that's why for now we don't use the index i0,j0 yet


void matrix_lw(il::Array<double> widthB,  double visco,
               Mesh mesh_total,int p,int i, int j, int i0, int j0){//widthB should be the ajusted-width profile at the previous moment
//the size is (N+1)*(N+1)

    Mesh meshC=take_submesh(i,j,mesh_total);
    int n=meshC.nelts();
    il::Array2D<double>m_lw{n+1,n+1,0.};
    for(int s=0;s<n+1;++s){//we set the structure like this, because it's tri-band
        //each line has three elements
        if(s==0){
            il::Array2D<double> xsm0{2,2,0.};
            take_submatrix(xsm0,0,1,0,1,meshC.Coor);
            SegmentCharacteristic segs0=get_segment_DD_characteristic(xsm0,p);
            double h0=segs0.size;
            m_lw(0,0)=pow(widthB[1],3.)/12./h0/visco;
            m_lw(0,1)=-pow((widthB[1]+widthB[3])/2.,3.)/12./h0/visco;
        }
        if(s==n){
            il::Array2D<double> xsmn{2,2,0.};
            take_submatrix(xsmn,n-1,n,0,1,meshC.Coor);
            SegmentCharacteristic segsn=get_segment_DD_characteristic(xsmn,p);
            double hn=segsn.size;
            m_lw(n,n-1)=-pow((widthB[4*n-3]+widthB[4*n-1])/2.,3.)/12./hn/visco;
            m_lw(n,n)=pow((widthB[4*n-3]+widthB[4*n-1])/2.,3.)/12./hn/visco;
        }
        else{
            double w_average1=(widthB[4*s-3]+widthB[4*s-1])/2.;//there's a problem, it should correspond to the fact that each element has two widths values
            il::Array2D<double> xsm1{2,2,0.};
            take_submatrix(xsm1,s-1,s,0,1,meshC.Coor);
            SegmentCharacteristic segs1=get_segment_DD_characteristic(xsm1,p);
            double h1=segs1.size;
            double w_average2=(widthB[4*s+1]+widthB[4*s+3])/2.;//here too, not correct. could be linked to the degree of freedom's matrix
            il::Array2D<double> xsm2{2,2,0.};
            take_submatrix(xsm2,s,s+1,0,1,meshC.Coor);
            SegmentCharacteristic segs2=get_segment_DD_characteristic(xsm2,p);
            double h2=segs2.size;
            m_lw(s,s-1)=-pow(w_average1,3.)/12./h1/visco;
            m_lw(s,s)=pow(w_average1,3.)/12./h1/visco+pow(w_average2,3.)/12./h2/visco;
            m_lw(s,s+1)=-pow(w_average2,3.)/12./h2/visco;
        }
    }
}

void matrix_edge_col(Mesh mesh_total,int i, int j, int i0, int j0,int p){//from edge to collocation points, the size is 4N*(N+1), because the width is 4*N(corresponding to the degrees of freedom)
//we base on the number of the elements
//presB:size N+1
    Mesh meshC=take_submesh(i,j,mesh_total);
    int n=meshC.nelts();
    il::Array2D<double> edge{4*n,n+1,0.};

    for(int s=0;s<n;++s){//s stands for the element number-1
        il::Array2D<double> xs{2,2,0.};
        take_submatrix(xs,s,s+1,0,1,meshC.Coor);
        SegmentCharacteristic segs=get_segment_DD_characteristic(xs,p);
        edge(4*s+1,s)=(xs(1,0)-segs.CollocationPoints(0,0))/segs.size;
        edge(4*s+1,s+1)=(segs.CollocationPoints(0,0)-xs(0,0))/segs.size;
        edge(4*s+3,s)=(xs(1,0)-segs.CollocationPoints(1,0))/segs.size;
        edge(4*s+3,s+1)=(segs.CollocationPoints(1,0)-xs(0,0))/segs.size; //we only care about the x coordinate
    }
}


void matrix_vp(il::Array<double> widthB,Mesh mesh_total,double beta,int i,int j, int i0, int j0,int p){//the size is (N+1)*(N+1), and it's a function of width profile,
// beta here is the compressibility coefficient, we consider it as a constant
    Mesh meshC=take_submesh(i,j,mesh_total);
    int n=meshC.nelts();
    il::Array2D<double> m_vp{n+1,n+1,0.};
    for(int s=1;s<n;++s){//s stands for the element number-1 here
        il::Array2D<double> xs1{2,2,0.};
        take_submatrix(xs1,s-1,s,0,1,meshC.Coor);
        il::Array2D<double> xs2{2,2,0.};
        take_submatrix(xs2,s,s+1,0,1,meshC.Coor);
        SegmentCharacteristic segs1=get_segment_DD_characteristic(xs1,p);
        SegmentCharacteristic segs2=get_segment_DD_characteristic(xs2,p);
        m_vp(s,s-1)=1./12./(segs1.size)*(widthB[4*s-1]+1./2.*widthB[4*s-3]);
        m_vp(s,s)=1./12./(segs1.size)*(7./2.*widthB[4*s-1]+widthB[4*s-3])
                  +1./12./(segs2.size)*(7./2.*widthB[4*s+1]+widthB[4*s+3]);
        m_vp(s,s+1)=1./12./(segs2.size)*(widthB[4*s+1]+1./2.*widthB[4*s+3]);
    }
    il::Array2D<double> xs0{2,2,0.};
    take_submatrix(xs0,0,1,0,1,meshC.Coor);
    SegmentCharacteristic seg0=get_segment_DD_characteristic(xs0,p);
    m_vp(0,0)=1./12./(seg0.size)*(7./2.*widthB[1]+widthB[3]);
    m_vp(0,1)=1./12./(seg0.size)*(widthB[1]+1./2.*widthB[3]);
    il::Array2D<double> xsn{2,2,0.};
    take_submatrix(xsn,n-1,n,0,1,meshC.Coor);
    SegmentCharacteristic segn=get_segment_DD_characteristic(xsn,p);
    m_vp(n,n-1)=1./12./(segn.size)*(widthB[4*n-1]+1./2.*widthB[4*n-3]);
    m_vp(n,n)=1./12./(segn.size)*(7./2.*widthB[4*n-1]+widthB[4*n-3]);

}

void matrix_vw(Mesh mesh_total,int i,int j, int i0, int j0,int p){//here in order to simplify the problem, we consider alpha here is a constant and alpha=1.
//and the size is (N+1)*4N
//we can build a whole matrix for the total mesh, and then for constructing the final matrix, we just take a sub-matrix.
    Mesh meshC=take_submesh(i,j,mesh_total);
    int n=meshC.nelts();
    il::Array2D<double> m_vw{n+1,4*n,0.};
    for(int s=1;s<n;++s){//s stands for the element number-1 here
        il::Array2D<double> xs1{2,2,0.};
        take_submatrix(xs1,s-1,s,0,1,meshC.Coor);
        il::Array2D<double> xs2{2,2,0.};
        take_submatrix(xs2,s,s+1,0,1,meshC.Coor);
        SegmentCharacteristic segs1=get_segment_DD_characteristic(xs1,p);
        SegmentCharacteristic segs2=get_segment_DD_characteristic(xs2,p);
        m_vw(s,4*s-3)=1./12./(segs1.size)*3./2.;
        m_vw(s,4*s-1)=1./12./(segs1.size)*9./2.;
        m_vw(s,4*s+1)=1./12./(segs2.size)*9./2.;
        m_vw(s,4*s+3)=1./12./(segs2.size)*3./2.;
    }
    il::Array2D<double> xs0{2,2,0.};
    take_submatrix(xs0,0,1,0,1,meshC.Coor);
    SegmentCharacteristic seg0=get_segment_DD_characteristic(xs0,p);
    m_vw(0,1)=1./12./(seg0.size)*9./2.;
    m_vw(0,3)=1./12./(seg0.size)*3./2.;
    il::Array2D<double> xsn{2,2,0.};
    take_submatrix(xsn,n-1,n,0,1,meshC.Coor);
    SegmentCharacteristic segn=get_segment_DD_characteristic(xsn,p);
    m_vw(n,4*n-3)=1./12./(segn.size)*3./2.;
    m_vw(n,4*n-1)=1./12./(segn.size)*9./2.;
}

void matrix_source(Mesh mesh_total,int i0,int j0,int i,int j,Initial_condition initial_condition){//Vector the size is N+1
//we consider i0 and j0 represent the element numbers and i0 and j0 begin from 1.2.3.4...
    int n=j-i+1;
    il::Array<double> source{n+1,0.};
    if(j0-i0==1){
        source[i0]=initial_condition.Q0;
    }
    else{
        int s=(j0-i0)/2;
        source[i0-1+s]=initial_condition.Q0;//needs to be improved here
    }
}