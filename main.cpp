
#include <iostream>
#include <string>
#include <cmath>
#include <complex>

#include <il/Array.h>
#include <il/math.h>

#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>


#include "Elasticity2D.h"
#include "Mesh.h"


void dofhandle_DG2D(il::Array2D<int>& dofhandle,Mesh mesh,int p)
{
// function creating a matrix of dof handle - for a piece-wise linear variation per element (Discontinous Galerkin type) on a 1d Mesh object for the case of 2 Degrees of Freedoms per node
  int ne=mesh.nelts();
  int ndof = ne*p*2*2;

//  il::Array2D<int> dofhandle{ne,2*p+2,0} ;

  int j ;

  for (int i = 0; i < ne; ++i) {
    j = i*(2*p+2) ;
    for (int k=0; k<2*p+2;++k) {
      dofhandle(i, k) = j + k;
    }
  }
//  return dofhandle; - starts at 0 for dof c++ style!
}


il::Array2D<double> rotation_matrix_2D(double theta)
{
  il::Array2D<double> R{2,2,0.};

  R(0,0)=cos(1.*theta);
  R(0,1)=-1.*sin(1.*theta);
  R(1,0)=sin(theta);
  R(1,1)=cos(theta);

  return R;
}

//void rotation_matrix_2D(il::Array2D<double>& R, double theta)
//{
//  R(0,0)=cos(theta);
//  R(0,1)=-1.*sin(theta);
//  R(1,0)=sin(theta);
//  R(0,1)=cos(theta);
//
//}

//----------------------------------------------------
// Segment Characteristics....
// use a struct ? (may be oject is better ? )
struct SegmentCharacteristic{
  double size;
  double theta; // angle w.r. to e_1
  il::Array2D<double> CollocationPoints;
  il::Array<double> n{2};
  il::Array<double> s{2};
  il::Array<double> Xmid{2}; //segment mid points coordinates.
};


SegmentCharacteristic get_segment_characteristic(const il::Array2D<double> Xs , int p){
  IL_ASSERT(Xs.size(0)==2);
  IL_ASSERT(Xs.size(1)==2);
  SegmentCharacteristic segment;

// compute element size
  il::Array<double> xdiff{2},s{2},n{2},xmean{2};
  il::Array2D<double> Xcol{p+1,2,0},R{2,2,0.};
  il::Array<double> xaux{2};

  xdiff[0]=Xs(1,0)-Xs(0,0);
  xdiff[1]=Xs(1,1)-Xs(0,1);

  segment.size=sqrt(pow(xdiff[0],2)+pow(xdiff[1],2));

  //s=xdiff; // tangent vector
  s[0]=xdiff[0]/segment.size;
  s[1]=xdiff[1]/segment.size;
  n[0]=-1.*s[1];n[1]=s[0]; //normal vector

  segment.s=s;
  segment.n =n;

  segment.theta = acos(s[0]/sqrt(pow(s[0],2)+pow(s[1],2)));
  if(s[1]<0) {segment.theta=-segment.theta;};

  xmean[0]=(Xs(1,0)+Xs(0,0))/2.;  xmean[1]=(Xs(1,1)+Xs(0,1))/2.;
  segment.Xmid=xmean;

  switch(p) {
    case 1 : {
      Xcol(0,0)=-1./sqrt(2.); Xcol(0,1)=0.;
      Xcol(1,0)=1./sqrt(2.);  Xcol(1,1)=0.;
    }; break;

    case 0 :{
    Xcol(0,0)=0.;Xcol(0,1)=0;
    };  break;
    default : std::cout << "error\n"; //  error
      break ;
    };

  R=rotation_matrix_2D(segment.theta);

  for (int i=0; i <p+1;++i){

    xaux[0] =(segment.size)*Xcol(i,0)/2.;
    xaux[1] =(segment.size)*Xcol(i,1)/2.;

    xaux=il::dot(R,xaux);

    Xcol(i,0)=xaux[0]+xmean[0];
    Xcol(i,1)=xaux[1]+xmean[1];

  }

  segment.CollocationPoints = Xcol;

  return segment; // return structure with all we need on the segment.
}
//----------------------------------------------------


void take_submatrix(il::Array2D<double>& sub, int i0, int i1, int j0, int j1, const il::Array2D<double>& A){
  IL_ASSERT((i1 - i0+1) == sub.size(0));
  IL_ASSERT((j1-j0+1) == sub.size(1));

  for(int i = i0; i <= i1;++i) {
    for (int j=j0; j<= j1;++j){
      sub(i-i0,j-j0)=A(i,j);

    }

  }
}

void set_submatrix(il::Array2D<double>& A, int i0, int i1, const il::Array2D<double>& B) {
  IL_ASSERT(i0 + B.size(0) <= A.size(0));
  IL_ASSERT(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
} // e.g. set_submatrix(A, 2, 3, B);


////////////////////////////////////////////////////////////////////////////////
int main() {

  int n=42, p=1;

  double h=2./(n-1) ; //  element size

  il::Array<double> x{n};

  il::Array2D<double> xy{n, 2, 0.0};
  il::Array2D<int> myconn{n-1, 2, 0.0};
  il::Array2D<int> id{n-1,4,0};

  int ndof=(n-1)*4;
  double Ep=1.; //Plane strain Young's modulus

  //  std::complex(double re = 0.0, double im = 0.0) myC2;
  //  myC.real(2.);
  // myC.imag(1.);

  // Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  // Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // 1D mesh ....
  for (int i = 0; i < xy.size(0); ++i) {
      xy(i,0)=-1.+i*h;
      xy(i,1)=0.;
   }
  for (int  i=0; i< myconn.size(0); ++i)
  {
      myconn(i,0)=i;
      myconn(i,1)=i+1;
  }

  Mesh mesh ;
  mesh.set_values(xy,myconn);

  dofhandle_DG2D(id,mesh,p); // dof handle for DDs

  il::Array2D<double> xe{2,2,0},xec{2,2,0};

  SegmentCharacteristic mysege,mysegc;

  il::Array2D<double> K{ndof,ndof,0.},R{2,2,0.};
  il::Array<int> dofe{2*(p+1),0},dofc{2*(p+1),0};
  il::Array<double> sec{2},nec{2},xcol{2};
  il::Array2D<double>  stnl{2,4,0.};
  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout << myconn.size(0)<< "\n";

  // loop on collocation points....
//  for (int e=0;e<mesh.nelts();++e) {
//    take_submatrix(xe,mesh.conn(e,0),mesh.conn(e,1),0,1,mesh.Coor); // take the coordinates of element e from the mesh object
//    mysege=get_segment_characteristic(xe,p); // get the segment characteristic.
//
//    std::cout << " \n";
//    for (int ic=0;ic<p+1;++ic) { // loop on collocation points
//      std::cout << " collocation pts elt e \n";
//      std::cout << mysege.CollocationPoints(ic, 0) << "  " << mysege.CollocationPoints(ic, 1) << "\n";
//    }
//    }

   std::cout << "------\n";
// double loop on elements to create the stiffness matrix...
  for (int e=0;e<mesh.nelts();++e){ // loop on all  elements

    take_submatrix(xe,mesh.conn(e,0),mesh.conn(e,1),0,1,mesh.Coor); // take the coordinates of element e from the mesh object
    mysege=get_segment_characteristic(xe,p); // get the segment characteristic.
    R=rotation_matrix_2D(mysege.theta);      //Rotation matrix of the element w.r. to x-axis.
    for (int i=0;i<2*(p+1);++i){ // vector of dof id of the element e
      dofe[i]=id(e,i);
    };

    for (int j=0;j<mesh.nelts();++j){// loop on all  elements

      take_submatrix(xec,mesh.conn(j,0),mesh.conn(j,1),0,1,mesh.Coor); // takes the coordinates of element j
      mysegc=get_segment_characteristic(xec,p);

      sec=il::dot(R,mysegc.s); // tangent of elt j
      nec=il::dot(R,mysegc.n); // normal of elt j

      for (int i=0;i<2*(p+1);++i){
        dofc[i]=id(j,i); // vector of dof id of the  element j
       };

      for (int ic=0;ic<p+1;++ic){ // loop on collocation points
         // we switch to the frame of element e
        for (int i=0;i<2;++i){
          xcol[i]=mysegc.CollocationPoints(ic,i)-mysege.Xmid[i]; //
        }
        xcol=il::dot(R,xcol);

        NormalShearStressKernel_LinearDD(stnl,xcol,mysege.size,sec,nec,Ep);
        //
        set_submatrix(K,dofc[2*ic],dofe[0],stnl);
      }
    }
  }
//  mysege=get_segment_characteristic(xs,p);

//  il::Array2D<double>  stnl2{2,4,0.};
//  il::Array<double> xcc{2};
//  xcc[0]=2;xcc[1]=0.;
//
//  std::cout <<" seg size " << mysege.size << " \n";

//  NormalShearStressKernel_LinearDD(stnl2,xcc,mysege.size,mysege.s,mysege.n,1.);
////
//  std::cout << stnl2(0,1)<< "\n";
//  std::cout << stnl2(0,1)<< "\n";
//  std::cout << stnl2(0,1)<< "\n";

//  il::Array2D<double> A{n, n, 1.0};
//  il::Array2D<double> B{n, 1, 1.0};
//  il::Array2D<double> Stress{2,4, 0.0};

//  double myh=0.3;
//  double xp =3.;
//  double yp = 5.;

//  Stress =StressesKernelLinearDD(myh,1.,xp,yp);
//  std::cout << Stress(0,0) ; std::cout <<  " ...\n";
//  std::cout << Stress(0,1) ; std::cout <<  " ...\n";
//  std::cout << Stress(0,2) ; std::cout <<  " ...\n";
//  std::cout << Stress(0,3) ; std::cout <<  " ...\n";
//  std::cout <<  " ...\n";
//
//  std::cout << Stress(1,0) ; std::cout <<  " ...\n";
//  std::cout << Stress(1,1) ; std::cout <<  " ...\n";
//  std::cout << Stress(1,2) ; std::cout <<  " ...\n";
//  std::cout << Stress(1,3) ; std::cout <<  " ...\n";

//  C = il::dot(A, B);
//  v2 = il::dot(A, v);
//
// std::cout << C ;


  std::string mystring;
  mystring = "This is the end of this code... ";
  std::cout << " mat entries \n" ;
  for (int i=0;i<8; ++i){
    std::cout << K(i,0) << " " << K(i,1) << " " << K(i,2) << " " << K(i,3) << " " << K(i,4) << " " << K(i,5) << " " << K(i,6)<< " " << K(i,7) ;
    std::cout <<  " ...\n";

  }

  return 0;

}



