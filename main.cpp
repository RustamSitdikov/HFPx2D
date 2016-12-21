
#include <iostream>
#include <string>
#include <cmath>
#include <complex>

#include <il/Array.h>
#include <il/math.h>

#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>


#include "Elasticity2D.h"
#include "Mesh.h"


il::Array2D<double> rotation_matrix_2D(double theta)
{
  il::Array2D<double> R{2,2,0.};

  R(0,0)=cos(1.*theta);
  R(0,1)=-1.*sin(1.*theta);
  R(1,0)=sin(theta);
  R(1,1)=cos(theta);

  return R;
}


//----------------------------------------------------
// Segment Characteristics....
// use a struct ? (may be an oject is better ? )
// should be migrated in a separate file
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

// some utilities.
void take_submatrix(il::Array2D<double>& sub, int i0, int i1, int j0, int j1, const il::Array2D<double>& A){
  IL_ASSERT((i1 - i0+1) == sub.size(0));
  IL_ASSERT((j1-j0+1) == sub.size(1));

  for(int i = i0; i <= i1;++i) {
    for (int j=j0; j<= j1;++j){
      sub(i-i0,j-j0)=A(i,j);

    }

  }
}
////////////////////////////////////////////////////////////////////////////////

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
// analytical solution of the griffith-crack (ct pressure)
il::Array<double> griffithcrack(il::Array<double>& x, double a,double Ep, double sig)
{
  double coef =4.*sig/(Ep);
  il::Array<double> wsol{x.size(),0.};

  for (int i=0; i<x.size();++i){
    if (abs(x[i])<a) {
      wsol[i]=coef*sqrt(pow(a,2)-pow(x[i],2));
    }
  }
  return wsol;
}

////////////////////////////////////////////////////////////////////////////////
int main() {

  int n=10, p=1;
  double h=2./(n-1) ; //  element size

  il::Array<double> x{n};

  il::Array2D<double> xy{n, 2, 0.0};
  il::Array2D<int> myconn{n-1, 2, 0.0};
  il::Array2D<int> id{n-1,4,0};

  int ndof=(n-1)*4;   // number of dofs
  double Ep=1.; //Plane strain Young's modulus

  //  std::complex(double re = 0.0, double im = 0.0) myC2;
  //  myC.real(2.);
  //  myC.imag(1.);
  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create 1D mesh ....
  for (int i = 0; i < xy.size(0); ++i) {
      xy(i,0)=-1.+i*h;
      xy(i,1)=0.;
   }

  for (int  i=0; i< myconn.size(0); ++i)
  {
      myconn(i,0)=i;
      myconn(i,1)=i+1;
  }

  // create mesh object
  Mesh mesh ;
  mesh.set_values(xy,myconn);

  dofhandle_DG2D(id,mesh,p); // dof handle for DDs

  // some definitions needed for matrix assembly
  il::Array2D<double> xe{2,2,0},xec{2,2,0};

  SegmentCharacteristic mysege,mysegc;

  il::Array2D<double> K{ndof,ndof,0.},R{2,2,0.};
  il::Array<int> dofe{2*(p+1),0},dofc{2*(p+1),0};
  il::Array<double> sec{2},nec{2},xcol{2};
  il::Array2D<double>  stnl{2,4,0.};
  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout <<  "Number of dofs :" << id.size(0)*id.size(1) << "---" << (n-1)*(p+1)*2 <<"---"<< ndof <<"\n";
  std::cout << myconn.size(0)<< "\n";


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

//  std::string mystring;
//  mystring = "This is the end of this code... ";
//  std::cout << " mat entries \n" ;
//  for (int i=0;i<8; ++i){
//    std::cout << K(i,0) << " " << K(i,1) << " " << K(i,2) << " " << K(i,3) << " " << K(i,4) << " " << K(i,5) << " " << K(i,6)<< " " << K(i,7) ;
//    std::cout <<  " ...\n";
//
//  }

// solve a constant pressurized crack problem...

  il::Array<double> f{ndof,-1.};

  // just opening dds - set others to zero
  for (int i=0; i<ndof/2;++i){
    f[2*i]=0;
  }

  il::Status status;
// example if LU decomposition
//  il::LU<il::Array2D<double>> lu_decomposition(K, il::io, status);
//  if (!status.ok()) {
//    // The matrix is singular to the machine precision. You should deal with the error.
//  }
// il::Array<double> dd = lu_decomposition.solve(f);

  // use a direct solver
  il::Array<double> dd = linear_solve(K,f,il::io,status);//lu_decomposition.solve(f);


//Analytical solution at nodes
  il::Array<double> thex{ndof/2,0},wsol{ndof/2,0} ;

  int i=0;
  for (int e=0; e<n-1;++e){   // this piece of codes gets 1D mesh of x doubling the nodes of adjacent elements (for comparison with analytical solution)
    thex[i]=mesh.Coor(mesh.conn(e,0),0);
    thex[i+1]=mesh.Coor(mesh.conn(e,1),0);
    i=i+2;
  }

  wsol = griffithcrack(thex,1.,1.,1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each nodes due to the piece-wise constant nature of the solution)...
  double rel_err ;
  for (int j=0; j<ndof/2;++j){

    rel_err=sqrt(pow(dd[j*2+1]-wsol[j],2))/wsol[j];

    std::cout << "x : " << thex[j] <<"..w anal:" << wsol[j] << " w num: " << dd[j*2+1]<<  " rel error: " << rel_err << "\n";
  }

  return 0;

}