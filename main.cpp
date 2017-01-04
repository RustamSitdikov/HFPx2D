
#include <iostream>
#include <string>
#include <cmath>
#include <complex>

#include <il/Array.h>
#include <il/math.h>
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include "Mesh.h"
#include "DOF_Handles.h"
#include "AssemblyDDM.h"


////////////////////////////////////////////////////////////////////////////////
// analytical solution of the griffith-crack (ct pressure)
il::Array<double> griffithcrack(il::Array<double>& x, double a,double Ep, double sig)
{
  double coef =4.*sig/(Ep);
  il::Array<double> wsol{x.size(),0.};

  for (int i=0; i<x.size();++i){
    if (std::abs(x[i])<a) {
      wsol[i]=coef*sqrt(pow(a,2)-pow(x[i],2));
    }
  }
  return wsol;
}

////////////////////////////////////////////////////////////////////////////////
int main() {

  int n=110, p=1;
  double h=2./(n-1) ; //  element size

  il::Array<double> x{n};

  il::Array2D<double> xy{n, 2, 0.0};
  il::Array2D<int> myconn{n-1, 2, 0.0};
  il::Array2D<int> id{n-1,4,0};

  int ndof=(n-1)*(p+1)*2;   // number of dofs
  double Ep=1.; //Plane strain Young's modulus

  //  std::complex(double re = 0.0, double im = 0.0) myC2;
  //  myC.real(2.);
  //  myC.imag(1.);
  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
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

  dofhandle_DG2D(id,2,n-1,p); // dof handle for DDs

  // some definitions needed for matrix assembly
  il::Array2D<double> xe{2,2,0},xec{2,2,0};

//  SegmentCharacteristic mysege,mysegc;

  il::Array2D<double> K{ndof,ndof,0.} ;

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout <<  "Number of dofs :" << id.size(0)*id.size(1) << "---" << (n-1)*(p+1)*2 <<"---"<< ndof <<"\n";
  std::cout << myconn.size(0)<< "\n";

  std::cout << "------\n";

  BasicAssembly(K, mesh, id, p , Ep ); // passing p could be avoided here.


// solve a constant pressurized crack problem...
  il::Array<double> f{ndof,-1.};
  // just opening dds - set shear loads to zero
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
  il::Array<double> dd = linear_solve(K,f,il::io,status); //lu_decomposition.solve(f);

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

