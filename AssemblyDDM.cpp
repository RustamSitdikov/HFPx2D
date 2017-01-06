//
// Created by Brice Lecampion on 03.01.17.
//

#include "AssemblyDDM.h"
#include "Mesh.h"
#include "Elasticity2D.h"
#include<il/linear_algebra.h>

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

void set_submatrix(il::Array2D<double>& A, int i0, int i1, const il::StaticArray2D<double,2,4>& B) {
  IL_ASSERT(i0 + B.size(0) <= A.size(0));
  IL_ASSERT(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
} // e.g. set_submatrix(A, 2, 3, B);



void BasicAssembly(il::Array2D<double> &Kmat, Mesh mesh, il::Array2D<int> id, int p , double Ep ){
  // Kmat : the stiffness matrix to assemble
  // mesh:: the Mesh object
  // id :: the DOF handle
  // p :: the interpolation order
  // Ep :: the Plane Strain Young's modulus
  IL_ASSERT(id.size(0) == mesh.nelts());
  IL_ASSERT(id.size(1) == 2*(p+1));
  IL_ASSERT(Kmat.size(0) ==Kmat.size(1));
  IL_ASSERT(Kmat.size(0) == id.size(0)*id.size(1) );


  il::Array2D<double> xe{2,2,0},xec{2,2,0};

  SegmentCharacteristic mysege,mysegc;

  il::StaticArray2D<double,2,2>  R ;
  il::Array<int> dofe{2*(p+1),0},dofc{2*(p+1),0};

  il::StaticArray2D<double,2,4>  stnl ;
  il::StaticArray<double,2> sec,nec,xcol;

// Brute Force assembly
// double loop on elements to create the stiffness matrix...
   for (int e=0;e<mesh.nelts();++e){ // loop on all  elements

    take_submatrix(xe,mesh.conn(e,0),mesh.conn(e,1),0,1,mesh.Coor); // take the coordinates of element e from the mesh object
    mysege=get_segment_DD_characteristic(xe,p); // get the segment characteristic.
    R=rotation_matrix_2D(mysege.theta);      //Rotation matrix of the element w.r. to x-axis.

    for (int i=0;i<2*(p+1);++i){ // vector of dof id of the element e
      dofe[i]=id(e,i);
    };

     for (int j=0;j<mesh.nelts();++j){// loop on all  elements

      take_submatrix(xec,mesh.conn(j,0),mesh.conn(j,1),0,1,mesh.Coor); // takes the coordinates of element j
      mysegc=get_segment_DD_characteristic(xec,p);

      sec=il::dot(R,mysegc.s); // tangent of elt j
      nec=il::dot(R,mysegc.n); // normal of elt j

      for (int i=0;i<2*(p+1);++i){
        dofc[i]=id(j,i); // vector of dof id of the  element j
      };

      for (int ic=0;ic<p+1;++ic){ // loop on collocation points of the target element
        // we switch to the frame of element e
        for (int i=0;i<2;++i){
          xcol[i]=mysegc.CollocationPoints(ic,i)-mysege.Xmid[i]; //
        }
        xcol=il::dot(R,xcol);

        NormalShearStressKernel_LinearDD(stnl,xcol,mysege.size,sec,nec,Ep);
        //

        set_submatrix(Kmat,dofc[2*ic],dofe[0],stnl);

      }
    }
  }



};
