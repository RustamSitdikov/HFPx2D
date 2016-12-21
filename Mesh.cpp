//
// Created by Brice Lecampion on 10.12.16.
//

#include "Mesh.h"

//mesh class
void Mesh::set_values(il::Array2D<double> xy, il::Array2D<int> ien)
{
  IL_ASSERT(xy.size(1) ==2); //check array dimensions ?
  IL_ASSERT(ien.size(1) ==2);//check array dimensions ??? -> this is only for 1D mesh so far

  Coor = xy;   // list of coordinates of points in the mesh
  conn = ien;  //  connectivity array -
}
// needs to add function to add one or more elements ... (needs to have active and passive elements to track active/passive fractures etc.)
//
/// // could provide a default constructor for a straight fracture ?








//  FUNCTION TO CREATE A DOF HANDLES -> Put in another file ?
void dofhandle_DG2D(il::Array2D<int> &dofhandle, Mesh mesh, int p) {
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
