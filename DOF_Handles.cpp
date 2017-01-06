//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <il/Array2D.h>


//  FUNCTION TO CREATE A DOF HANDLE
void dofhandle_DG2D(il::Array2D<int> &dofhandle, int dof_dim, int nelts, int p) {
// function creating a matrix of dof handle - for a piece-wise linear variation per element (Discontinous Galerkin type)
// on a 1d Mesh object for the case of dof_dim Degrees of Freedoms per node
// format of the handle : number of elements \times (p+1)*dof_dim
// dof_dim :: number of dof per node
// nelts :: number of elements in the mesh
// p :: interpolation order inside each element
// needs to be generalized ....

  int ndof = nelts*(p+1)*dof_dim;

  int j ;

  for (int i = 0; i < nelts; ++i) {
    j = i*dof_dim*(p+1) ;
    for (int k=0; k<dof_dim*(p+1);++k) {
      dofhandle(i, k) = j + k;
    }
  }
//  return dofhandle; - starts at 0 for dof c++ style!
}

