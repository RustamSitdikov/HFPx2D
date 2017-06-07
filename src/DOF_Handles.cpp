//
// HFPx2D project.
//
// Created by Brice Lecampion on 03.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

// Inclusion from Inside Loop library
#include <il/Array2D.h>

namespace hfp2d {

//  FUNCTION TO CREATE A DOF HANDLE
il::Array2D<int> dofhandle_dp(int dof_dim, il::int_t Nelts, int p, il::io_t) {
  // function creating a matrix of dof handle - for a piece-wise linear
  // variation per element of BOTH shear AND opening DDs
  // (Discontinous Galerkin type)
  // on a 1d Mesh object for the case of dof_dim Degrees of Freedoms per node
  // format of the handle : number of elements \times (p+1)*dof_dim
  // dof_dim :: number of dof per node
  // Nelts :: number of elements in the mesh
  // p :: interpolation order inside each element
  // io_t -> everything on the left of il::io_t is read-only and is not
  //         going to be mutated

  il::Array2D<int> Dof{Nelts, 2 * dof_dim, 0};

  for (int i = 0, j; i < Nelts; ++i) {
    j = i * dof_dim * (p + 1);
    for (int k = 0; k < dof_dim * (p + 1); ++k) {
      Dof(i, k) = j + k;
    }
  }

  return Dof;
}

il::Array2D<int> dofhandle_cg(int dof_dim, il::int_t Nelts, il::io_t) {
  // function creating a matrix of dof handle - for continuous linear
  // variation per element (Continuous Galerkin type)
  // on a 1d Mesh object for the case of dof_dim Degrees of Freedoms per node
  // format of the handle : number of elements \times (p+1)*dof_dim
  // dof_dim :: number of dof per node
  // Nelts :: number of elements in the mesh
  // io_t -> everything on the left of il::io_t is read-only and is not
  //         going to be mutated

  il::Array2D<int> Dofp{Nelts, dof_dim, 0};

  for (int i = 0; i < Dofp.size(0); ++i) {

    Dofp(i, 0) = i;
    Dofp(i, 1) = i + 1;
  }

  return Dofp;
}
}
