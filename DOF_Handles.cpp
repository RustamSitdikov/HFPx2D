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
il::Array2D<int> dofhandle_dg2d(const int dof_dim, const int Nelts,
                                const int p) {
  // function creating a matrix of dof handle - for a piece-wise linear
  // variation per element (Discontinous Galerkin type)
  // on a 1d Mesh object for the case of dof_dim Degrees of Freedoms per node
  // format of the handle : number of elements \times (p+1)*dof_dim
  // dof_dim :: number of dof per node
  // Nelts :: number of elements in the mesh
  // p :: interpolation order inside each element

  il::Array2D<int> Dof{Nelts, 2 * dof_dim, 0};

  int j;

  for (int i = 0; i < Nelts; ++i) {
    j = i * dof_dim * (p + 1);
    for (int k = 0; k < dof_dim * (p + 1); ++k) {
      Dof(i, k) = j + k;
    }
  }

  return Dof;
}
}
