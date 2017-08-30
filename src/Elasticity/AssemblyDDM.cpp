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
#include <il/linear_algebra.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "Simplified3D.h"
#include "PlaneStrainInfinite.h"
#include <src/Elasticity/ElasticProperties.h>

#include <src/core/Mesh.h>


namespace hfp2d {

// Some utilities //
void take_submatrix(il::Array2D<double> &sub, int i0, int i1, int j0, int j1,
                    const il::Array2D<double> &A) {
  IL_EXPECT_FAST((i1 - i0 + 1) == sub.size(0));
  IL_EXPECT_FAST((j1 - j0 + 1) == sub.size(1));

  for (int i = i0; i <= i1; ++i) {
    for (int j = j0; j <= j1; ++j) {
      sub(i - i0, j - j0) = A(i, j);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

void set_submatrix(il::Array2D<double> &A, int i0, int i1,
                   const il::StaticArray2D<double, 2, 4> &B) {
  IL_EXPECT_FAST(i0 + B.size(0) <= A.size(0));
  IL_EXPECT_FAST(i1 + B.size(1) <= A.size(1));

  for (int j1 = 0; j1 < B.size(1); ++j1) {
    for (int j0 = 0; j0 < B.size(0); ++j0) {
      A(i0 + j0, i1 + j1) = B(j0, j1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
il::Array2D<double> basic_assembly(Mesh &mesh, il::Array2D<int> &id, int p,
                                   ElasticProperties& elas){
  // Kmat : the stiffness matrix to assemble
  // mesh:: the Mesh object
  // id :: the DOF handle
  // p :: the interpolation order
  // Ep :: the Plane Strain Young's modulus
  IL_EXPECT_FAST(id.size(0) == mesh.nelts());
  IL_EXPECT_FAST(id.size(1) == 2 * (p + 1));


  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  hfp2d::SegmentData mysege, mysegc;

  il::StaticArray2D<double, 2, 2> R;
  il::Array<int> dofe{2 * (p + 1), 0}, dofc{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;
  il::StaticArray<double, 2> sec, nec, xcol;

  il::Array2D<double> Kmat{id.size(0) * id.size(1), id.size(0) * id.size(1)};

  // Brute Force assembly
  // double loop on elements to create the stiffness matrix ...

  #pragma omp parallel for

  for (int e = 0; e < mesh.nelts(); ++e) { // loop on all  elements

    //   get characteristic of element # e
    mysege = hfp2d::get_segment_DD_data(mesh, e, p);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotation_matrix_2D(mysege.theta);

    // vector of dof id of  element e
    for (int i = 0; i < 2 * (p + 1); ++i) {
      dofe[i] = id(e, i);
    };

    // loop on all  elements - to compute the effect of e on all other elements
    for (int j = 0; j < mesh.nelts(); ++j) {
      //   get characteristic of element # j
      mysegc = hfp2d::get_segment_DD_data(mesh, j, p);

      sec = il::dot(R, mysegc.s); // tangent of elt j in the frame of elt e
      nec = il::dot(R, mysegc.n); // normal of elt j in the frame of elt e

      for (int i = 0; i < 2 * (p + 1); ++i) {
        dofc[i] = id(j, i); // vector of dof id of the  element j
      };

      // loop on collocation points of the target element
      for (int ic = 0; ic < p + 1; ++ic) {
        // we switch to the frame of element e
        for (int i = 0; i < 2; ++i) {
          xcol[i] = mysegc.CollocationPoints(ic, i) - mysege.Xmid[i];
        }
        xcol = il::dot(R, xcol);

        // TODO: Make the call Kernel agnostic..... and add a virtual fction call
        // call kernel fction for the effect of element e on collocation ic of element j
        stnl = hfp2d::normal_shear_stress_kernel_dp1_dd(xcol, mysege.size, sec,
                                                        nec, elas.Ep());

        hfp2d::set_submatrix(Kmat, dofc[2 * ic], dofe[0], stnl);
      }
    }
  }
  return Kmat;
};
}
