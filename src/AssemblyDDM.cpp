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
#include <il/math.h>
// Inclusion from the project
#include "AssemblyDDM.h"
#include "Elasticity2D.h"
//#include "Mesh.h"

namespace hfp2d {

// some utilities.
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
}  // e.g. set_submatrix(A, 2, 3, B);

il::Array2D<double> basic_assembly(Mesh &mesh, il::Array2D<int> &id, int p,
                                   double Ep) {
  // Kmat : the stiffness matrix to assemble
  // mesh:: the Mesh object
  // id :: the DOF handle
  // p :: the interpolation order
  // Ep :: the Plane Strain Young's modulus
  IL_EXPECT_FAST(id.size(0) == mesh.nelts());
  IL_EXPECT_FAST(id.size(1) == 2 * (p + 1));
  //  IL_EXPECT_FAST(Kmat.size(0) == Kmat.size(1));
  //  IL_EXPECT_FAST(Kmat.size(0) == id.size(0) * id.size(1));

  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  hfp2d::SegmentCharacteristic mysege, mysegc;
  il::Array2D<double> Kmat{id.size(0) * id.size(1), id.size(0) * id.size(1)};

  il::StaticArray2D<double, 2, 2> R;
  il::Array<int> dofe{2 * (p + 1), 0}, dofc{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;
  il::StaticArray<double, 2> sec, nec, xcol;

  // Brute Force assembly
  // double loop on elements to create the stiffness matrix ...
  for (il::int_t e = 0; e < mesh.nelts(); ++e) {  // loop on all  elements

    //   get characteristic of element # e
    mysege = hfp2d::get_segment_DD_characteristic(mesh, e, p);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotation_matrix_2D(mysege.theta);

    for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
      // vector of dof id of the element e
      dofe[i] = id(e, i);
    };

    // loop on all  elements
    for (il::int_t j = 0; j < mesh.nelts(); ++j) {
      //   get characteristic of element # j
      mysegc = hfp2d::get_segment_DD_characteristic(mesh, j, p);

      sec = il::dot(R, mysegc.s);  // tangent of elt j
      nec = il::dot(R, mysegc.n);  // normal of elt j

      for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
        dofc[i] = id(j, i);  // vector of dof id of the  element j
      };
      // loop on collocation points of the target element
      for (il::int_t ic = 0; ic < p + 1; ++ic) {
        // we switch to the frame of element e
        for (il::int_t i = 0; i < 2; ++i) {
          xcol[i] = mysegc.CollocationPoints(ic, i) - mysege.Xmid[i];
        };

        xcol = il::dot(R, xcol);

        stnl = hfp2d::normal_shear_stress_kernel_dp1_dd(xcol, mysege.size, sec,
                                                        nec, Ep);
     //   hfp2d::set_submatrix(Kmat, dofc[2 * ic], dofe[0], stnl);

        for (il::int_t j1 = 0; j1 < 4; ++j1) {
          for (il::int_t j0 = 0; j0 < 2; ++j0) {
            Kmat(dofc[2 * ic] + j0, dofe[0] + j1) = stnl(j0, j1);
           }
        }

      }
    }
  }
  return Kmat;
};
}
