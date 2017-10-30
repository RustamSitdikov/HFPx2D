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
#include "PlaneStrainInfinite.h"
#include "Simplified3D.h"

#include <src/core/ElasticProperties.h>
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
///////////////////////////////////////////////////////////////////////////////

// todo need to write a function similar to Assembly for the addition of new
// rows and columms corresponding to the addition of new elements in the mesh !

///////////////////////////////////////////////////////////////////////////////
il::Array2D<double> basic_assembly(Mesh &mesh, ElasticProperties &elas,
                                   vKernelCall KernelCall, double ker_options) {
  // Kmat :: the stiffness matrix to assemble
  // mesh :: the Mesh object
  // id   :: the DOF handle
  // p    :: the interpolation order
  // elas :: the elastic properties object

  il::int_t p = mesh.interpolationOrder();

  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

//  hfp2d::SegmentData mysege, mysegc;

  il::StaticArray2D<double, 2, 2> R;
  il::Array<il::int_t> dofe{2 * (p + 1), 0}, dofc{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;

  il::Array2D<double> Kmat{mesh.numberOfDDDofs(),mesh.numberOfDDDofs()};

  // Brute Force assembly
  // double loop on elements to create the stiffness matrix ...

  for (il::int_t e = 0; e < mesh.numberOfElements();
       ++e) {  // loop on all  elements

    //   get characteristic of element # e
    hfp2d::SegmentData mysege = mesh.getElementData(e);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotation_matrix_2D(mysege.theta());

    // vector of dof id of  element e
    for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
      dofe[i] = mesh.dofDD(e, i);
    };

    // loop on all  elements - to compute the effect of e on all other elements
    for (il::int_t j = 0; j < mesh.numberOfElements(); ++j) {
      //   get characteristic of element # j
      hfp2d::SegmentData  mysegc = mesh.getElementData( j);

      for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
        dofc[i] = mesh.dofDD(j, i);  // vector of dof id of the  element j
      };

      // loop on collocation points of the target element
      for (il::int_t ic = 0; ic < p + 1; ++ic) {
        // call kernel fction for the effect of element e on collocation ic of
        // element j

        stnl = KernelCall(mysege, mysegc, ic, elas, ker_options);

        for (il::int_t j1 = 0; j1 < 2 * (p + 1); ++j1) {
          for (il::int_t j0 = 0; j0 < 2; ++j0) {
            Kmat(dofc[2 * ic] + j0, dofe[0] + j1) = stnl(j0, j1);
          }
        }
      }
    }
  }
  return Kmat;
};

// todo need to write a function similar to Assembly for the addition of new
// rows and columms corresponding to the addition of new elements in the mesh !
///////////////////////////////////////////////////////////////////////////////
void basic_assembly_add_elts(Mesh &new_mesh,il::int_t n_add,ElasticProperties &elas,
                             vKernelCall KernelCall, double ker_options,il::io_t,il::Array2D<double> &K){
// here newmesh containing the mesh on which the whole elasticity matrix will be built
// n_add : number of elements added - note that the new element added MUST be stored
// at the end of the newmesh connectivity table (consistency not checked)
// (or better have the new and old mesh ? )

  // resize K

  // loop on all element of new_mesh (including the added elt)

  // loop to calculate effect on the added elts.

};
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
//  tip correction for PO elements
void AddTipCorrectionP0(hfp2d::Mesh &mesh, const ElasticProperties &elas,
                        il::int_t tipElt, il::Array2D<double> &Kmat ) {

// getting the element size ;( -> cry for a method in mesh class !

//  correction factor from Ryder & Napier 1985.

  double correct =- elas.Ep()*(1. / 3.) / (4. * (mesh.elt_size(tipElt)));

  Kmat(mesh.dofDD(tipElt, 0), mesh.dofDD(tipElt, 0))+=correct;

  Kmat(mesh.dofDD(tipElt, 1), mesh.dofDD(tipElt, 1))+=correct;

}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// remove tip correction for PO elements
void RemoveTipCorrectionP0(hfp2d::Mesh &mesh, const ElasticProperties &elas,
                        il::int_t tipElt, il::Array2D<double> &Kmat ) {

//// getting the element size ;( -> cry for a method in mesh class !
//  il::StaticArray2D<double, 2, 2> Xs;
//  Xs(0, 0) = mesh.coordinates(mesh.connectivity(tipElt, 0), 0);
//  Xs(0, 1) = mesh.coordinates(mesh.connectivity(tipElt, 0), 1);
//  Xs(1, 0) = mesh.coordinates(mesh.connectivity(tipElt, 1), 0);
//  Xs(1, 1) = mesh.coordinates(mesh.connectivity(tipElt, 1), 1);
//
//  il::StaticArray<double, 2> xdiff;
//  xdiff[0] = Xs(1, 0) - Xs(0, 0);
//  xdiff[1] = Xs(1, 1) - Xs(0, 1);
//  double hx = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

//  correction factor
  double correct =- elas.Ep()*(1. / 3.) / (4. * (mesh.elt_size(tipElt)));

  Kmat(mesh.dofDD(tipElt, 0), mesh.dofDD(tipElt, 0))-=correct;

  Kmat(mesh.dofDD(tipElt, 1), mesh.dofDD(tipElt, 1))-=correct;

}


il::Array2D<double> ReArrangeKP0(const Mesh &mesh,il::Array2D<double> &Kmat) {
  // reorder K in the following blocks type
  //  Kss Ksn
  //  Kns Knn
  //  not used.
  IL_EXPECT_FAST(Kmat.size(0) == Kmat.size(1));
// test that it should even (/2)

  il::Array2D<double> Knew{Kmat.size(0), Kmat.size(1)};

  il::int_t k = 0; il::int_t l = 0;
  k=0;
  for (il::int_t i = 0; i < Kmat.size(0); i=i+2) {
    l=0;
    for (il::int_t j = 0; j < Kmat.size(1); j = j + 2){
      Knew(k,l) = Kmat(i,j);
      l++;
    }
    for (il::int_t j = 1; j < Kmat.size(1); j = j + 2){
      Knew(k,l) = Kmat(i,j);
      l++;
    }
    k++;
    };

  for (il::int_t i = 1; i < Kmat.size(0); i=i+2) {
    l=0;
    for (il::int_t j = 0; j < Kmat.size(1); j = j + 2){
      Knew(k,l) = Kmat(i,j);
      l++;
    }
    for (il::int_t j = 1; j < Kmat.size(1); j = j + 2){
      Knew(k,l) = Kmat(i,j);
      l++;
    }
    k++;
  };

  return Knew;


};

}
