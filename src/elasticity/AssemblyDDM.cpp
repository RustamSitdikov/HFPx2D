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
#include <src/core/Utilities.h>

namespace hfp2d {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// todo need to write a function similar to Assembly for the addition of new
// rows and columms corresponding to the addition of new elements in the mesh !

///////////////////////////////////////////////////////////////////////////////
il::Array2D<double> basic_assembly(Mesh &mesh, hfp2d::ElasticProperties &elas,
                                   vKernelCall KernelCall, double ker_options) {
  // Kmat :: the stiffness matrix to assemble
  // mesh :: the Mesh object
  // id   :: the DOF handle
  // p    :: the interpolation order
  // elas :: the elastic properties object

  il::int_t p = mesh.interpolationOrder();

  il::StaticArray2D<double, 2, 2> R;
  il::Array<il::int_t> dofe{2 * (p + 1), 0}, dofc{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;

  il::Array2D<double> Kmat{mesh.numberDDDofs(), mesh.numberDDDofs()};

  // Brute Force assembly
  // double loop on elements to create the stiffness matrix ...

  for (il::int_t e = 0; e < mesh.numberOfElts(); ++e) {  // loop on all elements

    //   get characteristic of element # e
    hfp2d::SegmentData mysege = mesh.getElementData(e);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotationMatrix2D(mysege.theta());

    // vector of dof id of  element e
    for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
      dofe[i] = mesh.dofDD(e, i);
    };

    // loop on all  elements - to compute the effect of e on all other elements
    for (il::int_t j = 0; j < mesh.numberOfElts(); ++j) {
      //   get characteristic of element # j
      hfp2d::SegmentData mysegc = mesh.getElementData(j);

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

// rows and columms corresponding to the addition of new elements in the mesh !
///////////////////////////////////////////////////////////////////////////////
void basic_assembly_add_elts(Mesh &new_mesh, il::int_t n_add,
                             hfp2d::ElasticProperties &elas,
                             vKernelCall KernelCall, double ker_options,
                             il::io_t, il::Array2D<double> &K) {
  // here newmesh containing the mesh on which the whole elasticity matrix will
  // be built
  // n_add : number of elements added - note that the new element added MUST be
  // stored
  // at the end of the newmesh connectivity table (consistency not checked)
  // (or better have the new and old mesh ? )

  il::int_t p = new_mesh.interpolationOrder();

  il::int_t newtot_dof =
      new_mesh.numberDDDofsPerElt() * new_mesh.numberOfElts();

  il::int_t old_nelts = new_mesh.numberOfElts() - n_add;
  il::StaticArray2D<double, 2, 2> R;
  il::Array<il::int_t> dofe{2 * (p + 1), 0}, dofc{2 * (p + 1), 0};

  il::StaticArray2D<double, 2, 4> stnl;

  // resize K
  K.resize(newtot_dof, newtot_dof);

  for (il::int_t e = 0; e < new_mesh.numberOfElts();
       ++e) {  // loop on all  elements

    //   get characteristic of element # e
    hfp2d::SegmentData mysege = new_mesh.getElementData(e);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotationMatrix2D(mysege.theta());

    // vector of dof id of  element e
    for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
      dofe[i] = new_mesh.dofDD(e, i);
    };

    // depending if the element is added or not
    // the inner loop is over all elements or only the added ones
    il::int_t nl = new_mesh.numberOfElts();
    il::int_t j_off = 0;

    if (e < old_nelts) {  // case of already existing elements-
      // we just compute their effects on the added elements
      il::int_t nl = n_add;
      il::int_t j_off = old_nelts;
    };
    // loop to calculate effect on the added elts.
    for (il::int_t j = 0; j < nl; ++j) {
      //   get characteristic of element # j
      hfp2d::SegmentData mysegc = new_mesh.getElementData(j_off + j);

      for (il::int_t i = 0; i < 2 * (p + 1); ++i) {
        dofc[i] = new_mesh.dofDD(j_off + j,
                                 i);  // vector of dof id of the  element j
      };

      // loop on collocation points of the target element
      for (il::int_t ic = 0; ic < p + 1; ++ic) {
        // call kernel fction for the effect of element e on collocation ic of
        // element j

        stnl = KernelCall(mysege, mysegc, ic, elas, ker_options);

        for (il::int_t j1 = 0; j1 < 2 * (p + 1); ++j1) {
          for (il::int_t j0 = 0; j0 < 2; ++j0) {
            K(dofc[2 * ic] + j0, dofe[0] + j1) = stnl(j0, j1);
          }
        }
      }
    }
  }
};
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  tip correction for PO elements
void AddTipCorrectionP0(hfp2d::Mesh &mesh, const hfp2d::ElasticProperties &elas,
                        il::int_t tipElt, il::Array2D<double> &Kmat) {
  // getting the element size ;( -> cry for a method in mesh class !

  //  correction factor from Ryder & Napier 1985.

  double correct = elas.Ep() * (1. / 3.) / (4. * (mesh.eltSize(tipElt)));

  Kmat(mesh.dofDD(tipElt, 0), mesh.dofDD(tipElt, 0)) += correct;

  Kmat(mesh.dofDD(tipElt, 1), mesh.dofDD(tipElt, 1)) += correct;
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// remove tip correction for PO elements
void RemoveTipCorrectionP0(hfp2d::Mesh &mesh,
                           const hfp2d::ElasticProperties &elas,
                           il::int_t tipElt, il::Array2D<double> &Kmat) {
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
  double correct = elas.Ep() * (1. / 3.) / (4. * (mesh.eltSize(tipElt)));

  Kmat(mesh.dofDD(tipElt, 0), mesh.dofDD(tipElt, 0)) -= correct;

  Kmat(mesh.dofDD(tipElt, 1), mesh.dofDD(tipElt, 1)) -= correct;
}

il::Array2D<double> ReArrangeKP0(const Mesh &mesh, il::Array2D<double> &Kmat) {
  // reorder K in the following blocks type
  //  Kss Ksn
  //  Kns Knn
  //  not used.
  IL_EXPECT_FAST(Kmat.size(0) == Kmat.size(1));
  // test that it should even (/2)

  il::Array2D<double> Knew{Kmat.size(0), Kmat.size(1)};

  il::int_t k = 0;
  il::int_t l = 0;
  k = 0;
  for (il::int_t i = 0; i < Kmat.size(0); i = i + 2) {
    l = 0;
    for (il::int_t j = 0; j < Kmat.size(1); j = j + 2) {
      Knew(k, l) = Kmat(i, j);
      l++;
    }
    for (il::int_t j = 1; j < Kmat.size(1); j = j + 2) {
      Knew(k, l) = Kmat(i, j);
      l++;
    }
    k++;
  };

  for (il::int_t i = 1; i < Kmat.size(0); i = i + 2) {
    l = 0;
    for (il::int_t j = 0; j < Kmat.size(1); j = j + 2) {
      Knew(k, l) = Kmat(i, j);
      l++;
    }
    for (il::int_t j = 1; j < Kmat.size(1); j = j + 2) {
      Knew(k, l) = Kmat(i, j);
      l++;
    }
    k++;
  };

  return Knew;
};

////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
il::Array2D<double> basic_assembly_nodal(Mesh &mesh,
                                         hfp2d::ElasticProperties &elas,
                                         vKernelCallNode KernelCallNodal,
                                         double ker_options) {
  // this matrix perform an assembly by NODE to by element
  // such that the call to the kernel function is for 2*2 block of the matrix
  // this sub-block call will be then use in the H-mat approx.

  // Kmat :: the stiffness matrix to assemble
  // mesh :: the Mesh object
  // id   :: the DOF handle
  // p    :: the interpolation order
  // elas :: the elastic properties object

  il::int_t p = mesh.interpolationOrder();

  il::StaticArray2D<double, 2, 2> R;
  il::Array<il::int_t> dofe{2, 0}, dofc{2, 0};

  il::StaticArray2D<double, 2, 2> stnl;
  il::int_t ndof = mesh.numberDDDofs();

  il::Array2D<double> Kmat{mesh.numberDDDofs(), mesh.numberDDDofs(), 0.};

  // Brute Force assembly
  // We perform a loop elts, not node here for simplicity !  note that the
  // discretization is piece-wise linear
  // so if a node is shared by 2 segments, it has 2*2=4 unknowns associated to
  // it

  for (il::int_t e = 0; e < mesh.numberOfElts(); ++e) {  // loop on all elements

    //   get characteristic of element # e
    hfp2d::SegmentData mysege = mesh.getElementData(e);
    // Rotation matrix of the element w.r. to x-axis.
    R = hfp2d::rotationMatrix2D(mysege.theta());

    // loop on the collocation point of that element
    for (il::int_t i_c = 0; i_c < (p + 1); i_c++) {
      // vector of dof id of  element e
      // vector of dof id of element e
      for (il::int_t i = 0; i < 2; ++i) {
        dofe[i] = mesh.dofDD(e, i + 2 * i_c);
      };

      // loop on all the other elements
      for (il::int_t r = 0; r < mesh.numberOfElts(); r++) {
        //   get characteristic of element # r
        hfp2d::SegmentData mysegc = mesh.getElementData(r);

        // loop on the collocation points
        for (il::int_t j_c = 0; j_c < (p + 1); j_c++) {
          for (il::int_t i = 0; i < 2; ++i) {
            dofc[i] = mesh.dofDD(
                r, i + 2 * j_c);  // vector of dof id of the  element j
          };

          stnl = KernelCallNodal(mysege, mysegc, i_c, j_c, elas, ker_options);

          for (il::int_t j1 = 0; j1 < 2; ++j1) {
            for (il::int_t j0 = 0; j0 < 2; ++j0) {
              Kmat(dofc[0] + j0, dofe[0] + j1) = stnl(j0, j1);
            }
          }
        }
      }
    }
  }
  return Kmat;
};
}
