//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include <src/core/Fluid.h>
#include <src/core_dev/FractureEvolution.h>
#include <src/devt/FiniteVolumeRoutines.h>

namespace hfp2d {

il::Array<double> edgeConductivitiesP1Newtonian(
    Mesh &theMesh, Fluid &FluidProperties,
    il::Array<double> &permeab_middle, const il::Array<double> &dilat_middle,
    const il::Array<double> &opening_middle) {

  il::Array<double> edge_cond_p1_newt{theMesh.numberOfElts(), 0};

  for (il::int_t i = 0; i < edge_cond_p1_newt.size(); ++i) {
    edge_cond_p1_newt[i] =
        ((FluidProperties.fluidDensity() *
          ((dilat_middle[i] + opening_middle[i]) * permeab_middle[i])) /
         theMesh.eltSize(i)) *
        (1 / (12 * FluidProperties.fluidViscosity()));
  }

  return edge_cond_p1_newt;
}

/////////////// **************************** ///////////////
///////////////         FVM routines         ///////////////
/////////////// **************************** ///////////////

// Functions for the coefficients of the Finite Difference Matrix "L"
// Output: array (vector) that contains all the coefficients for each element
il::Array<double> shearConductivitiesP1Newtonian(
    Mesh &theMesh, Fluid &FluidProperties,
    FractureEvolution &FractureEvolution, const il::Array<double> &slip,
    const il::Array<double> &opening) {
  il::Array<double> wm_mid{theMesh.numberOfElts(), 0.};
  il::Array<double> wh_mid{theMesh.numberOfElts(), 0.};
  il::Array<double> kf_mid{theMesh.numberOfElts(), 0.};
  il::Array<double> whi{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> kfi{2 * theMesh.numberOfElts(), 0.};

  wm_mid = average(opening);

  whi = FractureEvolution.linearDilatancy(slip, FractureEvolution);
  wh_mid = average(whi);

  kfi = FractureEvolution.linearPermeability(slip, FractureEvolution);
  kf_mid = average(kfi);

  auto Out = edgeConductivitiesP1Newtonian(theMesh, FluidProperties, kf_mid,
                                           wh_mid, wm_mid);
  return Out;
};

///
// Function that assemble the Finite Difference matrix "L"
il::Array2D<double> buildLMatrix(Mesh &theMesh, const il::Array<double> &slip,
                                 const il::Array<double> &opening,
                                 Fluid &FluidProperties,
                                 FractureEvolution &FractureEvolution,
                                 double TimeStep) {
  il::Array<double> Kk;
  Kk = shearConductivitiesP1Newtonian(theMesh, FluidProperties,
                                      FractureEvolution, slip, opening);

  il::Array2D<double> LL{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0};
  il::Array2D<il::int_t> ed;
  il::int_t ni;
  il::int_t ej;
  il::int_t dofj = 0;
  il::Array<il::int_t> t;
  il::Array2D<il::int_t> dofhandle_press{theMesh.numberOfElts(),
                                         theMesh.interpolationOrder() + 1, 0};
  for (il::int_t i = 0; i < dofhandle_press.size(0); ++i) {
    dofhandle_press(i, 0) = i;
    dofhandle_press(i, 1) = i + 1;
  }

  // Loop over the pressure nodes
  for (int i = 0; i < LL.size(0); ++i) {
    ed = position_2d_array(dofhandle_press, i);
    ni = ed.size(0);

    for (il::int_t j = 0; j < ni; ++j) {
      ej = ed(j, 0);
      t = row_selection(dofhandle_press, ej);

      for (il::int_t k = 0; k < t.size(); ++k) {
        if (t[k] != i) dofj = t[k];
      }

      LL(i, i) = LL(i, i) - Kk[ej];
      LL(i, dofj) = LL(i, dofj) + Kk[ej];
    }
  }

  // Finally multiply the finite difference matrix by TimeStep
  il::Array2D<double> L{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0};
  for (il::int_t k = 0; k < L.size(0); ++k) {
    for (il::int_t i = 0; i < L.size(1); ++i) {
      L(k, i) = TimeStep * LL(k, i);
    }
  }

  return L;
};

////
// Function for assembling the Pressure matrix "Vp" for piecewise LINEAR DDs
// (p = 1)
il::Array2D<double> buildVpMatrix(Mesh &theMesh,
                                  FractureEvolution &FractureEvolution,
                                  Fluid &FluidProperties,
                                  const il::Array<double> &slip) {
  // Create an auxiliary vector for the assembling
  il::Array2D<il::int_t> h{2 * theMesh.numberOfElts() + 1, 2, 0};
  for (il::int_t i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{theMesh.numberOfNodes(), 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // Create all the vectors for compressibility of fluid
  // Vector of compressibility of the fluid at nodal points
  il::Array<double> Cf{vertices.size(), FluidProperties.fluidCompressibility()};
  // Vector of compressibility of the fluid at the midpoints of each element
  il::Array<double> Cfmid{theMesh.numberOfElts(),
                          FluidProperties.fluidCompressibility()};
  // Vector of compressibility of the fluid at +/- 1/4 of each element
  il::Array<double> Cfquart{2 * theMesh.numberOfElts(),
                            FluidProperties.fluidCompressibility()};

  il::Array<double> whi{2 * theMesh.numberOfElts(), 0.};
  whi = FractureEvolution.linearDilatancy(slip, FractureEvolution);

  il::Array<double> whi_mid{theMesh.numberOfElts(), 0};
  whi_mid = average(whi);

  il::Array<double> wquart{2 * theMesh.numberOfElts(), 0};
  wquart = quarter(whi);

  // Assembling the matrix
  il::Array2D<double> Vp{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0};
  il::Array2D<il::int_t> ed;
  il::Array2D<il::int_t> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::int_t dofj = 0;
  il::Array<il::int_t> t;

  /// Assembling procedure ///
  for (int m = 0; m < Vp.size(0); ++m) {
    ed = position_2d_array(theMesh.connectivity(), m);
    hi = position_2d_array(h, (2 * m));
    ni = ed.size(0);

    for (il::int_t i = 0; i < ni; ++i) {
      ej = ed(i, 0);
      hj = hi(i, 0);
      t = row_selection(theMesh.connectivity(), ej);

      for (il::int_t j = 0; j < t.size(); ++j) {
        if (t[j] != vertices[m]) dofj = t[j];
      }

      Vp(vertices[m], vertices[m]) =
          Vp(vertices[m], vertices[m]) +
          (theMesh.eltSize(ej) / 12) *
              ((whi[hj] * Cf[vertices[m]]) + (0.5 * whi_mid[ej] * Cfmid[ej]) +
               (3 * wquart[hj] * Cfquart[hj]));
      Vp(vertices[m], dofj) =
          Vp(vertices[m], dofj) +
          (theMesh.eltSize(ej) / 12) *
              ((0.5 * whi_mid[ej] * Cfmid[ej]) + (wquart[hj] * Cfquart[hj]));
    }
  }

  return Vp;
};

///
// Function for assembling the Mass/Dilatancy matrix "Vd" for piecewise LINEAR
// DDs (p = 1)
il::Array2D<double> buildVdMatrix(Mesh &theMesh,
                                  FractureEvolution &FractureEvolution,
                                  Fluid &FluidProperties,
                                  const il::Array<double> &slip) {
  // Create an auxiliary vector for the assembling
  il::Array2D<il::int_t> h{2 * theMesh.numberOfElts() + 1, 2, 0};
  for (int i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{theMesh.numberOfNodes(), 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // Create the matrix of dof for just shear DD
  // Remember: 4 DOFs per element {shear_i, normal_i, shear_j, normal_j}
  il::Array2D<il::int_t> dofhandle_single_dd{theMesh.numberOfElts(), 2, 0};
  for (il::int_t k = 0; k < dofhandle_single_dd.size(0); ++k) {
    for (il::int_t i = 0, q = 0; i < dofhandle_single_dd.size(1);
         ++i, q = q + 2) {
      dofhandle_single_dd(k, i) = theMesh.dofDD(k, q);
    }
  }

  il::Array<double> rho{2 * theMesh.numberOfElts(),
                        FluidProperties.fluidDensity()};
  il::Array<double> rho_mid{theMesh.numberOfElts(),
                            FluidProperties.fluidDensity()};
  il::Array<double> rho_quart{2 * theMesh.numberOfElts(),
                              FluidProperties.fluidDensity()};

  il::Array<double> Bi{2 * theMesh.numberOfElts(), 0.};
  Bi = FractureEvolution.derivativeLinearDilatancy(slip, FractureEvolution);

  il::Array<double> Bi_mid{theMesh.numberOfElts(), 0.};
  Bi_mid = average(Bi);

  il::Array<double> Bi_quart{2 * theMesh.numberOfElts(), 0.};
  Bi_quart = quarter(Bi);

  // Initialization
  il::Array2D<double> Vd{theMesh.numberOfNodes(), theMesh.numberDDDofs(), 0};
  il::Array2D<il::int_t> ed;
  il::Array2D<il::int_t> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::Array<il::int_t> t;

  il::int_t dofwi;
  il::int_t dofwj = 0;

  /// Assembling procedure ///
  for (int i = 0; i < Vd.size(0); ++i) {
    ed = position_2d_array(theMesh.connectivity(), i);
    hi = position_2d_array(h, (2 * i));
    ni = ed.size(0);

    for (il::int_t j = 0; j < ni; ++j) {
      ej = ed(j, 0);
      hj = hi(j, 0);

      dofwi = dofhandle_single_dd(ej, ed(j, 1));
      t = row_selection(dofhandle_single_dd, ej);

      for (il::int_t k = 0; k < t.size(); ++k) {
        if (t[k] != dofwi) dofwj = t[k];
      }

      Vd(vertices[i], dofwi) =
          Vd(vertices[i], dofwi) +
          (theMesh.eltSize(ej) / 12) *
              ((rho[hj] * Bi[hj]) + (0.5 * rho_mid[ej] * Bi_mid[ej]) +
               (3 * rho_quart[hj] * Bi_quart[hj]));
      Vd(vertices[i], dofwj) =
          Vd(vertices[i], dofwj) +
          (theMesh.eltSize(ej) / 12) * ((0.5 * rho_mid[ej] * Bi_mid[ej]) +
                                        (rho_quart[hj] * Bi_quart[hj]));
    }
  }

  return Vd;
};
}