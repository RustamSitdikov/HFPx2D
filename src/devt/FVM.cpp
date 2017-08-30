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
#include "ConductivitiesNewtonian.h"
#include "src/core/DOF_Handles.h"
#include "EHLDs.h"

namespace hfp2d {

/////// Some utilities ///////

/// 1
// This function calculates the average between two values for each element
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector that contains the average values of each row (element)

il::Array<double> average(const il::Array2D<double> &d, il::io_t) {

  il::Array<double> Average{d.size(0), 0};
  for (il::int_t i = 0; i < d.size(0); ++i) {

    Average[i] = (d(i, 0) + d(i, 1)) / 2;
  }
  return Average;
};

/// 2
// This function calculates the slip/opening at +/- 1/4 -> the control volume is
// centered on the nodes!
// Input: matrix of slip/opening for each element -> size Nelts x 2
// Remember: piecewise linear variation over the element
// Output: vector -> {slip_+1/4 , slip_+3/4}

il::Array<double> quarter(const il::Array2D<double> &d, il::io_t) {

  il::Array<double> Quarter(2 * d.size(0), 0);
  for (il::int_t i = 0, j = 0; i < (d.size(0)); ++i, j = j + 2) {

    Quarter[j] = ((3 * d(i, 0)) + d(i, 1)) / 4;
    Quarter[j + 1] = (d(i, 0) + (3 * d(i, 1))) / 4;
  }

  return Quarter;
};

/// 3
// Function to find out the position of a value in a 2D array
// It returns 2x2 array with row&col of the seek value
// It is completely general in a sense that the output can be a vector or a
// matrix (2x2)
il::Array2D<int> position_2d_array(const il::Array2D<int> &arr2D, int seek,
                                   il::io_t) {

  il::Array2D<int> M{arr2D.size(1) * arr2D.size(0), 2, -1};
  int k = 0;

  for (int i = 0; i < arr2D.size(0); ++i) {

    for (int j = 0; j < arr2D.size(1); ++j) {

      if (arr2D(i, j) == seek) {

        M(k, 0) = i;
        M(k, 1) = j;

        k = k + 1;
      }
    }
  }

  il::Array2D<int> outp{k, 2, 0};

  for (int l = 0; l < k; ++l) {

    for (int j = 0; j < 2; ++j) {

      outp(l, j) = M(l, j);
    }
  }

  return outp;
};

/// 4
// Function to find out the position of a value in a 2D array in a more
// efficient way
// It is completely general in a sense that the output can be a vector or a
// matrix (2x2)
il::Array2D<int> search(const il::Array2D<int> &matrix, int x, il::io_t) {
  il::Array2D<int> ans{0, 2};
  il::int_t k = 0;

  for (int j = 0; j < matrix.size(1); ++j) {
    for (int i = 0; i < matrix.size(0); ++i) {
      if (matrix(i, j) == x) {
        ans.resize(k + 1, 2);
        ans(k, 0) = i;
        ans(k, 1) = j;
        ++k;
      }
    }
  }

  return ans;
}

/// 5
// Auxiliary function for assembly process
// It returns a given row (vector - specified by idx) of a 2D array
il::Array<int> row_selection(const il::Array2D<int> &arr, il::int_t idx,
                             il::io_t) {

  il::Array<int> vect{arr.size(1), 0};
  for (il::int_t i = 0; i < vect.size(); ++i) {

    vect[i] = arr(idx, i);
  }

  return vect;
};

/////////////// **************************** ///////////////
///////////////         FVM routines         ///////////////
/////////////// **************************** ///////////////

// Functions for the coefficients of the Finite Difference Matrix "L"
// Output: array (vector) that contains all the coefficients for each element
il::Array<double> shear_conductivities_newtonian(
    Parameters_fluid &fluid_parameters, Mesh mesh, const il::Array2D<double> &d,
    Parameters_dilatancy &dilat_parameters,
    Parameters_permeability &permeab_parameters, il::io_t) {

  // Inputs:
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - mesh -> mesh class
  //  - d -> matrix of slip at nodes {{d1_left,d1_right},{d2_left,d2_right}..}
  //  (size -> Nelts x 2)
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - permeab_parameters -> structure that contains all the permeability
  //  parameters
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> d_mid;
  il::Array<double> wh_mid;
  il::Array<double> rho_mid;
  il::Array<double> kf_mid;

  d_mid = average(d, il::io);
  wh_mid = dilatancy(dilat_parameters, d_mid, il::io);
  rho_mid = average(fluid_parameters.density, il::io);
  kf_mid = permeability(permeab_parameters, d_mid, il::io);

  // create the array of element size
  il::Array<double> EltSizes{mesh.nelts(), 0.};
  for (il::int_t i = 0; i < EltSizes.size(); ++i) {
    for (il::int_t j = 0; j < (mesh.conn()).size(1); ++j) {

      EltSizes[i] =
          fabs(fabs(EltSizes[i]) - fabs(mesh.node(mesh.connectivity(i, j), 0)));
    }
  }

  il::Array<double> Out{mesh.nelts(), 0};
  Out = conductivities_newtonian(rho_mid, wh_mid, EltSizes, fluid_parameters,
                                 kf_mid, il::io);

  return Out;
};

///
// Function that assemble the Finite Difference matrix "L"
il::Array2D<double> build_l_matrix(Mesh mesh, const il::Array2D<double> &d,
                                   Parameters_fluid &fluid_parameters,
                                   Parameters_dilatancy &dilat_parameters,
                                   const double &TimeStep,
                                   Parameters_permeability &permeab_parameters,
                                   il::io_t) {

  // Inputs:
  //  - mesh -> mesh class
  //  - d -> matrix of shear DD or opening DD at nodes
  //  {{d1_left,d1_right},{d2_left,d2_right}..} (size -> Nelts x 2)
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - TimeStep -> time step
  //  - permeab_paramters -> structure that contains all the permeability
  //  parameters
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  // Output:
  //  - Finite difference matrix "L" (size -> Nnodes x Nnodes)

  il::Array<double> Kk;
  Kk = shear_conductivities_newtonian(
      fluid_parameters, mesh, d, dilat_parameters, permeab_parameters, il::io);

  il::Array2D<double> LL{mesh.nelts() + 1, mesh.nelts() + 1, 0};
  il::Array2D<int> ed;
  il::int_t ni;
  il::int_t ej;
  il::int_t dofj;
  il::Array<int> t;
  il::Array2D<int> Dofp;
  Dofp = hfp2d::dofhandle_cg(2, mesh.nelts(), il::io);

  // Loop over the pressure nodes
  for (il::int_t i = 0; i < LL.size(0); ++i) {

    ed = position_2d_array(Dofp, ((double)i), il::io);
    ni = ed.size(0);

    for (il::int_t j = 0; j < ni; ++j) {

      ej = ed(j, 0);
      t = row_selection(Dofp, ej, il::io);

      for (il::int_t k = 0; k < t.size(); ++k) {

        if (t[k] != i)
          dofj = t[k];
      }

      LL(i, i) = LL(i, i) - Kk[ej];
      LL(i, dofj) = LL(i, dofj) + Kk[ej];
    }
  }

  // Finally multiply the finite difference matrix by TimeStep
  il::Array2D<double> L{mesh.nelts() + 1, mesh.nelts() + 1, 0};
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
il::Array2D<double> build_vp_matrix_p1(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array2D<double> &d, il::io_t) {

  // Inputs:
  //  - mesh -> mesh class
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - d -> matrix of shear DD at nodes {{d1_left, d1_right},{d2_left,
  //  d2_right}..} (size -> Nelts x 2)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  // Output:
  //  - Pressure matrix "Vp", whose size is Nnodes x Nnodes. It is a four banded
  //  diagonal matrix

  // Create an auxiliary vector for the assembling
  il::Array2D<int> h{2 * mesh.nelts() + 1, 2, 0};
  for (int i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{mesh.nelts() + 1, 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // create the array of element size
  il::Array<double> EltSizes{mesh.nelts(), 0};
  for (il::int_t i = 0; i < EltSizes.size(); ++i) {
    for (il::int_t j = 0; j < (mesh.conn()).size(1); ++j) {

      EltSizes[i] =
          fabs(fabs(EltSizes[i]) - fabs(mesh.node(mesh.connectivity(i, j), 0)));
    }
  }

  // Create all the vectors for compressibility of fluid
  // Vector of compressibility of the fluid at nodal points
  il::Array<double> Cf{vertices.size(), fluid_parameters.compressibility};
  // Vector of compressibility of the fluid at the midpoints of each element
  il::Array<double> Cfmid{mesh.nelts(), fluid_parameters.compressibility};
  // Vector of compressibility of the fluid at +/- 1/4 of each element
  il::Array<double> Cfquart{2 * mesh.nelts(), fluid_parameters.compressibility};

  il::Array<double> d_left{mesh.nelts(), 0};
  for (il::int_t k = 0; k < d_left.size(); ++k) {
    d_left[k] = d(k, 0);
  }

  il::Array<double> d_right{(mesh.conn()).size(0), 0};
  for (il::int_t j = 0; j < d_right.size(); ++j) {
    d_right[j] = d(j, 1);
  }

  il::Array<double> whi_left{mesh.nelts(), 0};
  il::Array<double> whi_right{mesh.nelts(), 0};

  whi_left = dilatancy(dilat_parameters, d_left, il::io);
  whi_right = dilatancy(dilat_parameters, d_right, il::io);

  il::Array2D<double> whi{mesh.nelts(), 2, 0};
  for (il::int_t l = 0; l < whi.size(0); ++l) {

    whi(l, 0) = whi_left[l];
    whi(l, 1) = whi_right[l];
  }

  il::Array<double> whi_mid{mesh.nelts(), 0};
  whi_mid = average(whi, il::io);

  il::Array<double> wquart{2 * mesh.nelts(), 0};
  wquart = quarter(whi, il::io);

  // Assembling the matrix
  il::Array2D<double> Vp{mesh.nelts() + 1, mesh.nelts() + 1, 0};
  il::Array2D<int> ed;
  il::Array2D<int> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::int_t dofj;
  il::Array<int> t;

  // Vector that we need for assembling process
  il::Array<double> Whi{2 * mesh.nelts(), 0};
  for (il::int_t n = 0, q = 0; n < whi_left.size(); ++n, q = q + 2) {

    Whi[q] = whi_left[n];
    Whi[q + 1] = whi_right[n];
  }

  /// Assembling procedure ///
  for (il::int_t m = 0; m < Vp.size(0); ++m) {

    ed = position_2d_array(mesh.conn(), ((double)m), il::io);
    hi = position_2d_array(h, ((double)(2 * m)), il::io);
    ni = ed.size(0);

    for (il::int_t i = 0; i < ni; ++i) {

      ej = ed(i, 0);
      hj = hi(i, 0);
      t = row_selection(mesh.conn(), ej, il::io);

      for (il::int_t j = 0; j < t.size(); ++j) {

        if (t[j] != vertices[m])
          dofj = t[j];
      }

      Vp(vertices[m], vertices[m]) =
          Vp(vertices[m], vertices[m]) +
          (EltSizes[ej] / 12) *
              ((Whi[hj] * Cf[vertices[m]]) + (0.5 * whi_mid[ej] * Cfmid[ej]) +
               (3 * wquart[hj] * Cfquart[hj]));
      Vp(vertices[m], dofj) =
          Vp(vertices[m], dofj) +
          (EltSizes[ej] / 12) *
              ((0.5 * whi_mid[ej] * Cfmid[ej]) + (wquart[hj] * Cfquart[hj]));
    }
  }

  return Vp;
};

///
// Function for assembling the Mass matrix "Vd" for piecewise LINEAR DDs (p = 1)
il::Array2D<double> build_vd_matrix_p1(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       il::Array2D<int> Dof,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array2D<double> &d, il::io_t) {

  // Inputs:
  //  - mesh -> mesh class
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - Dof -> matrix that contains all the degrees of freedom (Remember 4 DOFs
  //  per element {shear_i,normal_i,shear_j,normal_j}) (size -> Nelts x 4)
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - d -> matrix of shear DD at nodes {{d1_left, d1_right},{d2_left,
  //  d2_right}..} (size -> Nelts x 2)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  // Output:
  //  - Mass matrix "Vd". (size -> Nnodes x 2*Nelts). It is a four banded
  //  diagonal matrix

  // Create an auxiliary vector for the assembling
  il::Array2D<int> h{2 * mesh.nelts() + 1, 2, 0};
  for (int i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{mesh.nelts() + 1, 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // create the array of element size
  il::Array<double> EltSizes{mesh.nelts(), 0};
  for (il::int_t i = 0; i < EltSizes.size(); ++i) {
    for (il::int_t j = 0; j < (mesh.conn()).size(1); ++j) {

      EltSizes[i] =
          fabs(fabs(EltSizes[i]) - fabs(mesh.node(mesh.connectivity(i, j), 0)));
    }
  }

  // Create the matrix of dof for just shear DD
  // Remember: 4 DOFs per element {shear_i, normal_i, shear_j, normal_j}
  il::Array2D<int> dofw{mesh.nelts(), 2, 0};
  for (il::int_t k = 0; k < dofw.size(0); ++k) {
    for (il::int_t i = 0, q = 0; i < dofw.size(1); ++i, q = q + 2) {

      dofw(k, i) = Dof(k, q);
    }
  }

  il::Array<double> rho_mid{mesh.nelts(), 0};
  il::Array<double> rho_quart{2 * mesh.nelts(), 0};
  rho_mid = average(fluid_parameters.density, il::io);
  rho_quart = quarter(fluid_parameters.density, il::io);

  il::Array<double> d_left{mesh.nelts(), 0};
  for (il::int_t k = 0; k < d_left.size(); ++k) {
    d_left[k] = d(k, 0);
  }

  il::Array<double> d_right{mesh.nelts(), 0};
  for (il::int_t j = 0; j < d_right.size(); ++j) {
    d_right[j] = d(j, 1);
  }

  il::Array<double> Bi_left{mesh.nelts(), 0};
  il::Array<double> Bi_right{mesh.nelts(), 0};

  Bi_left = d_dilatancy(dilat_parameters, d_left, il::io);
  Bi_right = d_dilatancy(dilat_parameters, d_right, il::io);

  il::Array2D<double> Bi{mesh.nelts(), 2, 0};
  for (il::int_t l = 0; l < Bi.size(0); ++l) {
    Bi(l, 0) = Bi_left[l];
    Bi(l, 1) = Bi_right[l];
  }

  il::Array<double> Bi_mid{mesh.nelts(), 0};
  Bi_mid = average(Bi, il::io);

  il::Array<double> Bi_quart{2 * mesh.nelts(), 0};
  Bi_quart = quarter(Bi, il::io);

  // Initialization
  il::Array2D<double> Vd{mesh.nelts() + 1, 4 * mesh.nelts(), 0};
  il::Array2D<int> ed;
  il::Array2D<int> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::Array<int> t;

  il::int_t dofwi;
  il::int_t dofwj;

  // We need for the assembling
  il::Array<double> Rho{2 * mesh.nelts(), 0};
  for (il::int_t m = 0, q = 0; m < fluid_parameters.density.size(0);
       ++m, q = q + 2) {
    Rho[q] = fluid_parameters.density(m, 0);
    Rho[q + 1] = fluid_parameters.density(m, 1);
  }

  // We need for the assembling
  il::Array<double> BI{2 * mesh.nelts(), 0};
  for (il::int_t m = 0, q = 0; m < Bi.size(0); ++m, q = q + 2) {
    BI[q] = Bi(m, 0);
    BI[q + 1] = Bi(m, 1);
  }

  /// Assembling procedure ///
  for (il::int_t i = 0; i < Vd.size(0); ++i) {

    ed = position_2d_array(mesh.conn(), ((double)i), il::io);
    hi = position_2d_array(h, ((double)(2 * i)), il::io);
    ni = ed.size(0);

    for (il::int_t j = 0; j < ni; ++j) {

      ej = ed(j, 0);
      hj = hi(j, 0);

      dofwi = dofw(ej, ed(j, 1));
      t = row_selection(dofw, ej, il::io);

      for (il::int_t k = 0; k < t.size(); ++k) {

        if (t[k] != dofwi)
          dofwj = t[k];
      }

      Vd(vertices[i], dofwi) =
          Vd(vertices[i], dofwi) +
          (EltSizes[ej] / 12) *
              ((Rho[hj] * BI[hj]) + (0.5 * rho_mid[ej] * Bi_mid[ej]) +
               (3 * rho_quart[hj] * Bi_quart[hj]));
      Vd(vertices[i], dofwj) =
          Vd(vertices[i], dofwj) +
          (EltSizes[ej] / 12) * ((0.5 * rho_mid[ej] * Bi_mid[ej]) +
                                 (rho_quart[hj] * Bi_quart[hj]));
    }
  }

  //  // Taking finally only the shear related contributions...
  //  il::Array2D<double> VD{Vd.size(0), Vd.size(1)/2, 0.};
  //
  //  for (il::int_t n = 0; n < VD.size(0); ++n) {
  //
  //    for (il::int_t i = 0, q = 0; i < VD.size(1); ++i, q = q + 2) {
  //
  //      VD(n, i) = Vd(n, q);
  //    }
  //  }

  return Vd;
};

////
// Function for assembling the Pressure matrix "Vp" for piecewise CONSTANT DDs
// (p = 0)
il::Array2D<double> build_vp_matrix_p0(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - mesh -> mesh class
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - d -> vector of shear DD for each element (Remember: piecewise CONSTANT
  //  DDs, so d_{i} = d_{i+1/2} = d_{i+1/4})  (size -> Nelts)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  // Output:
  //  - Pressure matrix "Vp", whose size is Nnodes x Nnodes. It is a four banded
  //  diagonal matrix

  // Create an auxiliary vector for the assembling
  il::Array2D<int> h{2 * mesh.nelts() + 1, 2, 0};
  for (int i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{mesh.nelts() + 1, 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // create the array of element size
  il::Array<double> EltSizes{mesh.nelts(), 0};
  for (il::int_t i = 0; i < EltSizes.size(); ++i) {
    for (il::int_t j = 0; j < (mesh.conn()).size(1); ++j) {

      EltSizes[i] =
          fabs(fabs(EltSizes[i]) - fabs(mesh.node(mesh.connectivity(i, j), 0)));
    }
  }

  // Create all the vectors for compressibility of fluid
  // Vector of compressibility of the fluid at nodal points
  il::Array<double> Cf{vertices.size(), fluid_parameters.compressibility};
  // Vector of compressibility of the fluid at the midpoints of each element
  il::Array<double> Cfmid{mesh.nelts(), fluid_parameters.compressibility};
  // Vector of compressibility of the fluid at +/- 1/4 of each element
  il::Array<double> Cfquart{2 * mesh.nelts(), fluid_parameters.compressibility};

  il::Array<double> whi{mesh.nelts(), 0};

  // wh_{i} = wh_{imid} = wh_{iquart}, because p=0
  whi = dilatancy(dilat_parameters, d, il::io);

  // Initialization
  il::Array2D<double> Vp{mesh.nelts() + 1, mesh.nelts() + 1, 0};
  il::Array2D<int> ed;
  il::Array2D<int> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::int_t dofj;
  il::Array<int> t;

  /// Assembling procedure ///
  for (il::int_t m = 0; m < Vp.size(0); ++m) {

    ed = position_2d_array(mesh.conn(), ((double)m), il::io);
    hi = position_2d_array(h, ((double)(2 * m)), il::io);
    ni = ed.size(0);

    for (il::int_t i = 0; i < ni; ++i) {

      ej = ed(i, 0);
      hj = hi(i, 0);
      t = row_selection(mesh.conn(), ej, il::io);

      for (il::int_t j = 0; j < t.size(); ++j) {

        if (t[j] != vertices[m])
          dofj = t[j];
      }

      Vp(vertices[m], vertices[m]) =
          Vp(vertices[m], vertices[m]) +
          ((EltSizes[ej] * whi[ej]) / 12) *
              ((Cf[vertices[m]]) + (0.5 * Cfmid[ej]) + (3 * Cfquart[hj]));
      Vp(vertices[m], dofj) =
          Vp(vertices[m], dofj) +
          ((EltSizes[ej] * whi[ej]) / 12) * ((0.5 * Cfmid[ej]) + (Cfquart[hj]));
    }
  }

  return Vp;
};

///
// Function for assembling the Mass matrix "Vd" for piecewise CONSTANT DDs
// (p = 0)
il::Array2D<double> build_vd_matrix_p0(Mesh mesh,
                                       Parameters_dilatancy &dilat_parameters,
                                       il::Array2D<int> &Dof,
                                       Parameters_fluid &fluid_parameters,
                                       const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - mesh -> mesh class
  //  - dilat_parameters -> structure that contains all the dilatancy parameters
  //  - Dof -> matrix that contains all the degrees of freedom (Remember 4 DOFs
  //  per element {shear_i,normal_i,shear_j,normal_j}) (size -> Nelts x 4)
  //  - fluid_parameters -> structure that contains all the fluid parameters
  //  - d -> vector of shear DD for each element (Remember: piecewise CONSTANT
  //  DDs, so d_{i} = d_{i+1/2} = d_{i+1/4})  (size -> Nelts)
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  // Output:
  //  - Mass matrix "Vd". (size -> Nnodes x 2*Nelts). It is a four banded
  //  diagonal matrix

  // Create an auxiliary vector for the assembling
  il::Array2D<int> h{2 * mesh.nelts() + 1, 2, 0};
  for (int i = 0; i < h.size(0); ++i) {
    h(i, 0) = i;
    h(i, 1) = i + 1;
  }

  // mesh nodes {0,1,2,3..}
  il::Array<int> vertices{mesh.nelts() + 1, 0};
  for (int j = 0; j < vertices.size(); ++j) {
    vertices[j] = j;
  }

  // create the array of element size
  il::Array<double> EltSizes{mesh.nelts(), 0};
  for (il::int_t i = 0; i < EltSizes.size(); ++i) {
    for (il::int_t j = 0; j < (mesh.conn()).size(1); ++j) {
      EltSizes[i] =
          fabs(fabs(EltSizes[i]) - fabs(mesh.node(mesh.connectivity(i, j), 0)));
    }
  }

  // Create the matrix of dof for just shear DD
  // Remember: 4 DOFs per element {shear_i, normal_i, shear_j, normal_j}
  il::Array2D<int> dofw{mesh.nelts(), 2, 0};
  for (il::int_t k = 0; k < dofw.size(0); ++k) {
    for (il::int_t i = 0, q = 0; i < dofw.size(1); ++i, q = q + 2) {

      dofw(k, i) = Dof(k, q);
    }
  }

  il::Array<double> rho_mid{mesh.nelts(), 0};
  il::Array<double> rho_quart{2 * mesh.nelts(), 0};

  rho_mid = average(fluid_parameters.density, il::io);
  rho_quart = quarter(fluid_parameters.density, il::io);

  il::Array<double> Bi{mesh.nelts(), 0};

  Bi = d_dilatancy(dilat_parameters, d, il::io);

  // Initialization
  il::Array2D<double> Vd{mesh.nelts() + 1, 4 * mesh.nelts(), 0};
  il::Array2D<int> ed;
  il::Array2D<int> hi;
  il::int_t ni;
  il::int_t ej;
  il::int_t hj;
  il::Array<int> t;

  il::int_t dofwi;
  il::int_t dofwj;

  // We need for the assembling
  il::Array<double> Rho{2 * (mesh.conn()).size(0), 0};
  //  il::Array<double> BI{2 * (mesh.conn).size(0),0};

  for (il::int_t m = 0, q = 0; m < fluid_parameters.density.size(0);
       ++m, q = q + 2) {
    Rho[q] = fluid_parameters.density(m, 0);
    Rho[q + 1] = fluid_parameters.density(m, 1);
  }

  /// Assembling procedure ///

  for (il::int_t i = 0; i < Vd.size(0); ++i) {

    ed = position_2d_array(mesh.conn(), ((double)i), il::io);
    hi = position_2d_array(h, ((double)(2 * i)), il::io);
    ni = ed.size(0);

    for (il::int_t j = 0; j < ni; ++j) {

      ej = ed(j, 0);
      hj = hi(j, 0);

      dofwi = dofw(ej, ed(j, 1));

      t = row_selection(dofw, ej, il::io);

      for (il::int_t k = 0; k < t.size(); ++k) {

        if (t[k] != dofwi)
          dofwj = t[k];
      }

      Vd(vertices[i], dofwi) =
          Vd(vertices[i], dofwi) +
          ((EltSizes[ej] * Bi[ej]) / 12) *
              ((Rho[hj]) + (0.5 * rho_mid[ej]) + (3 * rho_quart[hj]));
      Vd(vertices[i], dofwj) = Vd(vertices[i], dofwj) +
                               ((EltSizes[ej] * Bi[ej]) / 12) *
                                   ((0.5 * rho_mid[ej]) + (rho_quart[hj]));
    }
  }

  //  // Taking finally only the shear related contributions...
  //  il::Array2D<double> VD{Vd.size(0), Vd.size(1) / 2, 0.};
  //
  //  for (il::int_t n = 0; n < VD.size(0); ++n) {
  //
  //    for (il::int_t i = 0, q = 0; i < VD.size(1); ++i, q = q + 2) {
  //
  //      VD(n, i) = Vd(n, q);
  //    }
  //  }

  return Vd;
};
}