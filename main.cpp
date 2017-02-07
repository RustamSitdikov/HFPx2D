//
// HFPx2D project.
//
// Created by Federico Ciardo on 09.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>


// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/math.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "ELHDs.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
#include "Mesh.h"
#include "TimeIncr.h"

// FUNCTION PROTOTYPE
il::Array<double> analytical_solution(double Dp, hfp2d::Mesh mesh);

double_t max_2d(il::Array2D<double> &arr2D);

double_t min_2d(il::Array2D<double> &arr2D);

////////////////////////////////////////////////////////////////////////////////

int main() {

  // Declare the input variables:
  // Ep -> Plain strain modulus
  // Density -> fluid density
  // Cohes -> cohesion (for M-C criterion)
  // CompressFluid -> compressibility of the fluid
  // Visc -> fluid viscosity
  // Init_dil -> initial value of dilatancy
  // Incr_dil -> Increment of dilatancy (difference between residual/peak
  //             dilatancy and initial dilatancy value)
  // d_wf -> slip dw for scaling (see exponential law in the report)
  // d_wd ->  slip dw for scaling (see dilatancy law in the report)
  // Peak_fric -> peak friction coefficient
  // Resid_fric -> residual friction coefficient
  // sigma_n0 -> Ambient normal stress
  // P0 -> Ambient pore pressure in the crack/fault
  // Source -> source term
  // l -> half crack/fault length
  double Ep, Density, Cohes, CompressFluid, Visc, Init_dil, Incr_dil, d_wf,
      d_wd, Peak_fric, Resid_fric, sigma_n0, P0, Source, l;
  int Nnodes;

  // Read the file that contains the input data
  std::ifstream read_file(
      "/Users/federicociardo/ClionProjects/HFPx2D/InputData.txt");

  // Check whether the input data file is opened or not
  IL_ASSERT(read_file.is_open());
  // Read the data
  read_file >> Nnodes >> Ep >> Density >> Cohes >> CompressFluid >> Visc >>
      Init_dil >> Incr_dil >> d_wf >> d_wd >> Peak_fric >> Resid_fric >>
      sigma_n0 >> P0 >> Source >> l;
  // Close the data file
  read_file.close();

  // Declare other variables
  // Number of elements
  int Nelts = Nnodes - 1;
  // Number of collocation points
  int NcollPoints = 2 * Nelts;
  // Degrees of freedom per node
  int dof_dim = 2;

  // Create the vector of density (rho), initial stress state (uniform for now,
  // Sigma0), initial/ambient pore pressure (Amb_press)

  // Fluid density matrix {{rho_1left, rho_2right},{rho_2left,rho_2right} ..}
  // Remember: continuous linear variation
  // Size -> Neltsx2
  il::Array2D<double> rho{Nelts,2, 0.};
  for (il::int_t i{0}; i < rho.size(0); ++i) {
      for (il::int_t j = 0; j < rho.size(1); ++j) {

          rho(i,j) = Density;
      }
  }

  // Declare some variables -> sigma_s0 is the ambient shear stress
  // (which is related to the ambient normal stress)
  // Dp is the constant overpressure
  double sigma_s0 = 0.55 * sigma_n0;
  double Dp = 0.5 * (sigma_n0 - P0);

  // Matrix -> {Ambient shear stress , Ambient normal stress}
  // for each collocation points
  il::Array2D<double> Sigma0{NcollPoints, 2, 0.};

  for (il::int_t i{0}; i < Sigma0.size(0); ++i) {
    Sigma0(i, 0) = sigma_s0;
    Sigma0(i, 1) = sigma_n0;
  }

  // Ambient pore pressure at nodal points {P0_1, P0_2, P0_3, ..}
  // Remember -> Pressure varies linearly and continuously over the elements
  // Size -> Nnodes
  il::Array<double> Amb_press{Nnodes, 0.};
  for (il::int_t k{0}; k < Amb_press.size(); ++k) {

    Amb_press[k] = P0;
  }

  // Interpolation order
  int p = 1;
  // element size -> for now UNIFORM mesh
  double h = (2 * l) / Nelts;

  /// Mesh ///
  // create a basic 1D mesh

  // xy -> matrix of mesh coordinates {{x1,y1},{x2,y2},{x3,y3}..}
  il::Array2D<double> xy{Nnodes, 2, 0.};
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -l + (i * h);
    xy(i, 1) = 0.;
  }

  // Connectivity matrix
  il::Array2D<int> myconn{Nelts, 2, 0.};
  for (il::int_t i{0}; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  }

  // Create mesh object
  hfp2d::Mesh mesh;
  mesh.set_values(xy, myconn);

  // Initial pore pressure profile (which must be added to the
  // ambient pressure in the crack/fault)
  // Size -> Nnodes
  il::Array<double> Pinit{Nnodes, 0.};
  Pinit = analytical_solution(Dp, mesh);

  // Find the position of the max value in Pinit
  il::int_t idx;
  idx = hfp2d::find(Pinit, hfp2d::max_1d(Pinit));

  // Source vector whose size is Nnodes
  // for constant overpressure -> source term is equal to zero
  il::Array<double> S{Nnodes, 0.};
  S[idx] = Source;

  // Total number of dofs
  int Ndof = (Nnodes - 1) * (p + 1) * 2;

  // Elasticity matrix assembling
  il::Array2D<double> kmat{Ndof, Ndof, 0.};
  hfp2d::basic_assembly(kmat, mesh, hfp2d::dofhandle_dg_full2d(dof_dim, Nelts, p), p,
                        Ep);

  /// Solution of fluid injection into frictional weakening dilatant fault ///
  hfp2d::time_incr(mesh, p, Cohes, kmat, Incr_dil, d_wd, rho, Init_dil,
                   CompressFluid, Visc, S, dof_dim, Peak_fric, Resid_fric, d_wf,
                   Sigma0, Amb_press, Pinit);


  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution for pore pressure distribution in a fracture/fault with
// constant permeability subjected to a constant overpressure Dp
// It return a vector whose size is Nnodes {Pinit_1, Pinit_2, Pinit_3 ...}
// Remember -> pressure varies linearly and continuously over the elements

il::Array<double> analytical_solution(double Dp, hfp2d::Mesh mesh) {

  // Inputs:
  //  - Dp -> Constant overpressure
  //  - mesh -> mesh class
  //  - myconn -> connectivity matrix

  il::Array<double> Pinit{mesh.nelts() + 1, 0.};

  for (il::int_t j = 0; j < Pinit.size(); ++j) {

    Pinit[j] = Dp * (erfc(fabs(mesh.Coor(j, 0) / sqrt(10 * 0.005)) / 2));
  }

  return Pinit;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the max value in an array2D.
// arr2D -> array2D in which we want to find the max

double_t max_2d(il::Array2D<double> &arr2D) {

  double_t max;
  max = arr2D(0., 0.);

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (max < arr2D(i, j)) {

        max = arr2D(i, j);
      }
    }
  }

  return max;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the min value in an array2D.
// arr2D -> array2D in which we want to find the min

double_t min_2d(il::Array2D<double> &arr2D) {

  double_t min;
  min = arr2D(0., 0.);

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (min > arr2D(i, j)) {

        min = arr2D(i, j);
      }
    }
  }

  return min;
}

