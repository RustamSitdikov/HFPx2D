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
//#include <string>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
//#include <il/math.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "ELHDs.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
//#include "Mesh.h"
#include "TimeIncr.h"
#include <il/io/yaml.h>

// FUNCTION PROTOTYPE
il::Array<double> analytical_solution(double Dp, hfp2d::Mesh mesh, il::io_t);

double_t max_2d(const il::Array2D<double> &arr2D, il::io_t);

double_t min_2d(const il::Array2D<double> &arr2D, il::io_t);

////////////////////////////////////////////////////////////////////////////////

int main() {

  il::Status status{};
  il::String filename =
      "/Users/federicociardo/ClionProjects/HFPx2D/InputData.yaml";
  il::Yaml config = il::load<il::Yaml>(filename, il::io, status);
  status.abort_on_error();

  il::int_t i;

  il::int_t Nnodes = 0;
  i = config.search("Number of nodes");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::integer);
    Nnodes = config.value_integer(i);
  }

  double Ep;
  i = config.search("Plain strain modulus");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Ep = config.value_floating_point(i);
  }

  double Density;
  i = config.search("Density");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Density = config.value_floating_point(i);
  }

  double Cohes;
  i = config.search("Cohesion");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Cohes = config.value_floating_point(i);
  }

  double CompressFluid;
  i = config.search("Fluid compressibility");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    CompressFluid = config.value_floating_point(i);
  }

  double Visc;
  i = config.search("Fluid viscosity");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Visc = config.value_floating_point(i);
  }

  double Init_dil;
  i = config.search("Initial value of dilatancy");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Init_dil = config.value_floating_point(i);
  }

  double Incr_dil;
  i = config.search("Increment of dilatancy");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Incr_dil = config.value_floating_point(i);
  }

  double d_wf;
  i = config.search("Slip dw for friction");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    d_wf = config.value_floating_point(i);
  }

  double d_wd;
  i = config.search("Slip dw for dilatancy");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    d_wd = config.value_floating_point(i);
  }

  double Peak_fric;
  i = config.search("Peak friction coefficient");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Peak_fric = config.value_floating_point(i);
  }

  double Resid_fric;
  i = config.search("Residual friction coefficient");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Resid_fric = config.value_floating_point(i);
  }

  double sigma_n0;
  i = config.search("Ambient normal stress");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    sigma_n0 = config.value_floating_point(i);
  }

  double P0;
  i = config.search("Ambient pore pressure");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    P0 = config.value_floating_point(i);
  }

  double Source;
  i = config.search("Source");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    Source = config.value_floating_point(i);
  }

  double l;
  i = config.search("Half fault length");
  if (config.found(i)) {
    IL_EXPECT_FAST(config.type(i) == il::Type::floating_point);
    l = config.value_floating_point(i);
  }

  // Declare other variables
  // Number of elements
  il::int_t Nelts = Nnodes - 1;
  // Number of collocation points
  il::int_t NcollPoints = 2 * Nelts;
  // Degrees of freedom per node
  int dof_dim = 2;

  // Directory for output results
  std::string Directory_results{
      "/Users/federicociardo/ClionProjects/HFPx2D/Results/"};

  // Create the matrix of density (rho), initial stress state (uniform for now,
  // Sigma0), initial/ambient pore pressure (Amb_press)

  // Fluid density matrix {{rho_1left, rho_2right},{rho_2left,rho_2right} ..}
  // Remember: continuous linear variation
  // Size -> Neltsx2
  il::Array2D<double> rho{Nelts, 2, 0};
  for (il::int_t ii = 0; ii < rho.size(0); ++ii) {
    for (il::int_t j = 0; j < rho.size(1); ++j) {

      rho(ii, j) = Density;
    }
  }

  // Declare some variables -> sigma_s0 is the ambient shear stress
  // (which is related to the ambient normal stress)
  // Dp is the constant overpressure
  double sigma_s0 = 0.7 * sigma_n0;
  double Dp = 0.5 * (sigma_n0 - P0);

  // Matrix -> {Ambient shear stress , Ambient normal stress}
  // for each collocation points
  il::Array2D<double> Sigma0{NcollPoints, 2, 0};

  for (il::int_t i3 = 0; i3 < Sigma0.size(0); ++i3) {
    Sigma0(i3, 0) = sigma_s0;
    Sigma0(i3, 1) = sigma_n0;
  }

  // Ambient pore pressure at nodal points {P0_1, P0_2, P0_3, ..}
  // Remember -> Pressure varies linearly and continuously over the elements
  // Size -> Nnodes
  il::Array<double> Amb_press{Nnodes, 0};
  for (il::int_t k = 0; k < Amb_press.size(); ++k) {

    Amb_press[k] = P0;
  }

  // Interpolation order
  int p = 1;

  // Element size -> for now UNIFORM mesh, but later it might not be uniform
  double h = (2 * l) / Nelts;

  /// Mesh ///
  // create a basic 1D mesh

  // xy -> matrix of mesh coordinates {{x1,y1},{x2,y2},{x3,y3}..}
  il::Array2D<double> xy{Nnodes, 2, 0};
  for (il::int_t i2 = 0; i2 < xy.size(0); ++i2) {
    xy(i2, 0) = -l + (i2 * h);
    xy(i2, 1) = 0.;
  }

  // Connectivity matrix
  il::Array2D<int> myconn{Nelts, 2, 0};
  for (il::int_t i4 = 0; i4 < myconn.size(0); ++i4) {
    myconn(i4, 0) = i4;
    myconn(i4, 1) = i4 + 1;
  }

  // Create mesh object
  hfp2d::Mesh mesh;
  mesh.set_values(xy, myconn);

  // Get collocation points' information (needed for calculation of slippage
  // length)
  il::Array<double> XColl{2 * mesh.nelts(), 0};
  il::Array2D<double> Xs{2, 2, 0};
  hfp2d::SegmentCharacteristic coo_coll;

  for (il::int_t m3 = 0, m4 = 0; m3 < mesh.nelts(); ++m3, m4 = m4 + 2) {
    for (il::int_t i = 0; i < Xs.size(1); ++i) {
      Xs(0, i) = mesh.Coor(
          hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io)(m3, 0), i);
      Xs(1, i) = mesh.Coor(
          hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io)(m3, 1), i);
    }

    coo_coll = hfp2d::get_segment_DD_characteristic(Xs, p);
    XColl[m4] = coo_coll.CollocationPoints(0, 0);
    XColl[m4 + 1] = coo_coll.CollocationPoints(1, 0);
  }

  // Initial pore pressure profile (which must be added to the
  // ambient pressure in the crack/fault)
  // Size -> Nnodes
  il::Array<double> Pinit{Nnodes, 0};
  Pinit = analytical_solution(Dp, mesh, il::io);

  // Find the position of the max value in Pinit
  int idx;
  idx = hfp2d::find(Pinit, hfp2d::max_1d(Pinit, il::io), il::io);

  // Source vector whose size is Nnodes
  // for constant overpressure -> source term is equal to zero
  il::Array<double> S{Nnodes, 0};
  S[idx] = Source;

  // Total number of dofs
  int Ndof = (Nnodes - 1) * (p + 1) * 2;

  // Elasticity matrix assembling
  il::Array2D<double> kmat{Ndof, Ndof, 0};
  hfp2d::basic_assembly(
      kmat, mesh, hfp2d::dofhandle_dg_full2d(dof_dim, Nelts, p, il::io), p, Ep);

  /// Solution of fluid injection into frictional weakening dilatant fault ///
  hfp2d::time_incr(mesh, p, Cohes, kmat, Incr_dil, d_wd, rho, Init_dil,
                   CompressFluid, Visc, S, dof_dim, Peak_fric, Resid_fric, d_wf,
                   Sigma0, Amb_press, Pinit, Directory_results, XColl, il::io);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution for pore pressure distribution in a fracture/fault with
// constant permeability subjected to a constant overpressure Dp
// It returns a vector whose size is Nnodes {Pinit_1, Pinit_2, Pinit_3 ...}
// Remember -> pressure varies linearly and continuously over the elements

il::Array<double> analytical_solution(double Dp, hfp2d::Mesh mesh, il::io_t) {

  // Inputs:
  //  - Dp -> Constant overpressure in the middle of the fracture/fault
  //  - mesh -> mesh class
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> Pinit{mesh.nelts() + 1, 0};

  for (il::int_t j = 0; j < Pinit.size(); ++j) {

    Pinit[j] = Dp * (erfc(fabs(mesh.Coor(j, 0) / sqrt(10 * 0.0005)) / 2));
  }

  return Pinit;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the max value in an array2D.
// arr2D -> array2D in which we want to find the max

double_t max_2d(const il::Array2D<double> &arr2D, il::io_t) {

  double_t max;
  max = arr2D(0, 0);

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

double_t min_2d(const il::Array2D<double> &arr2D, il::io_t) {

  double_t min;
  min = arr2D(0, 0);

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (min > arr2D(i, j)) {

        min = arr2D(i, j);
      }
    }
  }

  return min;
}
