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

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>

// Inclusion from the project
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "Dilatancy.h"
#include "FVM.h"
#include "Friction.h"
#include "FromEdgeToCol.h"
#include "MC_criterion.h"
#include "TimeIncr.h"
#include <il/Toml.h>

// FUNCTION PROTOTYPE
il::Array<double> analytical_solution(double Dp, double alpha, double t_0plus,
                                      hfp2d::Mesh mesh, il::io_t);

////////////////////////////////////////////////////////////////////////////////

int main() {

  /*  **** Read the input data fromm TOML input file **** */

  il::String filename =
      "/Users/federicociardo/ClionProjects/HFPx2D/InputData.toml";

  il::Status status{};
  auto config =
      il::load<il::HashMap<il::String, il::Dynamic>>(filename, il::io, status);
  status.abort_on_error();

  il::int_t i;

  il::int_t Nnodes = 0;
  i = config.search("Number_of_nodes");
  if (config.found(i) && config.value(i).is_integer()) {
    Nnodes = config.value(i).to_integer();
  }

  double Ep;
  i = config.search("Plain_strain_modulus");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Ep = config.value(i).to_floating_point();
  }

  double Density;
  i = config.search("Density");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Density = config.value(i).to_floating_point();
  }

  double Cohes;
  i = config.search("Cohesion");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Cohes = config.value(i).to_floating_point();
  }

  double CompressFluid;
  i = config.search("Fluid_compressibility");
  if (config.found(i) && config.value(i).is_floating_point()) {
    CompressFluid = config.value(i).to_floating_point();
  }

  double Visc;
  i = config.search("Fluid_viscosity");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Visc = config.value(i).to_floating_point();
  }

  double Init_dil;
  i = config.search("Initial_value_of_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Init_dil = config.value(i).to_floating_point();
  }

  double Incr_dil;
  i = config.search("Increment_of_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Incr_dil = config.value(i).to_floating_point();
  }

  double d_wfriction;
  i = config.search("Slip_dw_for_friction");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wfriction = config.value(i).to_floating_point();
  }

  double d_wdilatancy;
  i = config.search("Slip_dw_for_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wdilatancy = config.value(i).to_floating_point();
  }

  double Peak_fric;
  i = config.search("Peak_friction_coefficient");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Peak_fric = config.value(i).to_floating_point();
  }

  double Resid_fric;
  i = config.search("Residual_friction_coefficient");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Resid_fric = config.value(i).to_floating_point();
  }

  double sigma_n0;
  i = config.search("Ambient_normal_stress");
  if (config.found(i) && config.value(i).is_floating_point()) {
    sigma_n0 = config.value(i).to_floating_point();
  }

  double P0;
  i = config.search("Ambient pore pressure");
  if (config.found(i) && config.value(i).is_floating_point()) {
    P0 = config.value(i).to_floating_point();
  }

  double Source;
  i = config.search("Source");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Source = config.value(i).to_floating_point();
  }

  double l;
  i = config.search("Half_fault_length");
  if (config.found(i) && config.value(i).is_floating_point()) {
    l = config.value(i).to_floating_point();
  }

  // Set the structure members of friction
  hfp2d::Parameters_friction fric_parameters;
  fric_parameters.Peak_fric_coeff = Peak_fric;
  fric_parameters.Resid_fric_coeff = Resid_fric;
  fric_parameters.d_wf = d_wfriction;

  // Set the structure members of dilatancy
  hfp2d::Parameters_dilatancy dilat_parameters;
  dilat_parameters.Init_dil = Init_dil;
  dilat_parameters.Incr_dil = Incr_dil;
  dilat_parameters.d_wd = d_wdilatancy;

  // Declare other variables:
  // Number of elements
  il::int_t Nelts = Nnodes - 1;
  // Number of collocation points
  il::int_t NCollPoints = 2 * Nelts;
  // Degrees of freedom per node
  int dof_dim = 2;

  // Directory for output results
  std::string Directory_results{
      "/Users/federicociardo/ClionProjects/HFPx2D/Results/"};

  /* Create:
   *  - the matrix of density (rho)
   *  - initial stress state (uniform for now, Sigma0)
   *  - initial/ambient pore pressure (Amb_press)
   *  - the matrix of cohesion for each collocation point */

  // Fluid density matrix {{rho_1left, rho_1right},{rho_2left,rho_2right} ..}
  // Remember: continuous linear variation -> rho_2right = rho_2left and so on..
  il::Array2D<double> rho{Nelts, 2, 0};
  for (il::int_t ii = 0; ii < rho.size(0); ++ii) {
    for (il::int_t j = 0; j < rho.size(1); ++j) {
      rho(ii, j) = Density;
    }
  }

  // Set the structure members of fluid
  hfp2d::Parameters_fluid fluid_parameters;
  fluid_parameters.compressibility = CompressFluid;
  fluid_parameters.density = rho;
  fluid_parameters.viscosity = Visc;

  // Declare some variables:
  // sigma_s0 is the ambient shear stress (which is related to the ambient
  // normal stress)
  // Dp is the constant overpressure in the middle of the fault/fracture
  double sigma_s0 = 0.55 * sigma_n0;
  double Dp = 0.5 * (sigma_n0 - P0);

  // Matrix of initial stress state
  // {Ambient SHEAR stress , Ambient NORMAL stress} for each collocation points
  il::Array2D<double> Sigma0{NCollPoints, 2, 0};
  for (il::int_t i3 = 0; i3 < Sigma0.size(0); ++i3) {
    Sigma0(i3, 0) = sigma_s0;
    Sigma0(i3, 1) = sigma_n0;
  }

  // Ambient pore pressure at nodal points {P0_1, P0_2, P0_3, ..}
  // Remember -> Pressure varies linearly and continuously over the elements
  il::Array<double> Amb_press{Nnodes, 0};
  for (il::int_t k = 0; k < Amb_press.size(); ++k) {
    Amb_press[k] = P0;
  }

  // Matrix of cohesion at collocation points
  il::Array<double> cohes{NCollPoints, 0};
  for (il::int_t n = 0; n < cohes.size(); ++n) {
    cohes[n] = Cohes;
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

  // Get collocation points' information
  il::Array<double> XColl{2 * mesh.nelts(), 0};
  hfp2d::SegmentCharacteristic coo_coll;
  for (il::int_t m = 0, m4 = 0; m < mesh.nelts(); ++m, m4 = m4 + 2) {
    coo_coll = hfp2d::get_segment_DD_characteristic(mesh, m, p);
    XColl[m4] = coo_coll.CollocationPoints(0, 0);
    XColl[m4 + 1] = coo_coll.CollocationPoints(1, 0);
  }

  // Initial pore pressure profile needed to activate the shear crack (which
  // must be added to the ambient pressure in the crack/fault)
  // alpha -> initial fault/fracture diffusivity
  // t_0plus -> initial time (starting time)
  double alpha = 10;
  double t_0plus = 0.005;
  il::Array<double> Pinit{Nnodes, 0};
  Pinit = analytical_solution(Dp, alpha, t_0plus, mesh, il::io);

  // Find the injection point (position of the max value in Pinit)
  int inj_point;
  inj_point = hfp2d::find(Pinit, hfp2d::max_1d(Pinit, il::io), il::io);

  // Source vector whose size is Nnodes
  // (for constant overpressure, source term is equal to zero)
  il::Array<double> S{Nnodes, 0};
  S[inj_point] = Source;

  // Total number of dofs
  int Ndof = (Nnodes - 1) * (p + 1) * 2;

  // Initialization of  matrix to switch from nodal points to collocation points
  il::Array2D<double> Fetc{4 * mesh.nelts(), mesh.nelts() + 1, 0};
  Fetc = hfp2d::from_edge_to_col_cg(
      dof_dim, hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io),
      hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io), il::io);

  // Elasticity matrix assembling
  il::Array2D<double> kmat{Ndof, Ndof, 0};
  hfp2d::basic_assembly(
      kmat, mesh, hfp2d::dofhandle_dg_full2d(dof_dim, Nelts, p, il::io), p, Ep);

  /// Solution of fluid injection into frictional weakening dilatant fault ///
  hfp2d::time_incr(t_0plus, inj_point, NCollPoints, mesh, p, cohes, kmat,
                   fric_parameters, dilat_parameters, fluid_parameters, S,
                   dof_dim, Sigma0, Amb_press, Pinit, Directory_results, XColl,
                   Fetc, il::io);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution for evolving pore pressure distribution in a
// fracture/fault with constant permeability subjected to a constant
// overpressure Dp.
// It returns a vector whose size is Nnodes {Pinit_1, Pinit_2, Pinit_3 ...}. It
// represents the initial pore pressure distribution/perturbation needed to
// activate the shear crack.
// Remember -> pressure varies linearly and continuously over the elements

il::Array<double> analytical_solution(double Dp, double alpha, double t_0plus,
                                      hfp2d::Mesh mesh, il::io_t) {

  // Inputs:
  //  - Dp -> Constant overpressure in the middle of the fracture/fault
  //  - alpha -> Initial fault/fracture diffusivity (Lˆ2 T^-1)
  //  - t_0plus -> Initial time (starting time)
  //  - mesh -> mesh class
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> Pinit{mesh.nelts() + 1, 0};

  for (il::int_t j = 0; j < Pinit.size(); ++j) {
    Pinit[j] = Dp * (erfc(fabs((mesh.node(j, 0)) / sqrt(alpha * t_0plus)) / 2));
  }

  return Pinit;
}
