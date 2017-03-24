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
#include <sys/stat.h>
#include <sys/types.h>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/Timer.h>
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
il::Array<double> analytical_solution(double Dp, double alpha,
                                      hfp2d::time_parameters time_parameters,
                                      hfp2d::Mesh mesh, il::io_t);

////////////////////////////////////////////////////////////////////////////////

int main() {

  ///  **** Read the input data from TOML input file **** ///

  il::String filename =
      "/Users/federicociardo/ClionProjects/HFPx2D-Collscheme/InputData.toml";

  il::Status status{};
  auto config =
      il::load<il::HashMap<il::String, il::Dynamic>>(filename, il::io, status);
  status.abort_on_error();

  il::int_t i;

  il::int_t Nnodes = 0;
  i = config.search("Number_of_nodes");
  if (config.found(i) && config.value(i).is_integer()) {
    Nnodes = config.value(i).to_integer();
  } else {
    Nnodes = 0;
  }

  double Ep;
  i = config.search("Plain_strain_modulus");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Ep = config.value(i).to_floating_point();
  } else {
    Ep = 0;
  }

  double Density;
  i = config.search("Density");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Density = config.value(i).to_floating_point();
  } else {
    Density = 0;
  }

  double Cohes;
  i = config.search("Cohesion");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Cohes = config.value(i).to_floating_point();
  } else {
    Cohes = 0;
  }

  double CompressFluid;
  i = config.search("Fluid_compressibility");
  if (config.found(i) && config.value(i).is_floating_point()) {
    CompressFluid = config.value(i).to_floating_point();
  } else {
    CompressFluid = 0;
  }

  double Visc;
  i = config.search("Fluid_viscosity");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Visc = config.value(i).to_floating_point();
  } else {
    Visc = 0;
  }

  double Init_hydr_width;
  i = config.search("Initial_hydraulic_width");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Init_hydr_width = config.value(i).to_floating_point();
  } else {
    Init_hydr_width = 0;
  }

  double Incr_dil;
  i = config.search("Increment_of_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Incr_dil = config.value(i).to_floating_point();
  } else {
    Incr_dil = 0;
  }

  double d_wfriction_layer1;
  i = config.search("Slip_dw_for_friction_layer1");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wfriction_layer1 = config.value(i).to_floating_point();
  } else {
    d_wfriction_layer1 = 0;
  }

  double d_wfriction_layer2;
  i = config.search("Slip_dw_for_friction_layer2");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wfriction_layer2 = config.value(i).to_floating_point();
  } else {
    d_wfriction_layer2 = 0;
  }

  double d_wfriction_layer3;
  i = config.search("Slip_dw_for_friction_layer3");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wfriction_layer3 = config.value(i).to_floating_point();
  } else {
    d_wfriction_layer3 = 0;
  }

  double d_wdilatancy;
  i = config.search("Slip_dw_for_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wdilatancy = config.value(i).to_floating_point();
  } else {
    d_wdilatancy = 0;
  }

  double Peak_fric_layer1;
  i = config.search("Peak_friction_coefficient_layer1");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Peak_fric_layer1 = config.value(i).to_floating_point();
  } else {
    Peak_fric_layer1 = 0;
  }

  double Resid_fric_layer1;
  i = config.search("Residual_friction_coefficient_layer1");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Resid_fric_layer1 = config.value(i).to_floating_point();
  } else {
    Resid_fric_layer1 = 0;
  }

  double Peak_fric_layer2;
  i = config.search("Peak_friction_coefficient_layer2");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Peak_fric_layer2 = config.value(i).to_floating_point();
  } else {
    Peak_fric_layer2 = 0;
  }

  double Resid_fric_layer2;
  i = config.search("Residual_friction_coefficient_layer2");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Resid_fric_layer2 = config.value(i).to_floating_point();
  } else {
    Resid_fric_layer2 = 0;
  }

  double Peak_fric_layer3;
  i = config.search("Peak_friction_coefficient_layer3");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Peak_fric_layer3 = config.value(i).to_floating_point();
  } else {
    Peak_fric_layer3 = 0;
  }

  double Resid_fric_layer3;
  i = config.search("Residual_friction_coefficient_layer3");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Resid_fric_layer3 = config.value(i).to_floating_point();
  } else {
    Resid_fric_layer3 = 0;
  }

  double sigma_n0;
  i = config.search("Ambient_normal_stress");
  if (config.found(i) && config.value(i).is_floating_point()) {
    sigma_n0 = config.value(i).to_floating_point();
  } else {
    sigma_n0 = 0;
  }

  double sigma_s0;
  i = config.search("Ambient_shear_stress");
  if (config.found(i) && config.value(i).is_floating_point()) {
    sigma_s0 = config.value(i).to_floating_point();
  } else {
    sigma_s0 = 0;
  }

  double P0;
  i = config.search("Ambient pore pressure");
  if (config.found(i) && config.value(i).is_floating_point()) {
    P0 = config.value(i).to_floating_point();
  } else {
    P0 = 0;
  }

  double Dp;
  i = config.search("Constant_overpressure");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Dp = config.value(i).to_floating_point();
  } else {
    Dp = 0;
  }

  double Source;
  i = config.search("Source");
  if (config.found(i) && config.value(i).is_floating_point()) {
    Source = config.value(i).to_floating_point();
  } else {
    Source = 0;
  }

  double l;
  i = config.search("Half_fault_length");
  if (config.found(i) && config.value(i).is_floating_point()) {
    l = config.value(i).to_floating_point();
  } else {
    l = 0;
  }

  double t_max;
  i = config.search("t_max");
  if (config.found(i) && config.value(i).is_floating_point()) {
    t_max = config.value(i).to_floating_point();
  } else {
    t_max = 0;
  }

  double t_0plus;
  i = config.search("t_0plus");
  if (config.found(i) && config.value(i).is_floating_point()) {
    t_0plus = config.value(i).to_floating_point();
  } else {
    t_0plus = 0;
  }

  double TimeStep_max;
  i = config.search("TimeStep_max");
  if (config.found(i) && config.value(i).is_floating_point()) {
    TimeStep_max = config.value(i).to_floating_point();
  } else {
    TimeStep_max = 0;
  }

  double TimeStep_min;
  i = config.search("TimeStep_min");
  if (config.found(i) && config.value(i).is_floating_point()) {
    TimeStep_min = config.value(i).to_floating_point();
  } else {
    TimeStep_min = 0;
  }

  // Set the structure members of friction
  hfp2d::Parameters_friction fric_parameters;
  fric_parameters.Peak_fric_coeff_layer1 = Peak_fric_layer1;
  fric_parameters.Resid_fric_coeff_layer1 = Resid_fric_layer1;
  fric_parameters.d_wf_layer1 = d_wfriction_layer1;
  fric_parameters.Peak_fric_coeff_layer2 = Peak_fric_layer2;
  fric_parameters.Resid_fric_coeff_layer2 = Resid_fric_layer2;
  fric_parameters.d_wf_layer2 = d_wfriction_layer2;
  fric_parameters.Peak_fric_coeff_layer3 = Peak_fric_layer3;
  fric_parameters.Resid_fric_coeff_layer3 = Resid_fric_layer3;
  fric_parameters.d_wf_layer3 = d_wfriction_layer3;

  // Set the structure members of dilatancy
  hfp2d::Parameters_dilatancy dilat_parameters;
  dilat_parameters.Init_hydr_width = Init_hydr_width;
  dilat_parameters.Incr_dil = Incr_dil;
  dilat_parameters.d_wd = d_wdilatancy;

  // Set the structure members for time parameters
  hfp2d::time_parameters time_parameters;
  time_parameters.t_0plus = t_0plus;
  time_parameters.t_max = t_max;
  time_parameters.TimeStep_max = TimeStep_max;
  time_parameters.TimeStep_min = TimeStep_min;

  // Declare other variables:
  // Number of elements
  il::int_t Nelts = Nnodes - 1;
  // Number of collocation points
  il::int_t NCollPoints = 2 * Nelts;
  // Degrees of freedom per node
  int dof_dim = 2;

  // Creating a directory for output results
  // Remember always to delete the existing directory (if already created)!
  std::string Directory_results{"/Users/federicociardo/ClionProjects/"
                                "HFPx2D-Collscheme/Results/"
                                "Test1su100_NoDil_055(2)/"};

  if (mkdir(Directory_results.c_str(), 0777) == -1) {
    std::cerr << "Error in creating the output directory:  " << strerror(errno)
              << std::endl;
    exit(1);
  }

  /* Create:
   *  - the matrix of density (rho)
   *  - initial stress state Sigma0 (uniform for now)
   *  - initial/ambient pore pressure (Amb_press)
   *  - the matrix of cohesion for each collocation point */

  // Fluid density matrix {{rho_1left, rho_1right},{rho_2left,rho_2right} ..}
  // Remember: continuous linear variation -> rho_2right = rho_1left and so on..
  il::Array2D<double> rho{Nelts, 2, Density};

  // Set the structure members of fluid
  hfp2d::Parameters_fluid fluid_parameters;
  fluid_parameters.compressibility = CompressFluid;
  fluid_parameters.density = rho;
  fluid_parameters.viscosity = Visc;

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

  // Element size -> for now UNIFORM mesh, but later it might be non-uniform
  double h = (2 * l) / Nelts;

  /// ****** Mesh ****** ///
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
  // alpha -> initial fault/fracture diffusivity [Lˆ2/T]
  // t_0plus -> initial time (starting time)
  double alpha = 10;
  il::Array<double> Pinit{Nnodes, 0};
  Pinit = analytical_solution(Dp, alpha, time_parameters, mesh, il::io);

  // Find the injection point (position of the max value in Pinit)
  int inj_point;
  inj_point = hfp2d::find(Pinit, hfp2d::max_1d(Pinit, il::io), il::io);

  // Source vector whose size is Nnodes
  // (for constant overpressure -> source term is equal to zero)
  il::Array<double> S{Nnodes, 0};
  S[inj_point] = Source;

  // Total number of dofs
  int Ndof = (Nnodes - 1) * (p + 1) * 2;

  // Matrix to switch from nodal points to collocation points
  il::Array2D<double> Fetc{4 * mesh.nelts(), mesh.nelts() + 1, 0};
  Fetc = hfp2d::from_edge_to_col_cg(
      dof_dim, hfp2d::dofhandle_dg_full2d(dof_dim, mesh.nelts(), p, il::io),
      hfp2d::dofhandle_cg2d(dof_dim, mesh.nelts(), il::io), il::io);

  // Get the elasticity matrix
  il::Array2D<double> kmat{Ndof, Ndof, 0};
  hfp2d::basic_assembly(
      kmat, mesh, hfp2d::dofhandle_dg_full2d(dof_dim, Nelts, p, il::io), p, Ep);

  /// Solution of fluid injection into frictional weakening dilatant fault ///
  hfp2d::time_incr(inj_point, NCollPoints, mesh, p, cohes, kmat,
                   fric_parameters, dilat_parameters, fluid_parameters, S,
                   dof_dim, Sigma0, Amb_press, Pinit, Directory_results, XColl,
                   Fetc, h, time_parameters, il::io);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution for evolving pore pressure distribution in a
// fracture/fault with constant permeability subjected to a constant
// overpressure Dp.
// It returns a vector whose size is Nnodes {Pinit_1, Pinit_2, Pinit_3 ...}. It
// represents the initial pore pressure distribution/perturbation at t_0plus
// needed to activate the shear crack.
// Remember -> pressure varies linearly and continuously over the elements
////////////////////////////////////////////////////////////////////////////////

il::Array<double> analytical_solution(double Dp, double alpha,
                                      hfp2d::time_parameters time_parameters,
                                      hfp2d::Mesh mesh, il::io_t) {

  // Inputs:
  //  - Dp -> Constant overpressure in the middle of the fracture/fault
  //  - alpha -> Initial fault/fracture diffusivity (Lˆ2/T)
  //  - hfp2d::time_parameters time_parameters -> structure that contains all
  //    the time parameters
  //  - mesh -> mesh class
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> Pinit{mesh.nelts() + 1, 0};

  for (il::int_t j = 0; j < Pinit.size(); ++j) {
    Pinit[j] =
        Dp *
        (erfc(fabs((mesh.node(j, 0)) / sqrt(alpha * time_parameters.t_0plus)) /
              2));
  }

  return Pinit;
}
