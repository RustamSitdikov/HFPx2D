//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <il/Array.h>

#include <il/Array2D.h>

#include "DevLHFMSolver.h"

#include <src/core/ElasticProperties.h>
#include <src/core/Fluid.h>
#include <src/core/Mesh.h>

#include <src/Elasticity/AssemblyDDM.h>
#include <src/Elasticity/Simplified3D.h>
#include <src/FractureFluidFlow/ReynoldsP0.h>
#include <src/core/SimulationParameters.h>

namespace hfp2d {

int TwoParallelHFs(int nelts, double dist) {
  // Test for Reynolds solver
  // 2 fractures parallel separated by dist
  // nelts per frac (numberOfElements should be odd for symmetry)

  // for now for debug we hardcode numberOfElements = 3  so injection is in elt
  // 1 and elt
  // 4

  // step 1 create mesh
  int p = 0;
  double h = 2. / (nelts);  //  element size

  // il::Array<double> x{numberOfElements + 1}; // Not needed
  int Ntot = 2 * nelts;

  il::Array2D<double> xy{Ntot + 2, 2, 0.0};
  il::Array2D<il::int_t> myconn{Ntot, 2, 0};
  il::Array2D<il::int_t> id_DD{Ntot, 2 * (p + 1), 0};
  il::Array2D<il::int_t> id_press{Ntot, p + 1, 0};
  il::Array<il::int_t> fracID{Ntot, 1};
  il::Array<il::int_t> matID{Ntot, 1};
  il::Array<il::int_t> condID{Ntot, 1};  // not needed

  double Ep = 1.;  // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };
  for (il::int_t i = nelts + 1; i < Ntot + 2; ++i) {
    xy(i, 0) = -1. + (i - (nelts + 1)) * h;
    xy(i, 1) = dist;
  };

  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  for (il::int_t i = nelts; i < Ntot; ++i) {
    myconn(i, 0) = i + 1;
    myconn(i, 1) = i + 2;
  };

  for (il::int_t i = 0; i < Ntot; i++) {
    for (il::int_t j = 0; j < 2 * (p + 1); j++) {
      id_DD(i, j) = i * 2 * (p + 1) + j;
    }
  }

  for (il::int_t i = 0; i < Ntot; i++) {
    id_press(i, 0) = i;
  }

  // second fracture
  for (il::int_t i = nelts; i < Ntot; i++) {
    fracID[i] = 2;
  }

  hfp2d::Mesh mesh(xy, myconn, p);

//  hfp2d::Mesh mesh2(xy, myconn, p);

  //    mesh2.AddNTipElements(9,10,2,0.);

  //  il::Array2D<il::int_t> myed =
  //      hfp2d::GetNodalEltConnectivity(mesh.numberOfNodes(),mesh.connectivity());

  //  const il::Array2D<il::int_t> edge = hfp2d::GetEdgesSharing2_old(mesh);

  hfp2d::ElasticProperties myelas(1, 0.);

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
      1000.);  // large pseudo-heigth to reproduce plane-strain kernel

  AddTipCorrectionP0(mesh, myelas, 0, K);
  AddTipCorrectionP0(mesh, myelas, nelts - 1, K);

  AddTipCorrectionP0(mesh, myelas, nelts, K);
  AddTipCorrectionP0(mesh, myelas, Ntot - 1, K);

  //   RE_ARRANGE THE DOF FOR DD HERE !!! this is probably not needed.
  //  K = ReArrangeKP0(mesh, K);

  // initial stress field uniform to test.
  il::Array<double> sig_o{Ntot, 0.1}, tau_o{Ntot, 0.};

  // initial fluid pressure
  il::Array<double> pf_o{
      Ntot, 0.1 + 1.e-1};  // slightly above sig_o to have initial width

  // solve the initial elastic system
  il::Array<double> fini{2 * Ntot, 0.};
  il::int_t j = 0;
  std::cout << "e !!!" << (Ntot + Ntot) << "\n";
  for (il::int_t i = 0; i < (Ntot + Ntot); i = i + 2) {
    fini[i] = tau_o[j];
    fini[i + 1] = -(pf_o[j] - sig_o[j]);
    j++;
  }

  il::int_t ea = 2;

  std::cout << "elt size" << mesh.elt_size(ea) << "w " << h << "\n";

  std::cout << " size of K" << K.size(0) << " by " << K.size(1) << "\n";
  std::cout << " size of f" << fini.size() << "\n";

  il::Status status;
  // use a direct solver
  il::Array<double> dd_ini = il::linearSolve(K, fini, il::io, status);
  status.ok();
  std::cout << "solved \n";
  for (il::int_t i = 0; i < 2 * Ntot; i++) {
    std::cout << " dd " << dd_ini[i] << "\n";
  }

  il::Array<double> width{mesh.numberOfElements(), 0.},
      sheardd{mesh.numberOfElements(), 0.};
  for (il::int_t i = 0; i < mesh.numberOfElements(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
  }

  // create a solution at time t=0 object.
  hfp2d::SolutionAtT Soln =
      hfp2d::SolutionAtT(mesh, 0., width, sheardd, pf_o, sig_o, tau_o);

  // fluid properties.
  hfp2d::Fluid water(1., 0.1, 1.e-10);

  // create source obj. - hardcoded for now....
  il::Array<il::int_t> elt_source{2, 0};
  elt_source[0] = 1;
  elt_source[1] = 4;
  il::Array<double> Qo{2, 0.001};
  hfp2d::Sources the_source = Sources(elt_source, Qo);
  //  std::cout << elt_source[0] << " " << Qo[0] << "\n";

  // create rock properties obj
  il::Array<double> wh_o{1, 1.e-6}, toughness{1, 1.e6}, Carter{1, 0.};
  hfp2d::SolidProperties the_rock =
      SolidProperties(myelas, toughness, wh_o, Carter);

  // call to Reynolds
  double dt = 0.000001;

  hfp2d::SimulationParameters SimulParam;

  bool imp_tip = false;  //
  il::Array<il::int_t> tip_elt{4};
  il::Array<double> tip_width{4, 0.16};

  tip_elt = mesh.tip_elts();
// create a function to get ribbon elements.


  hfp2d::SolutionAtT Soln1 =
      ReynoldsSolverP0(Soln, K, water, the_rock, the_source, dt, imp_tip,
                       tip_elt, tip_width, SimulParam);

  std::cout << "now out of reynolds"
            << "\n";

  return 0;
};
}