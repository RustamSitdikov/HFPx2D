//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 15.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <fstream>
#include <iostream>

#include <il/base.h>
#include <il/norm.h>

#include <src/core/SimulationParameters.h>
#include <src/solvers/MultiFracsSolver.h>
#include <src/wellbore/LoadInputMultistage.h>
#include <src/wellbore/WellFlowP0.h>
#include <src/wellbore/WellSolution.h>

#include <src/solvers/MultiFracsSolution.h>

#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/Simplified3D.h>

#include <src/util/json.hpp>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
using json = nlohmann::json;
// write a routine - for well + n fracs benchmark
// with json inputs.

int MultipleFracsPropagation() {
  // routine for the propagation of multiple Heigth contained HFs from a
  // horizontal wellbore
  // DEBUGGING

  //

  std::string wellfilename = "../Debug/WellTest.json";

  std::ifstream input(wellfilename);  // ?
  json js;
  input >> js;

  if ((js.count("Wellbore mesh") != 1)) {
    std::cout << "No wellbore mesh in json input file ";
    il::abort();
  }
  json j_wmesh = js["Wellbore mesh"];

  if ((js.count("Model parameters") != 1)) {
    std::cout << "No parameters in input file ";
    il::abort();
  }
  json j_params = js["Model parameters"];

  if ((js.count("Simulation parameters") != 1)) {
    std::cout << "No Simulation parameters in input file ";
    il::abort();
  }
  json j_simul = js["Simulation parameters"];

  // start loading in our objects.
  hfp2d::WellMesh the_well = loadWellMesh(j_wmesh);
  hfp2d::WellInjection w_inj = loadWellParameters(j_params, the_well);

  if (j_params.count("Fluid properties") != 1) {
    std::cout << "No fluid properties input in  model parameters";
    il::abort();
  }

  json j_fluid = j_params["Fluid properties"];
  hfp2d::Fluid fracfluid = loadFluidProperties(j_fluid);

  if (j_params.count("Rock properties") != 1) {
    std::cout << "No rock properties input in  model parameters";
    il::abort();
  }
  json j_rock = j_params["Rock properties"];
  hfp2d::SolidProperties rock = loadSolidProperties(j_rock);

  if (j_params.count("Fracture height") != 1) {
    std::cout << "No Fracture height input in  model parameters";
    il::abort();
  }
  auto frac_heigth = j_params["Fracture height"].get<double>();

  // ----------
  // Create Initial fracture mesh
  // fracture are always assume to start at 90 from the well axis for now.
  //
  // -> need to have a vector of initial length of fracs
  //  and initial number of elements per frac.

  long nf = j_params.count("Initial fracture length");
  if (j_params.count("Initial fracture length") != 1) {
    std::cout << "No initial frac length in input file ";
    il::abort();
  }

  il::int_t nfracs = j_params["Initial fracture length"].size();
  // check consistency with number of perf.
  IL_EXPECT_FAST(nfracs == w_inj.coefPerf().size());
  // check consistency with number of perf.
  IL_EXPECT_FAST(nfracs == w_inj.hfLocation().size());

  // number of element per fracture.
  il::int_t nelts = 7;  // default value
  if (j_simul.count("Initial number of elements per fracture") == 1) {
    nelts = j_simul["Initial number of elements per fracture"];
    if (!(nelts %
          2)) {  // we ensure we have an odd number of elements per frac.
      nelts = nelts + 1;
    }
  }

  il::Array2D<double> xy{(nelts + 1) * nfracs, 2, 0.};
  il::Array2D<il::int_t> conn{nelts * nfracs, 2, 0};
  il::Array<il::int_t> source_loc_frac{nfracs, 0};

  double well_az =
      the_well.azimuth();  // is stored in degree not radians w.r. to y !
  // with respect to true North which is by definition the y-axis
  // double frac_az=the_well.azimuth()*(il::pi)/180.-il::pi/2.;

  double frac_az_x = the_well.azimuth() * (il::pi) / 180.;
  double cos_az = cos(frac_az_x);
  double sin_az = sin(frac_az_x);

  for (il::int_t f = 0; f < nfracs; f++) {
    double xo = the_well.coordinates(w_inj.hfLocation(f), 0);
    double yo = the_well.coordinates(w_inj.hfLocation(f), 1);
    double lo = j_params["Initial fracture length"][f];
    double h = 2 * lo / nelts;
    for (il::int_t i = (nelts + 1) * f; i < (nelts + 1) * (f + 1); i++) {
      xy(i, 0) = xo + ((i - (nelts + 1) * f) * h - lo) * cos_az;
      xy(i, 1) = yo + ((i - (nelts + 1) * f) * h - lo) * sin_az;
    }
    for (il::int_t i = nelts * f; i < (nelts) * (f + 1); i++) {
      conn(i, 0) = i + f;
      conn(i, 1) = i + 1 + f;
    }
    // floor 7/2 will give 3 which is the 4 elt in the mesh so it's ok
    source_loc_frac[f] = nelts / 2;
  }

  hfp2d::Mesh fracsMesh(xy, conn, 0);

  // zero initial rates entering the fracs
  il::Array<double> frac_rates{nfracs, 0.};
  hfp2d::Sources frac_inj(source_loc_frac, frac_rates);

  // in-situ conditions on current mesh .....
  //
  if (j_params.count("In-situ conditions") != 1) {
    std::cout << "No In-situ conditions input in  model parameters";
    il::abort();
  }
  json j_insitu = j_params["In-situ conditions"];
  hfp2d::InSituConditions inSituStress = loadInSitu(j_insitu);

  // Simulation Parameters....
  // default for now
  // for frac ehl flow
  hfp2d::SimulationParameters SimulFracParam;
  SimulFracParam.frac_front_max_its = 40;
  SimulFracParam.frac_front_tolerance = 1.e-3;
  SimulFracParam.ehl_relaxation = 0.95;
  SimulFracParam.ehl_tolerance = 1.e-6;

  // for pipe flow
  hfp2d::SimulationParameters WellFlowParam;
  WellFlowParam.ehl_relaxation = 1.0;
  WellFlowParam.ehl_tolerance = 1.e-6;

  // Create Initial state.....

  // initial well solution.
  hfp2d::WellSolution wellIniSol(the_well, w_inj, fracfluid);

  // in-situ traction on fracture mesh....
  il::Array<double> in_situ_traction =
      inSituStress.uniformAllInSituTractions(fracsMesh);

  il::Array<double> pc = wellIniSol.pressureAtClusters();

  il::Array<double> snc{nfracs, 0.};
  for (il::int_t i = 0; i < frac_inj.SourceElt().size(); i++) {
    snc[i] = in_situ_traction[frac_inj.SourceElt(i) * 2 + 1];
  }

  // as long as fluid pressure is lower than the in-situ stress we just loop on
  // the wb flow - pressurizing the well only.

  // one step.
  double dt = 0.1;

  int k = 0;
  hfp2d::WellSolution wellSol_n = wellIniSol;

  while (k < 10) {
    k++;

    hfp2d::WellSolution wellSol_n_p_1 =
        hfp2d::wellFlowSolverP0(wellSol_n, the_well, w_inj, fracfluid,
                                ffChurchill, dt, WellFlowParam, true);
    //
    wellSol_n = wellSol_n_p_1;

    pc = wellSol_n.pressureAtClusters();

    for (il::int_t i = 0; i < nfracs; i++) {
      std::cout << "------\n";
      std::cout << "pressure in well at cluster " << i << " is: " << pc[i]
                << " and normal in-situ stress " << snc[i] << "\n";
    }

    if (il::norm(pc, il::Norm::L2) > il::norm(snc, il::Norm::L2)) {
      std::cout << "fluid pressure just above in-situ stress ! \n";
      break;
    };
  }

  // compute elasticity matrix
  hfp2d::ElasticProperties elasprop = rock.ElasticProperties();

  il::Array2D<double> K = hfp2d::basic_assembly(
      fracsMesh, elasprop, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
      frac_heigth);

  // add tip correction for P0 for each tips in the mesh
  for (il::int_t i = 0; i < fracsMesh.tipElts().size(); i++) {
    AddTipCorrectionP0(fracsMesh, elasprop, fracsMesh.tipElts(i), K);
  }

  // initial net loading.
  il::Array<double> ftini{fracsMesh.numberDDDofs(), 0.},
      pfo{fracsMesh.numberDDDofs(), 0.};
  ftini = inSituStress.uniformAllInSituTractions(fracsMesh);
  for (il::int_t i = 0; i < nfracs; i++) {
    for (il::int_t j = 0; j < nelts; j++) {
      pfo[i * 2 * nelts + 2 * j + 1] = pc[i];
    }
  }

  il::blas(-1., ftini, il::io, pfo);

  std::cout << " initial elastic problem \n";
  il::Status status;
  // use a direct solver
  il::Array<double> dd_ini = il::linearSolve(K, pfo, il::io, status);
  status.ok();
  std::cout << "solved \n";
  for (il::int_t i = 0; i < fracsMesh.numberDDDofs(); i++) {
    std::cout << "traction " << pfo[i] << " dd " << dd_ini[i] << "\n";
  }

  il::Array<double> width{fracsMesh.numberOfElts(), 0.},
      sheardd{fracsMesh.numberOfElts(), 0.}, pf_o{fracsMesh.numberOfElts(), 0.};

  il::Array<double> sn_o{fracsMesh.numberOfElts(), 0.},
      tau_o{fracsMesh.numberOfElts(), 0.};

  for (il::int_t i = 0; i < fracsMesh.numberOfElts(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
    sn_o[i] = ftini[2 * i + 1];
    tau_o[i] = ftini[2 * i];
  }

  // create a solution at time t=0 object.
  hfp2d::Solution fracSol_n =
      hfp2d::Solution(fracsMesh, 0., width, sheardd, pf_o, sn_o, tau_o);
  il::Array<double> vel0{4, 0.};  // initial tip velocity
  fracSol_n.setTipsVelocity(vel0);

  il::Array<double> s0{4, 0.};  // initial ribbont tip distance
  for (il::int_t i = 0; i < fracsMesh.tipElts().size(); i++) {
    il::int_t e = fracsMesh.tipElts(i);
    s0[i] = 1.5 * fracsMesh.eltSize(e);
  }
  fracSol_n.setRibbonDistances(s0);

  // we know are ready
  il::Array<double> dpc{nfracs, 0.};

  hfp2d::MultiFracsSolution completeSol_n(fracSol_n, wellSol_n, frac_inj, dpc, 0.,
                                         0);


  return 0;
}

////////////////////////////////////////////////////////////////////////////////
hfp2d::MultiFracsSolution wellHFsSolver(hfp2d::MultiFracsSolution &Sol_n,
                                        double dt, hfp2d::WellMesh &w_mesh,
                                        hfp2d::WellInjection &w_inj,
                                        hfp2d::SimulationParameters &well_solver,
                                        hfp2d::SimulationParameters &frac_solver,
                                        il::io_t,
                                        il::Array2D<double> &K) {

  // We solve here the coupling between wellbore flow, fluid partitioning
  // between fractures
  // and hydraulic fracture propagation over one time-step.

  // Fixed point iteration scheme




  //  hfp2d::MultiFracsSolution newSol();
  return Sol_n;

}
}
