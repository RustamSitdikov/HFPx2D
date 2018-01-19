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
#include <src/input/json/LoadInputMultistage.h>
#include <src/solvers/MultiFracsSolver.h>
#include <src/wellbore/WellFlowP0.h>
#include <src/wellbore/WellSolution.h>

#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/Simplified3D.h>
#include <src/solvers/HFPropagationP0.h>
#include <src/solvers/MultiFracsSolution.h>

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

  // ...

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
  hfp2d::WellMesh well_mesh = loadWellMesh(j_wmesh);
  hfp2d::WellInjection w_inj = loadWellParameters(j_params, well_mesh);

  if (j_params.count("Fluid properties") != 1) {
    std::cout << "No fluid properties input in  model parameters";
    il::abort();
  }

  hfp2d::Sources well_sources = loadWellSource(j_params, well_mesh);

  // json j_fluid = j_params["Fluid properties"];
  hfp2d::Fluid fracfluid = loadFluidProperties(j_params["Fluid properties"]);

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
  //  IL_EXPECT_FAST(nfracs == w_inj.hfLocation().size());

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
      well_mesh.azimuth();  // is stored in degree not radians w.r. to y !
  // with respect to true North which is by definition the y-axis
  // double frac_az=well_mesh.azimuth()*(il::pi)/180.-il::pi/2.;

  double frac_az_x = well_mesh.azimuth() * (il::pi) / 180.;
  double cos_az = cos(frac_az_x);
  double sin_az = sin(frac_az_x);

  for (il::int_t f = 0; f < nfracs; f++) {
    double xo = well_mesh.coordinates(well_sources.SourceElt(f), 0);
    double yo = well_mesh.coordinates(well_sources.SourceElt(f), 1);
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
    source_loc_frac[f] = nelts / 2 + f * nelts;
  }

  hfp2d::Mesh fracsMesh(xy, conn, 0);

  // zero initial rates entering the fracs
  il::Array<double> frac_rates{nfracs, 0.};
  hfp2d::Sources frac_sources(source_loc_frac, frac_rates);

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
  hfp2d::SimulationParameters SimulWellFlowParam;
  SimulWellFlowParam.ehl_relaxation = 1.0;
  SimulWellFlowParam.ehl_tolerance = 1.e-6;

  // Create Initial state.....

  // initial well solution.
  hfp2d::WellSolution wellIniSol(well_mesh, w_inj, fracfluid);

  // in-situ traction on fracture mesh....
  il::Array<double> in_situ_traction =
      inSituStress.uniformAllInSituTractions(fracsMesh);

  il::Array<il::int_t> cluster_loc = well_sources.SourceElt();
  il::Array<double> pc = wellIniSol.pressureAtElts(cluster_loc);

  il::Array<double> snc{nfracs, 0.};
  for (il::int_t i = 0; i < frac_sources.SourceElt().size(); i++) {
    snc[i] = in_situ_traction[frac_sources.SourceElt(i) * 2 + 1];
  }

  // as long as fluid pressure is lower than the in-situ stress we just loop on
  // the wb flow - pressurizing the well only.

  // one step.
  double dt = 0.05;

  int k = 0;
  hfp2d::WellSolution wellSol_n = wellIniSol;

  while (k < 100) {
    k++;

    hfp2d::WellSolution wellSol_n_p_1 = hfp2d::wellFlowSolverP0(
        wellSol_n, well_mesh, w_inj, well_sources, ffChurchill, dt,
        SimulWellFlowParam, true, fracfluid);
    //
    wellSol_n = wellSol_n_p_1;

    pc = wellSol_n.pressureAtElts(cluster_loc);

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

  il::Array<double> pnet0=pfo;

  il::blas(-1., ftini, il::io, pnet0);

  std::cout << " initial elastic problem \n";
  il::Status status;
  // use a direct solver
  il::Array<double> dd_ini = il::linearSolve(K, pnet0, il::io, status);
  status.ok();
  std::cout << "solved \n";
  for (il::int_t i = 0; i < fracsMesh.numberDDDofs(); i++) {
    std::cout << "net load " << pnet0[i] << " dd " << dd_ini[i] << "\n";
  }

  il::Array<double> width{fracsMesh.numberOfElts(), 0.},
      sheardd{fracsMesh.numberOfElts(), 0.};

  il::Array<double> sn_o{fracsMesh.numberOfElts(), 0.},
      tau_o{fracsMesh.numberOfElts(), 0.},fluid_pressure_o{fracsMesh.numberOfElts(),0.};

  for (il::int_t i = 0; i < fracsMesh.numberOfElts(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
    sn_o[i] = ftini[2 * i + 1];
    tau_o[i] = ftini[2 * i];
    fluid_pressure_o[i]=pfo[2*i+1];
  }

  // create a solution at the current time of the wellsolution  object.
  hfp2d::Solution fracSol_n = hfp2d::Solution(
      fracsMesh, wellSol_n.time(), width, sheardd, fluid_pressure_o, sn_o, tau_o);
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

  hfp2d::MultiFracsSolution completeSol_n(fracSol_n, wellSol_n, frac_sources,
                                          well_sources, 0, dpc, 0.);

  // time step loop !!!

  il::int_t jt = 0;
  il::int_t nsteps = 1;

  dt = 0.001;

  while (jt < nsteps) {
    jt++;

    MultiFracsSolution completeSol_n_1 = wellHFsSolver(
        completeSol_n, dt, well_mesh, w_inj, fracfluid, rock, SimulFracParam,
        SimulWellFlowParam, frac_heigth, false, il::io, K);

    std::cout << "----------\n";
    std::cout << " Step " << jt << " time: " << completeSol_n_1.time();
    for (il::int_t i = 0; i < nfracs; i++) {
      std::cout << "influx in frac: " << i << " = "
                << completeSol_n_1.clusterFluxes(i) << "\n";
    }
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

// ENTRY FRICTION Residuals fction
double entryFrictionResiduals(double s,IFParamEntryFriction &params) {

  double res=params.dp-(params.fp)*(s*s)-(params.ft)*(pow(abs(s),params.beta_t));
  return res;
}




////////////////////////////////////////////////////////////////////////////////
// todo : add Sources in the Solution class ?
// todo : add wellbore mesh in well solution class ? for simplicity ?
// such that we have simpler api

////////////////////////////////////////////////////////////////////////////////
hfp2d::MultiFracsSolution wellHFsSolver(
    hfp2d::MultiFracsSolution &Sol_n, double dt, hfp2d::WellMesh &w_mesh,
    hfp2d::WellInjection &w_inj, hfp2d::Fluid &fracfluid,
    hfp2d::SolidProperties &rock, hfp2d::SimulationParameters &frac_solver_p,
    hfp2d::SimulationParameters &well_solver_p, double &frac_heigth, bool mute,
    il::io_t, il::Array2D<double> &K) {
  // We solve here the coupling between wellbore flow, fluid partitioning
  // between fractures
  // and hydraulic fracture propagation over one time-step.
  //
  // we start with rates at the previous time-step for the trial solution
  il::int_t nclusters = Sol_n.fracFluxes().size();

  hfp2d::Solution fracSol_n = Sol_n.fracSolution();
  hfp2d::WellSolution wellSol_n = Sol_n.wellSolution();

  hfp2d::Solution fracSol_k = Sol_n.fracSolution();
  hfp2d::WellSolution wellSol_k = Sol_n.wellSolution();

  hfp2d::Sources frac_sources_k = Sol_n.fracSources();
  hfp2d::Sources well_sources_k = Sol_n.wellSources();

  il::Array<il::int_t> frac_inj_loc = frac_sources_k.SourceElt();
  il::Array<il::int_t> well_out_loc = well_sources_k.SourceElt();

  il::Array<double> Q_in_k = Sol_n.fracFluxes();
  il::Array<double> Q_old = Q_in_k;
  il::Array<double> pc_w_k{nclusters, 0.}, pc_f_k{nclusters, 0.};
  il::Array<double> dpc_k{nclusters, 0.}, errQ{nclusters, 0.};

//  il::Array<double> allpf=fracSol_n.pressure();

  IFParamEntryFriction entry_struct;

  double pump_rate = w_inj.wellInjRate();
  double rr;


  std::cout << "+++++++++++++++++++++++++";

  // Fixed point iteration scheme
  //
  //
  double rela_flux = 0.01;  // pass as arguments ?
  double err = 1.;
  double Tolerance = 1.e-3;
  int k = 0;
  int kmax = 20;
  // note if all the fluxes are zero do not solve for frac flux, just the
  // wellbore
  while ((k < kmax) && (err > Tolerance)) {
    k++;

    frac_sources_k.setInjectionRates(Q_in_k);
    well_sources_k.setInjectionRates(Q_in_k);

    // solve for wellbore flow
    wellSol_k = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                        well_sources_k, ffChurchill, dt,
                                        well_solver_p, true, fracfluid);

    // solve for fracture propagation with given flux.
    //if (il::norm(Q_in_k, il::Norm::L2) != 0.) {
      fracSol_k = hfp2d::FractureFrontLoop(fracSol_n, fracfluid, rock,
                                           frac_sources_k, frac_heigth, dt,
                                           frac_solver_p, true, il::io_t(), K);
    //}
    // would need to have some checks for convergences of the solvers.....

    // extract pressure at clusters
    // wellbore side
    pc_w_k = wellSol_k.pressureAtElts(well_out_loc);

    // fracture side and compute new flow rates.
    for (il::int_t i = 0; i < nclusters; i++) {
      pc_f_k[i] = fracSol_k.pressure(frac_inj_loc[i]);
      dpc_k[i] = pc_w_k[i] - pc_f_k[i];
      // compute new flow rates.... the best would be to root find here ?

      entry_struct.dp=abs(dpc_k[i]);
      entry_struct.fp=w_inj.coefPerf(i);
      entry_struct.ft=w_inj.coefTort(i);
      entry_struct.beta_t=w_inj.betaTort(i);
//      std::cout << " Dp " << dpc_k[i] << "\n";
//      std::cout << entryFrictionResiduals(0.,entry_struct) <<"\n";
//      std::cout << entryFrictionResiduals(1.*pump_rate,entry_struct) <<"\n";

      rr = imf::brent((imf::ImplicitFunction)entryFrictionResiduals , entry_struct, 0.0, 1000.*pump_rate,
                      1e-6, 100, false);
//      std::cout << " flux " << rr <<  "Dp-F(Q) " <<  entryFrictionResiduals(rr,entry_struct) <<"\n";

      Q_in_k[i]=rr*(dpc_k[i]/abs(dpc_k[i]));

//      if (Q_in_k[i] == 0.) {  // to ensure we don t have infinite flux.
//        if (w_inj.coefTort(i) != 0.) {
//          Q_in_k[i] = std::pow(std::abs(dpc_k[i]), 1. / w_inj.betaTort(i)) /
//                      w_inj.coefTort(i);
//        } else {
//          Q_in_k[i] = std::sqrt(std::abs(dpc_k[i])) /
//                      w_inj.coefPerf(i);  // ok provided coefPerf is non-zero.
//        }
//      }
//
//      Q_in_k[i] = dpc_k[i] / (w_inj.coefPerf(i) * abs(Q_in_k[i]) +
//                              w_inj.coefTort(i) *
//                                  std::pow(abs(Q_in_k[i]), w_inj.betaTort(i) - 1.0));
    }

    if (!mute) {
      std::cout << "-------\n";
      std::cout << "Its " << k << " DP : ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << " : " << dpc_k[i] << " ";
      }
      std::cout << "error: " << err << "\n";
    }

    // underelaxation of the flux....
    il::blas((1. - rela_flux), Q_old, rela_flux, il::io, Q_in_k);
    // compute successive relative difference, L2 norm

    for (il::int_t i = 0; i < nclusters; i++) {
      errQ[i] = abs((Q_in_k[i] - Q_old[i]) / Q_in_k[i]);
    }
    Q_old = Q_in_k;

    err = il::norm(errQ, il::Norm::L2);
    if (!mute) {
      std::cout << "Its " << k << " new fluxes: ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << " : " << Q_in_k[i] << " ";
      }
      std::cout << "error: " << err << "\n";
    }
  }

  if (!mute) {
    if (err < Tolerance) {
      std::cout << "Well - HFs  converged after " << k << " Its\n";
    } else
      std::cout << "Well - HFs NOT converged after " << k
                << " Its, error is: " << err << "\n";
  }

  hfp2d::MultiFracsSolution newSol(fracSol_k, wellSol_k, frac_sources_k,
                                   well_sources_k, k, dpc_k, err);

  return newSol;
}
}
