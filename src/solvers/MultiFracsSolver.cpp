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
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra/dense/factorization/linearSolve.h>
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

  hfp2d::Sources well_sources = loadWellSource(j_params, well_mesh);

  if (j_params.count("Fluid properties") != 1) {
    std::cout << "No fluid properties input in  model parameters";
    il::abort();
  }

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
    if (!(nelts % 2)) {
      // we ensure we have an odd number of elements per frac.
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
  il::Array<double> pnet0 = pfo;
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
      tau_o{fracsMesh.numberOfElts(), 0.},
      fluid_pressure_o{fracsMesh.numberOfElts(), 0.};

  for (il::int_t i = 0; i < fracsMesh.numberOfElts(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
    sn_o[i] = ftini[2 * i + 1];
    tau_o[i] = ftini[2 * i];
    fluid_pressure_o[i] = pfo[2 * i + 1];
  }

  // create a solution at the current time of the wellsolution  object.
  hfp2d::Solution fracSol_n =
      hfp2d::Solution(fracsMesh, wellSol_n.time(), width, sheardd,
                      fluid_pressure_o, sn_o, tau_o);
  il::Array<double> vel0{4, 0.};  // initial tip velocity
  fracSol_n.setTipsVelocity(vel0);

  il::Array<double> s0{4, 0.};  // initial ribbon tip distance
  for (il::int_t i = 0; i < fracsMesh.tipElts().size(); i++) {
    il::int_t e = fracsMesh.tipElts(i);
    s0[i] = 1.5 * fracsMesh.eltSize(e);
  }
  fracSol_n.setRibbonDistances(s0);

  // we know are ready
  il::Array<double> dpc{nfracs, 0.};

  hfp2d::MultiFracsSolution completeSol_n(fracSol_n, wellSol_n, frac_sources,
                                          well_sources, 0, dpc, 0.);

  double max_time = 1.;
  if ( j_simul.count("Maximum time") ==1 ) {
    max_time = j_simul["Maximum time"].get<double>();
  }
  il::int_t max_steps = 150;
  if ( j_simul.count("Maximum number of steps") ==1 ) {
    max_steps = j_simul["Maximum number of steps"].get<long>();
  }

  dt = 0.002;
  if ( j_simul.count("Time step") ==1 ) {
    dt = j_simul["Time step"].get<double>();
  }

  double dt_min = 0.00001;
  if ( j_simul.count("Minimum time step") ==1 ) {
    dt_min = j_simul["Minimum time step"].get<double>();
  }

  il::int_t jt = 0;

  // time step loop !!!

  while ( (jt < max_steps)  && (completeSol_n.time() < max_time) ) {
    jt++;

    MultiFracsSolution completeSol_n_1 = wellHFsSolver_fixedpts(
        completeSol_n, dt, well_mesh, w_inj, fracfluid, rock, SimulFracParam,
        SimulWellFlowParam, frac_heigth, false, il::io, K);

    std::cout << "----------" << std::endl;
    std::cout << "Step " << jt << "; time: " << completeSol_n_1.time()
              << std::endl;
    for (il::int_t i = 0; i < nfracs; i++) {
      std::cout << "influx in frac. " << i << ": "
                << completeSol_n_1.clusterFluxes(i) << "; DP in frac. " << i
                << ": " << completeSol_n_1.dpEntries(i) << std::endl;
    }
    std::cout << "----------" << std::endl;
    completeSol_n = completeSol_n_1;
    // todo time step acceptance test
    // todo adaptative time-step
    // saving of solution.

  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////////

// ENTRY FRICTION Residuals
double entryFrictionResiduals(double s, IFParamEntryFriction &params) {
  double sgn_s = ((s < 0) ? -1.0 : double((s > 0)));
  double res = params.dp -
               sgn_s * ((params.fp) * (s * s) +
                        (params.ft) * (std::pow(std::fabs(s), params.beta_t)));
  return res;
}

////////////////////////////////////////////////////////////////////////////////
// todo : add Sources in the Solution class ?
// todo : add wellbore mesh in well solution class ? for simplicity ?
// such that we have simpler api

////////////////////////////////////////////////////////////////////////////////
hfp2d::MultiFracsSolution wellHFsSolver_fixedpts(
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

  hfp2d::Solution fracSol_var = Sol_n.fracSolution();
  hfp2d::WellSolution wellSol_var = Sol_n.wellSolution();

  hfp2d::Sources frac_sources_k = Sol_n.fracSources();
  hfp2d::Sources well_sources_k = Sol_n.wellSources();

  hfp2d::Sources frac_sources_var = Sol_n.fracSources();
  hfp2d::Sources well_sources_var = Sol_n.wellSources();

  il::Array<il::int_t> frac_inj_loc = frac_sources_k.SourceElt();
  il::Array<il::int_t> well_out_loc = well_sources_k.SourceElt();

  il::Array<double> Q_in_k = Sol_n.fracFluxes();
  il::Array<double> Q_old = Q_in_k;
  il::Array<double> pc_w_k{nclusters, 0.}, pc_f_k{nclusters, 0.};
  il::Array<double> dpc_k{nclusters, 0.};
  il::Array<double> errQ{nclusters, 0.}, errDP{nclusters, 0.};
  double err = 1.;
  il::Array<double> res_v{nclusters, 0.};
  il::Array<double> dQ_v{nclusters, 0.};

  //  il::Array<double> allpf=fracSol_n.pressure();

  IFParamEntryFriction entry_struct;

  double pump_rate = w_inj.wellInjRate();
  //  double rr;

  il::Array2D<double> Jacob{nclusters, nclusters, 0.};
  il::Status status;
  //  il::LU<il::Array2D<double>> Jacob_LU();
  il::Array2D<double> Jacob_inv;

  // remember that the rate entering the fracture is rate / fracture heigth
  il::Array<double> rates_per_height{nclusters,0.};

  if (!mute) {
    std::cout << "+++++++++++++++++++++++++" << std::endl;
  }

  // todo: start with non-zero fluxes - solve the ENTRY FRICTION Residuals
  // eqn(?)

  // Quasi-Newton iteration scheme
  //
  //
  // todo: pass as arguments OR move to numerical parameters file
  double dQn = 2.0e-8;
  int num_J_reuse = 1;
  double rela_flux = 1.;  // 0.01;
  double Tolerance = 1.e-6;
  int kmax = 20;


  int k = 0;
  // note: if all the fluxes are zero do not solve for frac flux,
  // just the wellbore
  while ((k < kmax) && (err > Tolerance)) {
    k++;

    for (il::int_t i=0;i<nclusters;i++){
      rates_per_height[i]=Q_in_k[i]/frac_heigth;
    }
    frac_sources_k.setInjectionRates(rates_per_height);

    well_sources_k.setInjectionRates(Q_in_k);

    if (!mute) {
      std::cout << "-------" << std::endl;
    }

    // solve for wellbore flow
    wellSol_k = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                        well_sources_k, ffChurchill, dt,
                                        well_solver_p, true, fracfluid);

    // solve for fracture propagation with given flux.
    // if (il::norm(Q_in_k, il::Norm::L2) != 0.) {
    fracSol_k = hfp2d::FractureFrontLoop(fracSol_n, fracfluid, rock,
                                         frac_sources_k, frac_heigth, dt,
                                         frac_solver_p, true, il::io_t(), K);
    //}
    // would need to have some checks for convergences of the solvers.....

    // extract pressure at clusters
    // wellbore side
    pc_w_k = wellSol_k.pressureAtElts(well_out_loc);

    // fracture side
    for (il::int_t i = 0; i < nclusters; i++) {
      pc_f_k[i] = fracSol_k.pressure(frac_inj_loc[i]);
      dpc_k[i] = pc_w_k[i] - pc_f_k[i];
    }

    // loop over cluster w. small variations of Q_in to estimate Jacobian
    if ((k - 1) % num_J_reuse == 0) {
      il::Array<double> Q_in_var = Q_in_k;
      il::Array<double> pc_w_var{nclusters, 0.}, pc_f_var{nclusters, 0.};
      for (il::int_t i = 0; i < nclusters; i++) {
        double abs_Q = std::fabs(Q_in_k[i]);
        // todo: scale it properly
        double scl_Q =
            ((Q_in_k[i] == 0.) ? std::fabs(pump_rate) / nclusters : abs_Q);
        // sign(Q_in_k[i])
        // double sgn_Q = ((Q_in_k[i] < 0) ? -1.0 : double((Q_in_k[i] > 0)));
        double sgn_Q = ((Q_in_k[i] < 0) ? -1.0 : 1.0);
        double dQi = std::fabs(dQn * scl_Q) * sgn_Q;

        Q_in_var[i] = Q_in_k[i] + dQi;

        frac_sources_var.setInjectionRates(Q_in_var);
        well_sources_var.setInjectionRates(Q_in_var);

        // solve for wellbore flow
        wellSol_var = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                              well_sources_var, ffChurchill, dt,
                                              well_solver_p, true, fracfluid);

        pc_w_var = wellSol_var.pressureAtElts(well_out_loc);

        // solve for fracture propagation with given flux.
        fracSol_var = hfp2d::FractureFrontLoop(
            fracSol_n, fracfluid, rock, frac_sources_var, frac_heigth, dt,
            frac_solver_p, true, il::io_t(), K);

        // estimate Jacobian
        for (il::int_t j = 0; j < nclusters; j++) {
          pc_f_var[j] = fracSol_var.pressure(frac_inj_loc[j]);
          Jacob(j, i) = ((pc_w_var[j] - pc_f_var[j]) - (dpc_k[j])) / (dQi);
        }
        // Well - HF coupling (ENTRY FRICTION derivatives)
        // todo: check the sign!
        if (Q_in_k[i] != 0.) {
          Jacob(i, i) -= 2. * w_inj.coefPerf(i) * abs_Q +
                         w_inj.coefTort(i) * sgn_Q *
                             std::pow(abs_Q, w_inj.betaTort(i) - 1.);
        }

        Q_in_var[i] = Q_in_k[i];
      }

      // LU decomposition of Jacobian
      il::LU<il::Array2D<double>> Jacob_LU(Jacob, il::io, status);
      status.abortOnError();

      // echo...
      if (!mute) {
        std::cout << "-------" << std::endl;
        std::cout << "Det(J): " << Jacob_LU.determinant() << std::endl;
        //          std::cout << "Cond(J): "
        //                    << Jacob_LU.conditionNumber(il::Norm::L1, ) <<
        //                    std::endl;
      }

      // inverse Jacobian
      Jacob_inv = Jacob_LU.inverse();
    }

    // residuals...
    for (il::int_t i = 0; i < nclusters; i++) {
      entry_struct.dp = dpc_k[i];
      entry_struct.fp = w_inj.coefPerf(i);
      entry_struct.ft = w_inj.coefTort(i);
      entry_struct.beta_t = w_inj.betaTort(i);
      // todo: check the sign!
      res_v[i] = -entryFrictionResiduals(Q_in_k[i], entry_struct);
    }

    // solution...
    // (use LU decomposition or inverse of Jacobian)
    //    dQ_v = Jacob_LU.solve(res_v);
    dQ_v = il::dot(Jacob_inv, res_v);
    //    dQ_v = il::linearSolve(Jacob, res_v, il::io, status);
    //    status.abortOnError();

    // compute new flow rates...
    for (il::int_t i = 0; i < nclusters; i++) {
      Q_in_k[i] += dQ_v[i] * rela_flux;
    }

    // echo...
    if (!mute) {
      std::cout << "-------" << std::endl;
      std::cout << "Its " << k << "; DP: ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << ": " << dpc_k[i] << "; ";
      }
      std::cout << "prev. error: " << err << std::endl;
    }

    //    // under-relaxation of the flow rates...
    //    il::blas((1. - rela_flux), Q_old, rela_flux, il::io, Q_in_k);

    // compute successive relative difference, L2 norm
    for (il::int_t i = 0; i < nclusters; i++) {
      errQ[i] = abs((Q_in_k[i] - Q_old[i]) / Q_in_k[i]);
      errDP[i] = abs(res_v[i] / dpc_k[i]);
    }
    err = il::norm(errQ, il::Norm::L2);  // + il::norm(errDP, il::Norm::L2);

    // echo...
    if (!mute) {
      std::cout << "old fluxes: ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << ": " << Q_old[i] << "; ";
      }
      std::cout << std::endl << "res. of DP: ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << ": " << res_v[i] << "; ";
      }
      std::cout << std::endl << "new fluxes: ";
      for (il::int_t i = 0; i < nclusters; i++) {
        std::cout << i << ": " << Q_in_k[i] << "; ";
      }
      std::cout << "error: " << err << std::endl;
    }

    // resume...
    Q_old = Q_in_k;
  }

  if (!mute) {
    if (err < Tolerance) {
      std::cout << "Well - HFs coupling converged after " << k << " its"
                << std::endl;
    } else
      std::cout << "Well - HFs coupling NOT converged after " << k
                << " its; error is: " << err << std::endl;
  }

  hfp2d::MultiFracsSolution newSol(fracSol_k, wellSol_k, frac_sources_k,
                                   well_sources_k, k, dpc_k, err);

  return newSol;
}
}
