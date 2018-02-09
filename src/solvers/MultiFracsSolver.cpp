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

int MultipleFracsPropagation(std::string &filename) {
  // routine for the propagation of multiple Heigth contained HFs from a
  // horizontal wellbore
  // DEBUGGING

  // ...

  std::string wellfilename = filename.empty()
                             ? "../Debug/WellTest.json"
                             : filename;

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
  auto frac_height = j_params["Fracture height"].get<double>();

  // ----------
  // Create Initial fracture mesh
  // fracture are always assume to start at 90 from the well axis for now.
  //
  // -> need to have a vector of initial length of fracs
  //  and initial number of elements per frac.

  if (j_params.count("Clusters MD") != 1) {
    std::cout << "No clusters MD input in  model parameters";
    il::abort();
  }

//  long nf = j_params.count("Initial fracture length");
  if (j_params.count("Initial fracture length") != 1) {
    std::cout << "No initial frac length in input file ";
    il::abort();
  }

  il::int_t nfracs = j_params["Initial fracture length"].size();
  // check consistency with number of perf.
  IL_EXPECT_FAST(nfracs == w_inj.coefPerf().size());
  IL_EXPECT_FAST(nfracs == j_params["Clusters MD"].size());
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

  il::Array<il::int_t> fracID{nfracs*nelts,0};

  double well_az = well_mesh.azimuth();
  // is stored in degree not radians w.r. to y !
  // with respect to true North which is by convention the y-axis
  // double frac_az=well_mesh.azimuth()*(il::pi)/180.-il::pi/2.;

  double frac_az_x = well_az * (il::pi) / 180.;
  double cos_az = cos(frac_az_x);
  double sin_az = sin(frac_az_x);

  // todo: "Clusters MD" should be used to avoid fracs overlapping!!!!!
  for (il::int_t f = 0; f < nfracs; f++) {
    // well mesh element containing f-th cluster
    il::int_t well_cluster_elt = well_sources.SourceElt(f);
    // closest upstream well mesh node
    il::int_t ref_well_node = well_mesh.connectivity(well_cluster_elt, 0);
    // distance from this node to the f-th cluster
    double so = (double)j_params["Clusters MD"][f]
                - well_mesh.md(ref_well_node);
    // location of the f-th cluster in Cartesian coordinates
    double xo = well_mesh.coordinates(ref_well_node, 0) + so * sin_az;
    double yo = well_mesh.coordinates(ref_well_node, 1) + so * cos_az;
    // frac initial length and element size
    double lo = j_params["Initial fracture length"][f];
    double h = 2 * lo / nelts;
    // initializing the f-th frac mesh
    for (il::int_t i = (nelts + 1) * f; i < (nelts + 1) * (f + 1); i++) {
      xy(i, 0) = xo + ((i - (nelts + 1) * f) * h - lo) * cos_az;
      xy(i, 1) = yo + ((i - (nelts + 1) * f) * h - lo) * sin_az;
    }
    for (il::int_t i = nelts * f; i < (nelts) * (f + 1); i++) {
      fracID[i]=f;
      conn(i, 0) = i + f;
      conn(i, 1) = i + 1 + f;
    }
    // floor 7/2 will give 3 which is the 4 elt in the mesh so it's ok
    source_loc_frac[f] = nelts / 2 + f * nelts;
  }

  il::Array<il::int_t> MatID{nfracs*nelts,0};

  hfp2d::Mesh fracsMesh(xy, conn,MatID,fracID, 0);

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
  SimulFracParam.frac_front_max_its = 50;
  SimulFracParam.frac_front_tolerance = 1.e-3;
  SimulFracParam.ehl_relaxation = 0.95;
  SimulFracParam.ehl_tolerance = 1.e-6;

  // for pipe (well) flow
  hfp2d::SimulationParameters SimulWellFlowParam;
  SimulWellFlowParam.ehl_relaxation = 1.0;
  SimulWellFlowParam.ehl_tolerance = 1.e-6;

  // for well-frac coupling
  hfp2d::SimulationParameters SimulCouplParam;
  SimulCouplParam.ehl_max_its = 20;
  SimulCouplParam.ehl_relaxation = 1.0;
  SimulCouplParam.ehl_tolerance = 1.e-6;

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
      frac_height);

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

  std::cout << "initial elastic problem ";
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
  il::Array<double> vel0{2 * nfracs, 0.};  // initial tip velocity
  fracSol_n.setTipsVelocity(vel0);

  il::Array<double> s0{2 * nfracs, 0.};  // initial ribbon tip distance
  for (il::int_t i = 0; i < fracsMesh.tipElts().size(); i++) {
    il::int_t e = fracsMesh.tipElts(i);
    s0[i] = 1.5 * fracsMesh.eltSize(e);
  }
  fracSol_n.setRibbonDistances(s0);

  // we know are ready
  il::Array<double> dpc{nfracs, 0.};

  hfp2d::MultiFracsSolution completeSol_n(fracSol_n, wellSol_n,
                                          frac_sources, well_sources,
                                          frac_height,
                                          0, dpc, 0., 0.);

  double max_time = 1.;
  if ( j_simul.count("Maximum time") == 1 ) {
    max_time = j_simul["Maximum time"].get<double>();
  }
  il::int_t max_steps = 150;
  if ( j_simul.count("Maximum number of steps") == 1 ) {
    max_steps = j_simul["Maximum number of steps"].get<long>();
  }

  dt = 0.01;
  if ( j_simul.count("Time step") == 1 ) {
    dt = j_simul["Time step"].get<double>();
  }

  double dt_min = 0.00001;
  if ( j_simul.count("Minimum time step") == 1 ) {
    dt_min = j_simul["Minimum time step"].get<double>();
  }

  il::int_t jt = 0; // ea = 0;

  std::string basefilename = "HFs_well_coupling_"; // default name
  if (js.count("Results files name core") == 1){
    basefilename = js["Results files name core"].get<std::string>();
  }

  std::string resfilename;

  double max_tip_v, min_ds; il::int_t num_f_its;

  // step acceptance (on FF loop convergence)
  bool accept = true;

  // time step loop !!!
  while ( jt < max_steps && completeSol_n.time() < max_time ) {
    jt++;

    MultiFracsSolution completeSol_n_1 = wellHFsSolver_fixedpts(
        completeSol_n, dt, well_mesh, w_inj, fracfluid, rock,
        SimulFracParam, SimulWellFlowParam, SimulCouplParam,
        frac_height,
        true, il::io, K, accept);

    if (accept) {
      std::cout << "-------------------------" << std::endl;
      std::cout << "Step # " << jt << "; time: " << completeSol_n_1.time()
                << std::endl;
      for (il::int_t i = 0; i < nfracs; i++) {
        std::cout << "influx in frac. " << i+1 << ": "
                  << completeSol_n_1.clusterFluxes(i)
                  << "; press.drop in frac. " << i+1 << ": "
                  << completeSol_n_1.dpEntries(i)
                  << std::endl;
      }
      std::cout << "coupling error: " << completeSol_n_1.errFluxes()
                << "; coupling its.: " << completeSol_n_1.itsCoupling()
                << "; DP error: " << completeSol_n_1.errDps()
                << std::endl;
      std::cout << "frac.tip error: " << completeSol_n_1.fracSolution().errFront()
                << "; tip loop its.: " << completeSol_n_1.fracSolution().frontIts()
                << "; EHL error (w): " << completeSol_n_1.fracSolution().errOpening()
                << "; new number of elts.: " << completeSol_n_1.fracSolution().currentMesh().numberOfElts()
                << std::endl;

      completeSol_n = completeSol_n_1;

      for (il::int_t i = 0;
           i < completeSol_n.fracSolution().currentMesh().tipElts().size();
           i++) {
        std::cout << "tip " << i+1 << " velocity: "
                  << completeSol_n.fracSolution().tipsVelocity()[i] << "; ";
      }
      std::cout << std::endl;

      // saving solution.
      resfilename = basefilename + std::to_string(jt) + ".json";
      completeSol_n.writeToFile(resfilename);

      // todo: re-start

      // adaptive time-step
      max_tip_v = il::max(completeSol_n.fracSolution().tipsVelocity());
      min_ds = il::min(completeSol_n.fracSolution().currentMesh().allEltSize());
      num_f_its = completeSol_n.fracSolution().frontIts();

      if ((max_tip_v > 0.0)) {  //&& (fracSol_n.tipsLocation()(1,0)>1.5)
        double dt_new = 0.5 * min_ds / max_tip_v;
        if (num_f_its > 10) {dt_new *= 5.0 / (double)num_f_its;}
        // modify to more clever ?
        if (dt_new > 1.25 * dt) {
          dt = 1.25 * dt;
        } else {
          if (dt_new < dt * 0.8) {
            dt = std::max(0.8 * dt, dt_min);
          } else {
            dt = std::max(dt_new, dt_min);
          }
        }
      }

    } else {
      // reject time step
      std::cout << "-------------------------" << std::endl;
      std::cout << "step # " << jt << " rejected; time: "
                << completeSol_n_1.time()
                << "; time step: " << dt << std::endl;
      std::cout << "Error on frac front " << completeSol_n_1.fracSolution().errFront()
                << " after " << completeSol_n_1.fracSolution().frontIts() << " its" << std::endl;

      // reduce time step
      if (dt / 2. >= dt_min) {
        dt = dt / 2.;
        std::cout << "Reducing time step to "
                  << dt << "s. in order to re-try \n";
        jt--;
      } else {
        std::cout << "Error on frac. front too large with small time steps. - "
                "stop simulation" << std::endl;
        break;
      }
    }
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
// todo: add Sources in the Solution class ?
// todo: add wellbore mesh in well solution class ? for simplicity ?
// such that we have simpler api

////////////////////////////////////////////////////////////////////////////////
hfp2d::MultiFracsSolution wellHFsSolver_fixedpts(
    hfp2d::MultiFracsSolution &Sol_n, double dt,
    hfp2d::WellMesh &w_mesh,
    hfp2d::WellInjection &w_inj,
    hfp2d::Fluid &fracfluid,
    hfp2d::SolidProperties &rock,
    hfp2d::SimulationParameters &frac_solver_p,
    hfp2d::SimulationParameters &well_solver_p,
    hfp2d::SimulationParameters &coupling_p,
    double frac_height,
    bool mute,
    il::io_t, il::Array2D<double> &K, bool &accept) {
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

  hfp2d::Solution fracSol_var; // = fracSol_k;
  hfp2d::WellSolution wellSol_var; // = wellSol_k;
  il::Array2D<double> K_pre = K;

  hfp2d::Sources frac_sources_k = Sol_n.fracSources();
  hfp2d::Sources well_sources_k = Sol_n.wellSources();

  hfp2d::Sources frac_sources_var = Sol_n.fracSources();
  hfp2d::Sources well_sources_var = Sol_n.wellSources();

  il::Array<il::int_t> frac_inj_loc = frac_sources_k.SourceElt();
  il::Array<il::int_t> well_out_loc = well_sources_k.SourceElt();

  // current well-to-fracs fluxes (in m^3/s)
  il::Array<double> rates_cur = Sol_n.clusterFluxes();
  // this is used in iteration loop to estimete errors
  il::Array<double> rates_old = rates_cur;
  // rate entering the fracture = rate escaping the well / fracture height
  il::Array<double> rates_per_height{nclusters, 0.};
  // pressures on well and frac side
  il::Array<double> pc_w_k{nclusters, 0.}, pc_f_k{nclusters, 0.};
  // pressure drops between well and fracs
  il::Array<double> dpc_k{nclusters, 0.};
  // errors of fluxes and pressure drops
  il::Array<double> errQ{nclusters, 0.}, errDP{nclusters, 0.};
  double err_f = 1., err_p = 1.;
  accept = true;
  // residuals
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

//  if (!mute) {
//    std::cout << "+++++++++++++++++++++++++" << std::endl;
//  }

  // Quasi-Newton iteration scheme
  //
  //
  // todo: relate it to the machine precision
  double dQn = 1.49e-8; // or >= sqrt of machine precision
  // todo: pass as argument OR move to numerical parameters file
  int num_J_reuse = 2; // save time for Jacobian re-calculation

  double rela_flux = coupling_p.ehl_relaxation;
  double coupl_tolerance = coupling_p.ehl_tolerance;
  il::int_t max_iters = coupling_p.ehl_max_its;


  il::int_t k = 0;
  // note: if all the fluxes are zero do not solve for frac flux,
  // just the wellbore
  while ((k < max_iters) && (err_f > coupl_tolerance)) {
    k++;

    K = K_pre;

    for (il::int_t i = 0; i < nclusters; i++) {
      rates_per_height[i] = rates_cur[i] / frac_height;
    }
    // rate entering the fracture = rate escaping the well / fracture height
    frac_sources_k.setInjectionRates(rates_per_height);
    well_sources_k.setInjectionRates(rates_cur);

    // solve for wellbore flow
    wellSol_k = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                        well_sources_k, ffChurchill, dt,
                                        well_solver_p, true, fracfluid);

    // solve for fracture propagation with given flux
    // if (il::norm(rates_cur, il::Norm::L2) != 0.) {
    fracSol_k = hfp2d::FractureFrontLoop(fracSol_n, fracfluid, rock,
                                         frac_sources_k, frac_height, dt,
                                         frac_solver_p, true, il::io, K);
    //}
    // would need to have some checks for convergences of the solvers....
    if (fracSol_k.errFront() >= frac_solver_p.frac_front_tolerance) {
      accept = false;
    }

    // echo...
    if (!mute) {
      std::cout << "-------------------------" << std::endl;
      std::cout << "coupling iter-n " << k << std::endl;
      std::cout << "-------" << std::endl;
      if (!accept) {
        std::cout << "fracture front loop not converged" << std::endl;
      }
      std::cout << "base FF loop its: " << fracSol_k.frontIts()
                << "; front error: " << fracSol_k.errFront()
                << "; EHL error (w): "<< fracSol_k.errOpening()
                << "; new number of elts: " << fracSol_k.currentMesh().numberOfElts()
                << std::endl;
    }

    if (!accept) {
      break;
    }

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
      il::Array<double> rates_var = rates_cur;
      il::Array<double> pc_w_var{nclusters, 0.}, pc_f_var{nclusters, 0.};

      for (il::int_t i = 0; i < nclusters; i++) {
        // reset the elasticity matrix to the previous time step
        K = K_pre;
        //
        double abs_Q = std::fabs(rates_cur[i]);
        // scale it properly
        double scl_Q =
            ((rates_cur[i] == 0.) ? std::fabs(pump_rate) / nclusters : abs_Q);
        // sign(rates_cur[i])
        // double sgn_Q = ((rates_cur[i] < 0) ? -1.0 : double((rates_cur[i] > 0)));
        double sgn_dQ = ((rates_cur[i] < 0) ? -1.0 : 1.0);
        double dQi = std::fabs(dQn * scl_Q) * sgn_dQ;

        // introduce variation of i-th flux
        rates_var[i] = rates_cur[i] + dQi;

        for (il::int_t j = 0; j < nclusters; j++) {
          rates_per_height[j] = rates_var[j] / frac_height;
        }
        // rate entering the fracture = rate escaping the well / fracture height
        frac_sources_var.setInjectionRates(rates_per_height);
        well_sources_var.setInjectionRates(rates_var);

        // solve for wellbore flow with perturbed i-th flux
        wellSol_var = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                              well_sources_var, ffChurchill, dt,
                                              well_solver_p, true, fracfluid);

        pc_w_var = wellSol_var.pressureAtElts(well_out_loc);

        // solve for fracture propagation with perturbed i-th flux
        fracSol_var = hfp2d::FractureFrontLoop(
            fracSol_n, fracfluid, rock, frac_sources_var, frac_height, dt,
            frac_solver_p, true, il::io, K);

        if (fracSol_var.errFront() >= frac_solver_p.frac_front_tolerance) {
          accept = false;
        }

        // echo...
        if (!mute) {
          if (!accept) {
            std::cout << "trial fracture front loop not converged" << std::endl;
          }
          std::cout << "trial FF loop its: " << fracSol_var.frontIts()
                    << "; front error: " << fracSol_var.errFront()
                    << "; EHL error (w): "<< fracSol_var.errOpening()
                    << "; new number of elts: " << fracSol_var.currentMesh().numberOfElts()
                    << std::endl;
        }

        if (accept) {
          // estimate Jacobian
          for (il::int_t j = 0; j < nclusters; j++) {
            pc_f_var[j] = fracSol_var.pressure(frac_inj_loc[j]);
            Jacob(j, i) = ((pc_w_var[j] - pc_f_var[j]) - (dpc_k[j])) / (dQi);
          }
          // Well - HF coupling (ENTRY FRICTION) derivatives
          // check the sign!
          if (rates_cur[i] != 0.) {
            Jacob(i, i) -= 2. * w_inj.coefPerf(i) * abs_Q
                           + w_inj.coefTort(i) * sgn_dQ *
                             std::pow(abs_Q, w_inj.betaTort(i) - 1.);
          }
        }

        // reset...
        rates_var[i] = rates_cur[i];
      }

      if (accept) {
        // LU decomposition of Jacobian
        il::LU<il::Array2D<double>> Jacob_LU(Jacob, il::io, status);
        status.abortOnError();

        // echo...
        if (!mute) {
//        if (!accept) {
//          std::cout << "fracture front loop not converged" << std::endl;
//        }
//        std::cout << "fr. front loop its: " << fracSol_k.frontIts()
//                  << "; front error: " << fracSol_k.errFront()
//                  << "; new number of elts: " << fracSol_k.currentMesh().numberOfElts()
//                  << "; error w EHL "<< fracSol_k.errOpening() << std::endl;
          std::cout << "-------" << std::endl;
          std::cout << "Det(J): " << Jacob_LU.determinant() << std::endl;
          //          std::cout << "Cond(J): "
          //                    << Jacob_LU.conditionNumber(il::Norm::L1, ) <<
          //                    std::endl;
        }

        // inverse Jacobian
        Jacob_inv = Jacob_LU.inverse();
      } else { break; }
    }

    if (accept) {
      // residuals...
      for (il::int_t i = 0; i < nclusters; i++) {
        entry_struct.dp = dpc_k[i];
        entry_struct.fp = w_inj.coefPerf(i);
        entry_struct.ft = w_inj.coefTort(i);
        entry_struct.beta_t = w_inj.betaTort(i);
        // check the sign!
        res_v[i] = -entryFrictionResiduals(rates_cur[i], entry_struct);
      }

      // solution...
      // (use LU decomposition or inverse of Jacobian)
      //    dQ_v = Jacob_LU.solve(res_v);
      dQ_v = il::dot(Jacob_inv, res_v);
      //    dQ_v = il::linearSolve(Jacob, res_v, il::io, status);
      //    status.abortOnError();

      // compute new flow rates...
      for (il::int_t i = 0; i < nclusters; i++) {
        rates_cur[i] += dQ_v[i] * rela_flux; // with kinda under-relaxation
      }

      // echo...
      if (!mute) {
        std::cout << "-------" << std::endl;
        // std::cout << "coupling iter-n " << k << "; "
        std::cout << "DP: ";
        for (il::int_t i = 0; i < nclusters; i++) {
          std::cout << i+1 << ": " << dpc_k[i] << "; ";
        }
        std::cout << "prev. error: " << err_f << std::endl;
      }

//    // under-relaxation of the flow rates...
//    il::blas((1. - rela_flux), rates_old, rela_flux, il::io, rates_cur);

      // compute successive relative difference, L2 norm
      for (il::int_t i = 0; i < nclusters; i++) {
        errQ[i] = abs((rates_cur[i] - rates_old[i]) / rates_cur[i]);
        errDP[i] = abs(res_v[i] / dpc_k[i]);
      }
      err_f = il::norm(errQ, il::Norm::L2);
      err_p = il::norm(errDP, il::Norm::L2);

      // echo...
      if (!mute) {
        std::cout << "old fluxes: ";
        for (il::int_t i = 0; i < nclusters; i++) {
          std::cout << i+1 << ": " << rates_old[i] << "; ";
        }
        std::cout << std::endl << "res. of DP: ";
        for (il::int_t i = 0; i < nclusters; i++) {
          std::cout << i+1 << ": " << res_v[i] << "; ";
        }
        std::cout << "error: " << err_p;
        std::cout << std::endl << "new fluxes: ";
        for (il::int_t i = 0; i < nclusters; i++) {
          std::cout << i+1 << ": " << rates_cur[i] << "; ";
        }
        std::cout << "error: " << err_f << std::endl;
      }

      // resume...
      rates_old = rates_cur;
    }
  }

  if (accept) {
    if (err_f < coupl_tolerance) {
      if (!mute) {
        std::cout << "-------" << std::endl;
        std::cout << "Well - HFs coupling converged after " << k << " its"
                  << std::endl;
      }
    } else {
      if (!mute) {
        std::cout << "-------" << std::endl;
        std::cout << "Well - HFs coupling NOT converged after " << k
                  << " its; error is: " << err_f << std::endl;
      }
    }

    // run it one more time w. updated fluxes...
    K = K_pre;
    for (il::int_t i = 0; i < nclusters; i++) {
      rates_per_height[i] = rates_cur[i] / frac_height;
    }
    // rate entering the fracture = rate escaping the well / fracture height
    frac_sources_k.setInjectionRates(rates_per_height);
    well_sources_k.setInjectionRates(rates_cur);
    // solve for wellbore flow
    wellSol_k = hfp2d::wellFlowSolverP0(wellSol_n, w_mesh, w_inj,
                                        well_sources_k, ffChurchill, dt,
                                        well_solver_p, true, fracfluid);
    // solve for fracture propagation with given flux
    fracSol_k = hfp2d::FractureFrontLoop(fracSol_n, fracfluid, rock,
                                         frac_sources_k, frac_height, dt,
                                         frac_solver_p, true, il::io, K);
    // would need to have some checks for convergences of the solvers....
    if (fracSol_k.errFront() >= frac_solver_p.frac_front_tolerance) {
      accept = false;
    }

    hfp2d::MultiFracsSolution newSol(fracSol_k, wellSol_k, frac_sources_k,
                                     well_sources_k, frac_height, k, dpc_k,
                                     err_f, err_p);
    return newSol;
  } else {
    if (!mute) {
      std::cout << "The step must be rejected; FF loop NOT converged. "
                // << std::endl
                << "Try to reduce the time step" << std::endl;
    }

    hfp2d::MultiFracsSolution newSol(fracSol_k, wellSol_k, frac_sources_k,
                                     well_sources_k, frac_height, k, dpc_k,
                                     err_f, err_p);
    return newSol;
//    return Sol_n;
  }
}
}
