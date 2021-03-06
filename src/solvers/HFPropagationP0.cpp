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
#include <il/linear_algebra/dense/norm.h>

#include <src/core/DomainMesh.h>
#include <src/core/ElasticProperties.h>
#include <src/core/Fluid.h>
#include <src/core/Mesh.h>

#include <src/core/SimulationParameters.h>
#include <src/ehlsolvers/ReynoldsP0.h>
#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/Simplified3D.h>
#include <src/input/json/LoadInputMultistage.h>
#include <src/input/json/loadJsonMesh.h>
#include <src/solvers/HFPropagationP0.h>
#include <src/tip/tipAsymptote.h>
#include <src/util/json.hpp>

namespace hfp2d {

using json = nlohmann::json;

int ParallelHFs(std::string &filename) {
  std::string wellfilename = "../Debug/ParallelHFTestsMvertex.json";

  std::ifstream input(wellfilename);  // ?
  json js;
  input >> js;

  if ((js.count("Fractures mesh") != 1)) {
    std::cout << "No Fractures mesh in json input file ";
    il::abort();
  }
  json j_fmesh = js["Fractures mesh"];

  hfp2d::Mesh mesh = loadJsonMesh(j_fmesh);

  std::cout << " Mesh loaded " << mesh.numberOfElts() << "\n";

  if ((js.count("Model parameters") != 1)) {
    std::cout << "No parameters in input file ";
    il::abort();
  }
  json j_params = js["Model parameters"];


  if (j_params.count("Fluid properties") != 1) {
    std::cout << "No fluid properties input in  model parameters";
    il::abort();
  }

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

  if (j_params.count("Number of fractures") != 1) {
    std::cout << "No Number of fractures  input in  model parameters";
    il::abort();
  }
  auto nfracs = j_params["Number of fractures"].get<long>();

  il::int_t nelts = mesh.numberOfElts() / nfracs;
  // we must ensure that the number of elements is odd.
  if (!(nelts % 2)) {
    std::cout << "  Number of elts per fracture is not odd !";
    il::abort();
  }

  if (j_params.count("In-situ conditions") != 1) {
    std::cout << "No In-situ conditions input in  model parameters";
    il::abort();
  }
  json j_insitu = j_params["In-situ conditions"];
  hfp2d::InSituConditions inSituStress = loadInSitu(j_insitu);


  if (j_params.count("Injection rate") != 1) {
    std::cout << "No Injection rate  input in  model parameters";
    il::abort();
  }
  auto Qo = j_params["Injection rate"].get<double>();

  il::Array<double> Rates{nfracs, Qo};

  il::Array<il::int_t> elt_source{nfracs, 0};
  for (il::int_t i = 0; i < nfracs; i++) {
    elt_source[i] = i * nelts + (nelts - 1) / 2;
    std::cout <<  elt_source[i] << " "  << Rates[i] <<"\n";
  }

  hfp2d::Sources the_source = Sources(elt_source, Rates);

  // background wellMesh 2 Quads
  il::Array2D<double> nodesB{6, 2, 0.};
  nodesB(0, 0) = -100.;
  nodesB(0, 1) = -100.;
  nodesB(1, 0) = -100.;
  nodesB(1, 1) = 100.;
  nodesB(2, 0) = 0.;
  nodesB(2, 1) = 100.;
  nodesB(3, 0) = 0.;
  nodesB(3, 1) = -100.;
  nodesB(4, 0) = 100.;
  nodesB(4, 1) = -100.;
  nodesB(5, 0) = 100.;
  nodesB(5, 1) = 100.;
  il::Array2D<il::int_t> connB{2, 4, 0};
  connB(0, 0) = 0;
  connB(0, 1) = 1;
  connB(0, 2) = 2;
  connB(0, 3) = 3;
  connB(1, 0) = 2;
  connB(1, 1) = 3;
  connB(1, 2) = 4;
  connB(1, 3) = 5;

  il::Array<il::int_t> matidB{2, 0};
  matidB[1] = 1;
  hfp2d::DomainMesh bckmesh(nodesB, connB, matidB);
  il::Array<double> myxt(2, 0.);
  myxt[0] = 1.;
  myxt[1] = 0.;
  il::int_t kk = bckmesh.locate(myxt);
  auto node_CONN = mesh.nodeEltConnectivity();

  // Parameters
  il::int_t Ntot = mesh.numberOfElts();

  // create source obj. - hardcoded for now....
  //  std::cout << elt_source[0] << " " << Qo[0] << "\n";

 // build matrix.
  hfp2d::ElasticProperties myelas = rock.ElasticProperties();

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, frac_height);

  // add tip correction for P0 for each tips in the mesh
  for (il::int_t i = 0; i < mesh.tipElts().size(); i++) {
    AddTipCorrectionP0(mesh, myelas, mesh.tipElts(i), K);
  }

  // initial net loading.
  il::Array<double> ftini{mesh.numberDDDofs(), 0.},
      pfo{mesh.numberDDDofs(), 0.};
  ftini = inSituStress.uniformAllInSituTractions(mesh);
  for (il::int_t i = 0; i < nfracs; i++) {
    for (il::int_t j = 0; j < nelts; j++) {
      pfo[i * 2 * nelts + 2 * j + 1] =  ftini[1]*(1.005)  ;
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
  for (il::int_t i = 0; i < mesh.numberDDDofs(); i++) {
    std::cout << "net load " << pnet0[i] << " dd " << dd_ini[i] << "\n";
  }

  il::Array<double> width{mesh.numberOfElts(), 0.},
      sheardd{mesh.numberOfElts(), 0.};

  il::Array<double> sn_o{mesh.numberOfElts(), 0.},
      tau_o{mesh.numberOfElts(), 0.},
      fluid_pressure_o{mesh.numberOfElts(), 0.};

  for (il::int_t i = 0; i < mesh.numberOfElts(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
    sn_o[i] = ftini[2 * i + 1];
    tau_o[i] = ftini[2 * i];
    fluid_pressure_o[i] = pfo[2 * i + 1];
  }

  // create a solution at the current time of the wellsolution  object.
  hfp2d::Solution fracSol_n =
      hfp2d::Solution(mesh, 0., width, sheardd,
                      fluid_pressure_o, sn_o, tau_o);
  il::Array<double> vel0{nfracs*2, 0.};  // initial tip velocity
  fracSol_n.setTipsVelocity(vel0);


  il::Array<double> s0{2 * nfracs, 0.};  // initial ribbont tip distance
  for (il::int_t i = 0; i < mesh.tipElts().size(); i++) {
    il::int_t e = mesh.tipElts(i);
    s0[i] = 1.5 * mesh.eltSize(e);
  }
  fracSol_n.setRibbonDistances(s0);

  il::int_t ea = 0;

  hfp2d::SimulationParameters SimulParam; // hardcoded here ....
  SimulParam.frac_front_max_its = 50;
  SimulParam.frac_front_tolerance = 1.e-3;
  SimulParam.ehl_relaxation = 0.95;
  SimulParam.ehl_tolerance = 1.e-6;

  if ((js.count("Simulation parameters") != 1)) {
    std::cout << "No Simulation parameters in input file ";
    il::abort();
  }
  json j_simul = js["Simulation parameters"];

  double max_time = 1.;
  if ( j_simul.count("Maximum time") ==1 ) {
    max_time = j_simul["Maximum time"].get<double>();
  }
  il::int_t max_steps = 150;
  if ( j_simul.count("Maximum number of steps") ==1 ) {
    max_steps = j_simul["Maximum number of steps"].get<long>();
  }

  double dt = 0.002;
  if ( j_simul.count("Time step") ==1 ) {
    dt = j_simul["Time step"].get<double>();
  }

  double dt_min = 0.00001;
  if ( j_simul.count("Minimum time step") ==1 ) {
    dt = j_simul["Minimum time step"].get<double>();
  }

  std::string basefilename = "HFs_given_rate"; // default name
  if (js.count("Results files name core")==1){
    basefilename =js["Results files name core"].get<std::string>();
  }

  std::string resfilename;

  double mean_tip_v;

  double Mbar = std::pow(myelas.Ep(), 3) * (12. * fracfluid.fluidViscosity()) *
                Qo / (std::pow((std::sqrt(32. / il::pi) * rock.KIc(ea)), 4.));

  std::cout << "Dimensionless Viscosity " << Mbar << "\n";

  il::int_t jt = 0;
  while ( (jt < max_steps)  && (fracSol_n.time() < max_time) ) {
    jt++;
    //
    //    Solution Soln1=hfp2d::ReynoldsSolverP0(fracSol_n, K, water, the_rock,
    //    the_source,
    //                                         dt, false, tip_region_elt_k,
    //                                         tip_region_width_k,
    //                                         SimulParam,true);

    Solution Soln1 =
        hfp2d::FractureFrontLoop(fracSol_n, fracfluid, rock, the_source, frac_height,
                                 dt, SimulParam, true, il::io, K);

    // accept time step if the error is below 0.01 (hardcoded value for now, should be named relaxed_tolerance ?)
    if (Soln1.errFront() < 0.01) {
      fracSol_n = Soln1;

      std::cout << " Steps # " << jt << " Time: " << fracSol_n.time() << " Time step: " << fracSol_n.timestep() << "\n";
      std::cout << " Error on frac front " << fracSol_n.errFront() << " after " << fracSol_n.frontIts() << " iterations" << "\n";
      std::cout << " P at source " << fracSol_n.pressure()[the_source.SourceElt(0)]
                << "\n";
      std::cout << " w at source " << fracSol_n.openingDD()[the_source.SourceElt(0)]
                << "\n";

      std::cout << " n elts " << fracSol_n.currentMesh().numberOfElts() << "\n";
//      std::cout << " nn " << fracSol_n.currentMesh().connectivity().size(0) <<"\n";

      resfilename =   basefilename + std::to_string(jt) + ".json";

      fracSol_n.writeToFile(resfilename);

      // adjust time step
//      for (il::int_t i=0;i< mesh.tipElts().size();i++){
//        std::cout << " tip " << i << " velocity " << fracSol_n.tipsVelocity()[i];
//      }
//      std::cout << "\n";

      std::cout << " ----++++-----++++-------\n";

      mean_tip_v = il::norm(fracSol_n.tipsVelocity(), il::Norm::L2);

      if ((mean_tip_v > 0.0)) {  //&& (fracSol_n.tipsLocation()(1,0)>1.5)
        double dt_new = 1. * mesh.eltSize(ea) / mean_tip_v;
        // modify to more clever ?
        if (dt_new > 2.5 * dt) {
          dt = 2.5 * dt;
        } else {
          if (dt_new < dt * 0.8) {
            dt = 0.8 * dt;
          } else {
            dt = dt_new;
          };
        }
      }

    } else {

      // reject time step
      std::cout << "Reject time step; non-convergence on fracture fronts \n";
      std::cout << " steps # " << jt << " time  " << Soln1.time() << "Time step: " << Soln1.timestep() << "\n";
      std::cout << " Error on frac front " << Soln1.errFront() << " after " << Soln1.frontIts() << " its" << "\n";

      if (dt / 2. >= dt_min) {
        dt = dt / 2.;
        std::cout << "Reducing time step to "
                  << dt << "s. in order to re-try \n";
        jt--;
      } else {
        std::cout << "Error on frac. front too large with small time steps. - "
                     "stop simulation \n";
        break;
      }
    }
  }

  std::cout << "now out of Time step loop "
            << "\n";

  return 0;
};

////////////////////////////////////////////////////////////////////////////////
// for a routine - for well + n fracs benchmark
// with json inputs, see MultiFracsSolver.cpp
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Fracture Front Loop - solve for one time step from tn to tn+timestep
hfp2d::Solution FractureFrontLoop(const hfp2d::Solution &Sol_n,
                                  hfp2d::Fluid &fluid,
                                  const hfp2d::SolidProperties &rock,
                                  const hfp2d::Sources &source,
                                  double frac_height, double timestep,
                                  hfp2d::SimulationParameters &simulParams,
                                  bool mute, il::io_t,
                                  il::Array2D<double> &ElasMat) {
  // INPUTS
  // Sol_n :: solution object containing the solution at time tn
  // ElasMat :: elasticity matrix on the current wellMesh  (might be modified)
  // fluid :: fluid object containing the fluid properties
  // rock :: solid properties object
  // source :: injection object
  // frac_height :: value of the constant fracture height
  // timestep :: double, value of current time step
  // simulParams :: structure containing the solvers parameters

  // notation
  // _n : solution at time tn (beginning of step)
  // _k : current iteration
  // _k_1 : previous iteration

  hfp2d::Mesh mesh_n = Sol_n.currentMesh();
  // might evolve during the iterative process.
  hfp2d::Mesh mesh_k = Sol_n.currentMesh();

  hfp2d::Solution Sol_n_k = Sol_n;  // Solution at time tn
  // might evolve during the iterative process (just because of element number
  // at tn+1 might increase)

  il::Array<double> Wn = Sol_n.openingDD();
  il::Array<double> Vn = Sol_n.shearDD();
  il::Array<double> Pn = Sol_n.pressure();

  il::Array<double> sig0 = Sol_n.sigma0();
  il::Array<double> tau0 = Sol_n.tau0();

  hfp2d::Solution Soln_p_1_k;  // Solution at time tn+1

  // preparation for the LHFM tip inversion for propagation

  // ribbon elts (fix during one time step)
  il::Array<il::int_t> ribbon_elt = mesh_n.getRibbonElements();

  // get initial distance ribbon-cell tip... s_o
  il::Array<double> s_o = Sol_n.ribbonsDistance();
  //
  // array of ribbon elt tip distance will change during its
  il::Array<double> s_t_k = Sol_n.ribbonsDistance();
  // array of ribbon tip distance during its at the previous its.
  il::Array<double> s_t_k_1 = Sol_n.ribbonsDistance();

  // new tip velocity for all tips -
  il::Array<double> v_tip_k = Sol_n.tipsVelocity();
  // the new s_o for the new ribbon cell
  // will be stored for the final solution
  il::Array<double> s_o_new = s_o;

  // previous time step tip nodes
  il::Array<il::int_t> tip_nodes_n = mesh_n.tipNodes();
  il::Array<il::int_t> tip_elt_n = mesh_n.tipElts();  // and element
  // number of tips that can propagate during the time step
  il::int_t n_tips = tip_elt_n.size();

  // only the tip elements - fixed size, content may evolve during iterations
  il::Array<il::int_t> tips_elt_k = mesh_n.tipElts();

  il::Array<il::int_t> tip_region_elt_k =
      mesh_n.tipElts();  // size might evolve during the iterative process.,
  il::int_t ntip_r_elt_k = tip_region_elt_k.size();
  il::int_t ntip_r_elt_k_1 = tip_region_elt_k.size();

  il::Array2D<double> tip_Loc_k{n_tips, 2, 0.};  // new (x,y) tips location

  // number of element that is added to the wellMesh in this time step - array -
  // for
  // each tips
  il::Array<il::int_t> n_add_elt_tip{n_tips, 0};

  double ribbon_width, h_ribbon;  // for ribbon width and ribbon elt size.

  bool imp_tip = false;  // boolean for ELH with imposed or not tip width.

  // an array containing the imposed tip width
  il::Array<double> tip_region_width_k{tip_region_elt_k.size(), 0.};

  hfp2d::ElasticProperties elasrock = rock.ElasticProperties();

  // tip parameters structure
  tip::TipParameters tipstruct;
  tipstruct.e_p = rock.ElasticProperties().Ep();
  tipstruct.k1c = rock.KIc(0);
  tipstruct.mu = fluid.fluidViscosity();
  tipstruct.cl = rock.Cl(0);
  tipstruct.dt = timestep;

  il::int_t k = 0;

  double errorF = 1.;
  bool firstime_elas_mod = true;

  if (!mute) {
    std::cout << "++++ Fracture Front loop ++++ \n";
  }

  while (((errorF > simulParams.frac_front_tolerance)) &&
         (k < simulParams.frac_front_max_its)) {
    k++;

    if (!mute) {
      std::cout << "its " << k <<  " error " << errorF << "n elts " << Sol_n_k.currentMesh().numberOfElts() <<  " \n";
    }

    // solution of ELH at fixed fracture front
    Soln_p_1_k = hfp2d::ReynoldsSolverP0(Sol_n_k, ElasMat, fluid, rock, source,
                                         timestep, imp_tip, tip_region_elt_k,
                                         tip_region_width_k, simulParams, mute);

    // Invert tip asymptotes, for all tips ...

    // we always restart from the  wellMesh at tn
    // tip-regions array, and corresponding tip width array evllces...
    // always restart from previous time step wellMesh
    mesh_k = mesh_n;
    tip_region_elt_k = mesh_n.tipElts();
    ntip_r_elt_k = tip_region_elt_k.size();
    tip_region_width_k.resize(ntip_r_elt_k);

    il::int_t n_elt_k = mesh_n.numberOfElts();  // this may increase

    // loop on the different tips / ribbon elts
    for (il::int_t i = 0; i < n_tips; i++) {
      ribbon_width = Soln_p_1_k.openingDD()[ribbon_elt[i]];
      h_ribbon = mesh_n.eltSize(ribbon_elt[i]);
      tipstruct.s0 = s_o[i];
      tipstruct.wa = ribbon_width;
      //      std::cout << "ribbon opg " << i << " = " << ribbon_width << "is
      //      prop ? "<< tip::isPropagating(s_o[i], tipstruct) << "\n";
      // invert tip asymptote for that ribbon elt
      tip::tipInversion((imf::ImplicitFunction)tip::res_u_0_m,
                        tipstruct, 100 * h_ribbon, 1e-6, 200,
                        mute);

      s_t_k[i] = tipstruct.st;
      v_tip_k[i] = tipstruct.vt;
      // add element in tip regions if needed (always start from tip_elt_n)
      n_add_elt_tip[i] = 0;  // only 1 tip elt (by default) , check now if
                             // fracture has grown in other elements
      //      double aux = (tipstruct.st - h_ribbon / 2.) - h_ribbon;

      if (((tipstruct.st - h_ribbon / 2.) > h_ribbon) && v_tip_k[i] > 0.) {
        n_add_elt_tip[i] =
                (il::int_t)
                        std::floor((tipstruct.st - h_ribbon / 2.) / h_ribbon);
        // get the current number of tips element before addition
        ntip_r_elt_k = tip_region_elt_k.size();
        tip_region_elt_k.resize(ntip_r_elt_k + n_add_elt_tip[i]);  // add space.
        tip_region_width_k.resize(ntip_r_elt_k + n_add_elt_tip[i]);
        for (il::int_t et = 0; et < n_add_elt_tip[i]; et++) {
          tip_region_elt_k[ntip_r_elt_k + et] = n_elt_k + et;
        };
        // add these new elements in the mesh_k
        mesh_k.addNTipElts(tip_elt_n[i], tip_nodes_n[i], n_add_elt_tip[i], 0.);
        n_elt_k = mesh_k.numberOfElts();  //  update the total size
      };

      // compute corresponding tip_volume to impose tip_widths
      il::Array<double> tipVol{n_add_elt_tip[i] + 1, 0.},
          tipW{n_add_elt_tip[i] + 1, 0.};

      double sc;
      // compute the volume from the edge of the element up to the
      // tip.
      for (il::int_t et = 0; et < n_add_elt_tip[i] + 1; et++) {
        sc = (tipstruct.st - h_ribbon / 2.) - h_ribbon * et;
        tipVol[et] = tip::moment0(sc, tipstruct);
      };
      //    prepare opening to impose as volume / wellMesh size.
      // be careful, do things in reverse from the new tip location at the end
      // of the vector...
      for (il::int_t et = n_add_elt_tip[i]; et >= 0; et--) {
        if (et == n_add_elt_tip[i]) {
          tipW[et] = tipVol[et] / h_ribbon;  // near tip element
        } else {
          tipW[et] = (tipVol[et] - tipVol[et + 1]) / h_ribbon;
        }
      };
      // create the proper tip_region_width_k contribution for that tip
      tip_region_width_k[i] = tipW[0];  // the pre-existing tip-elt width
      // if more elts are now present (compared to time tn) for that tip
      // index starting at current value of ntip_r_elt_k
      for (il::int_t et = 0; et < n_add_elt_tip[i]; et++) {
        tip_region_width_k[ntip_r_elt_k + et] = tipW[1 + et];
      }

      // compute  the location of the new current tip location estimate
      // (to store it when cvged directly)...
      il::int_t ti_e = mesh_n.tipElts(i);
      il::int_t ti_b = mesh_n.tipNodes(i);
      s_o_new[i] = s_t_k[i] - h_ribbon * n_add_elt_tip[i];
      double fill_f = (s_o_new[i] - h_ribbon / 2.) / h_ribbon;

      if (n_add_elt_tip[i] > 0) {
        // in that case,
        //  we know that the corresponding new tip element is just the last one
        // inserted in mesh_k, i.e. the one that has been just added
        ti_e = mesh_k.tipElts(i);
        ti_b = mesh_k.tipNodes(i);
      }
      il::int_t ti_a;
      if (mesh_k.connectivity(ti_e, 0) == ti_b) {
        ti_a = mesh_k.connectivity(ti_e, 1);
      } else {
        if (mesh_k.connectivity(ti_e, 1) == ti_b) {
          ti_a = mesh_k.connectivity(ti_e, 0);
        } else {
          il::abort();
        }
      };
      // get tip nodes b, and the other node a of the tip element (see in
      // Mesh)
      // get their coordinates,
      // do   b - (b - a)*(1 - fill_f) to get the tip location.
      il::StaticArray<double, 2> a = mesh_k.coordinates(ti_a);
      il::StaticArray<double, 2> b = mesh_k.coordinates(ti_b);

      tip_Loc_k(i, 0) = b[0] - (b[0] - a[0]) * (1. - fill_f);
      tip_Loc_k(i, 1) = b[1] - (b[1] - a[1]) * (1. - fill_f);
      // tip location done

    }  // end of loops on all tips.

    // ---
    // we now re-compute part of the Elastic stiffness matrix (if necessary)
    // current number of elts in tip regions
    ntip_r_elt_k = tip_region_elt_k.size();
    n_elt_k = mesh_k.numberOfElts();  // update this - effect of last tip reg.

    if (ntip_r_elt_k_1 != ntip_r_elt_k) {  // if the number of element in tip
      // region  has changed compared to the previous
      // iteration  we resize the vector of the solution at the previous
      // time step ...
      Wn.resize(n_elt_k);
      Pn.resize(n_elt_k);
      Vn.resize(n_elt_k);
      sig0.resize(n_elt_k);
      tau0.resize(n_elt_k);

      // if the number of elements have increased compared to the solution at
      // time tn
      // note that this include also possibly recession between subsequent
      // iterations
      if (mesh_k.numberOfElts() > mesh_n.numberOfElts()) {
        for (il::int_t i = mesh_n.numberOfElts(); i < mesh_k.numberOfElts();
             i++) {
          Wn[i] = 0.;
          Vn[i] = 0.;
          sig0[i] = sig0[0];  // constant only ->
          // todo here would need to call the initial stress fction
          tau0[i] = tau0[0];
          Pn[i] = 0.;
        }
        // From initial wellMesh
        il::int_t ntot_add = mesh_k.numberOfElts() - mesh_n.numberOfElts();

        // resize elasticity matrix to previous time-step dofs size
        ElasMat.resize(mesh_n.numberDDDofs(), mesh_n.numberDDDofs());

        // loop to remove tip correction ONLY  if it is the first time that
        // the elasticity matrix has to be modified.
        if (firstime_elas_mod) {
          for (il::int_t i = 0; i < ribbon_elt.size(); i++) {
            hfp2d::RemoveTipCorrectionP0(mesh_n, elasrock, tip_elt_n[i],
                                         ElasMat);
          };
          firstime_elas_mod = false;
        }

        hfp2d::basic_assembly_add_elts(
            mesh_k, ntot_add, elasrock,
            hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, frac_height, il::io,
            ElasMat);
      } else
      // case where the number of elements recede from the previous iteration
      // to the number of elements at the previous time step
      {
        if (mesh_k.numberOfElts() == mesh_n.numberOfElts()) {
          ElasMat.resize(mesh_n.numberDDDofs(), mesh_n.numberDDDofs());
        } else {  // cannot have less element than in the wellMesh of the
                  // solution
                  // at the previous time step
          std::cout << "incorrect - back propagation ! ";
          il::abort();
        }
      }

      // adapted solution at time tn (previous step), padded with zero(s) on the
      // added tip region elements
      Sol_n_k = hfp2d::Solution(mesh_k, Sol_n.time(), Wn, Vn, Pn, sig0, tau0);

      // loop on the actual tip element to add tip correction.
      tips_elt_k = mesh_k.tipElts();
      for (il::int_t i = 0; i < ribbon_elt.size(); i++) {
        hfp2d::AddTipCorrectionP0(mesh_k, elasrock, tips_elt_k[i], ElasMat);
      };
    }
    // current number of elements in tip regions at that current iteration
    // (to check if it will be different at the next iteration)
    ntip_r_elt_k_1 = ntip_r_elt_k;

    // next iterations, impose width in tip regions ....
    imp_tip = true;

    // compute error on fracture front position .... on the distance ribbon-tip
    // for each tip
    errorF = 0.;
    for (il::int_t i = 0; i < n_tips; i++) {
      errorF += std::abs(1. - (s_t_k_1[i]) / (s_t_k[i]));  // use L2 norm
    }
    //    errorF=std::sqrt(errorF);
    s_t_k_1 = s_t_k;  // switch new estimate into previous estimate for the next
                      // iterations
  };

  // last call to Reynolds with the converged position of the front...
  // if it has changed !
  if (errorF != 0) {
    //    std::cout << "last reynolds call \n";
    Soln_p_1_k = hfp2d::ReynoldsSolverP0(Sol_n_k, ElasMat, fluid, rock, source,
                                         timestep, imp_tip, tip_region_elt_k,
                                         tip_region_width_k, simulParams, mute);
  };

  // end of iteration of fracture front position

  if (!mute){
  std::cout << "end of fracture front loop iterations after " << k
            << " its, errorF=" << errorF << " new number of elts "
            << mesh_k.numberOfElts() <<   " error w EHL "<< Soln_p_1_k.errOpening() <<"\n";
  }

  // update last pieces of the solution object...

  Soln_p_1_k.setRibbonDistances(s_o_new);  //

  Soln_p_1_k.setTipsVelocity(v_tip_k);

  Soln_p_1_k.setTipsLocation(tip_Loc_k);  //

  Soln_p_1_k.setErrorFront(errorF);
  Soln_p_1_k.setItsFront(k);

  //
  return Soln_p_1_k;
}
}