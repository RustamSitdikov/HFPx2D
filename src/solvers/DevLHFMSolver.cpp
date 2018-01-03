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

#include <src/core/ElasticProperties.h>
#include <src/core/Fluid.h>
#include <src/core/Mesh.h>
#include <src/core/DomainMesh.h>

#include <src/core/SimulationParameters.h>
#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/Simplified3D.h>
#include <src/elhsolvers/ReynoldsP0.h>
#include <src/solvers/DevLHFMSolver.h>

#include <src/tip/tipAsymptote.h>


namespace hfp2d {

int TwoParallelHFs(int nelts, double dist) {
  // Test for Reynolds solver
  // 2 fractures parallel separated by dist
  // nelts per frac (numberOfElts should be odd for symmetry)

  // for now for debug we hardcode numberOfElts = 3  so injection is in elt
  // 1 and elt
  // 4

  // step 1 create mesh
  int p = 0;
  double h = 2. / (nelts);  //  element size

  // il::Array<double> x{numberOfElts + 1}; // Not needed
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

  hfp2d::Mesh mesh(xy, myconn, p);

  // background mesh 2 Quads
   il::Array2D<double> nodesB{6,2,0.};
    nodesB(0,0)=-100.;
    nodesB(0,1)=-100.;
    nodesB(1,0)=-100.;
    nodesB(1,1)=100.;
    nodesB(2,0)=0.;
    nodesB(2,1)=100.;
    nodesB(3,0)=0.;
    nodesB(3,1)=-100.;
    nodesB(4,0)=100.;
    nodesB(4,1)=-100.;
    nodesB(5,0)=100.;
    nodesB(5,1)=100.;
  il::Array2D<il::int_t> connB{2,4,0};
  connB(0,0)=0;
  connB(0,1)=1;
  connB(0,2)=2;
  connB(0,3)=3;
  connB(1,0)=2;
  connB(1,1)=3;
  connB(1,2)=4;
  connB(1,3)=5;

  il::Array<il::int_t> matidB{2,0};
  matidB[1]=1;

  hfp2d::DomainMesh bckmesh(nodesB,connB,matidB);

  il::Array<double> myxt(2,0.);
  myxt[0]=1.;myxt[1]=0.;

  il::int_t kk = bckmesh.locate(myxt);

  auto node_CONN = mesh.nodeEltConnectivity();

  // Parameters

  hfp2d::ElasticProperties myelas(10.e9, 0.);
  // initial stress field uniform to test.
  il::Array<double> sig_o{Ntot, 1e6}, tau_o{Ntot, 0.};

  // initial fluid pressure
  il::Array<double> pf_o{Ntot, 1.2e6};  // slightly above sig_o to have initial width
  // create source obj. - hardcoded for now....
  il::Array<double> Qo{2, 0.0001};

  il::Array<il::int_t> elt_source{2, 0};
  elt_source[0] = (nelts - 1) / 2;
  elt_source[1] = nelts + (nelts - 1) / 2;

  hfp2d::Sources the_source = Sources(elt_source, Qo);
  //  std::cout << elt_source[0] << " " << Qo[0] << "\n";

  // fluid properties.
  double mu = 0.1;
  hfp2d::Fluid water(1., mu, 5.e-10);
  // create rock properties obj
  double k1c = 0.5e6;
  il::Array<double> wh_o{1, 1.e-9}, toughness{1, k1c}, Carter{1, 0.};
  hfp2d::SolidProperties the_rock =
      SolidProperties(myelas, toughness, wh_o, Carter);

  double frac_heigth = 1000.;

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
      frac_heigth);  // large pseudo-heigth to reproduce plane-strain kernel

  AddTipCorrectionP0(mesh, myelas, 0, K);
  AddTipCorrectionP0(mesh, myelas, nelts - 1, K);

  AddTipCorrectionP0(mesh, myelas, nelts, K);
  AddTipCorrectionP0(mesh, myelas, Ntot - 1, K);

  // solve the initial elastic system
  il::Array<double> fini{2 * Ntot, 0.};
  il::int_t j = 0;
  std::cout << "e !!!" << (Ntot + Ntot) << "\n";
  for (il::int_t i = 0; i < (Ntot + Ntot); i = i + 2) {
    fini[i] = tau_o[j];
    fini[i + 1] = (pf_o[j] - sig_o[j]);
    j++;
  }

  il::int_t ea = 2;

  std::cout << "elt size" << mesh.eltSize(ea) << "w " << h << "\n";

  std::cout << " size of K" << K.size(0) << " by " << K.size(1) << "\n";
  std::cout << " size of f" << fini.size() << "\n";

  il::Status status;
  // use a direct solver
  il::Array<double> dd_ini = il::linearSolve(K, fini, il::io, status);
  status.ok();
  std::cout << "solved \n";

  double ki=(myelas.Ep()*dd_ini[3]/std::sqrt(3*h/2.)/std::sqrt(32./il::pi)) ;


  //  for (il::int_t i = 0; i < 2 * Ntot; i++) {
  //    std::cout << " dd " << dd_ini[i] << "\n";
  //  }

  il::Array<double> width{mesh.numberOfElts(), 0.},
      sheardd{mesh.numberOfElts(), 0.};
  for (il::int_t i = 0; i < mesh.numberOfElts(); i++) {
    sheardd[i] = dd_ini[2 * i];
    width[i] = dd_ini[2 * i + 1];
  }

  // create a solution at time t=0 object.
  hfp2d::Solution Soln =
      hfp2d::Solution(mesh, 0., width, sheardd, pf_o, sig_o, tau_o);

  il::Array<double> s0{4, 3 * h / 2.};

  il::Array<double> vel0{4, 0.};

  Soln.setTipsVelocity(vel0);
  Soln.setRibbonDistances(s0);

  hfp2d::SimulationParameters SimulParam;

  SimulParam.frac_front_max_its = 40;
  SimulParam.frac_front_tolerance = 1.e-3;
  SimulParam.ehl_relaxation = 0.95;
  SimulParam.ehl_tolerance = 1.e-6;

  double dt = 0.05;
  double dt_min = 0.0001;

  il::int_t jt = 0;
  il::int_t nsteps = 300;

  //  il::Array<il::int_t> tip_region_elt_k=mesh.tipElts();
  //  il::Array<double>    tip_region_width_k{4,0.};

  std::string dir = "../Results/";
  std::string basefilename = "KGD-2HF-M-18.5-";
  std::string filename;

  double mean_tip_v;
  double Mbar = std::pow(myelas.Ep(), 3.) * (12 * mu) * Qo[1] /
                (std::pow((std::sqrt(32. / il::pi) * k1c), 4.));

  std::cout << "Dimensionless Viscosity " << Mbar << "\n";

  while (jt < nsteps) {
    jt++;
    //
    //    Solution Soln1=hfp2d::ReynoldsSolverP0(Soln, K, water, the_rock,
    //    the_source,
    //                                         dt, false, tip_region_elt_k,
    //                                         tip_region_width_k,
    //                                         SimulParam,true);

    Solution Soln1 = hfp2d::FractureFrontLoop(
        Soln, K, water, the_rock, the_source, frac_heigth, dt, SimulParam, true);

    // accept time steps ?

    if (Soln1.errFront()<0.01)
    {
      Soln = Soln1;
      std::cout << " steps # " << jt << " time  " << Soln.time() << "\n";
      std::cout << " P at source " << Soln.pressure()[the_source.SourceElt(0)]
                << "\n";
      std::cout << " w at source " << Soln.openingDD()[the_source.SourceElt(0)]
                << "\n";

      //    std::cout << "size of K: " << K.size(0) << " by " << K.size(1) <<
      //    "\n";

      std::cout << "n elts " << Soln.currentMesh().numberOfElts() << "\n";
      std::cout << " ---------\n";
      filename = dir + basefilename + std::to_string(jt) + ".json";
      Soln.writeToFile(filename);

      mean_tip_v = il::norm(Soln.tipsVelocity(), il::Norm::L2);
      if ( (mean_tip_v > 0.0) ) { //&& (Soln.tipsLocation()(1,0)>1.5)
        double dt_new = 1.25  * h / mean_tip_v;
        // modify to more clever ?
        if (dt_new > 3. * dt) {
          dt = 3. * dt;
        } else {
          if (dt_new < dt * 0.9) {
            dt = 0.9 * dt;
          } else {
            dt = dt_new;
          };
        }
      }
    } else {  // reject time step
      std::cout << "Reject time step - non-convergence on fracture fronts \n";
        if (dt/2.>= dt_min){
        dt=dt/2.;
          std::cout << "Reduce time steps. ";
        jt--;
        }else
        {
          std::cout << "Error on frac. front too large with small time steps. - stop simulation \n";
          break;
        }
    }

  }

  std::cout << "now out of Time step loop " << "\n";

  return 0;
};

////////////////////////////////////////////////////////////////////////////////
// Fracture Front Loop - solve for one time step from tn to tn+timestep
hfp2d::Solution FractureFrontLoop(
    const hfp2d::Solution &Sol_n, il::Array2D<double> &ElasMat,
    hfp2d::Fluid &fluid, const hfp2d::SolidProperties &rock,
    const hfp2d::Sources &source, double frac_height, double timestep,
    hfp2d::SimulationParameters &simulParams, bool mute) {
  // INPUTS
  // Sol_n :: solution object containing the solution at time tn
  // ElasMat :: elasticity matrix on the current mesh  (might be modified)
  // fluid :: fluid object containing the fluid properties
  // rock :: solid properties object
  // source :: injection object
  // frac_heigth :: value of the constant fracture heigth
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

  // number of element that is added to the mesh in this time step - array - for
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

  while (((errorF > simulParams.frac_front_tolerance)) &&
         (k < simulParams.frac_front_max_its)) {

    k++;

    // solution of ELH at fixed fracture front
    Soln_p_1_k = hfp2d::ReynoldsSolverP0(Sol_n_k, ElasMat, fluid, rock, source,
                                         timestep, imp_tip, tip_region_elt_k,
                                         tip_region_width_k, simulParams, mute);

    // Invert tip asymptotes, for all tips ...

    // we always restart from the  mesh at tn
    // tip-regions array, and corresponding tip width array evllces...
    // always restart from previous time step mesh
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
      tip::tipInversion(tip::res_u_0_m, tipstruct, 100 * h_ribbon, 1e-6, 200,
                        mute);

      s_t_k[i] = tipstruct.st;
      v_tip_k[i] = tipstruct.vt;
      // add element in tip regions if needed (always start from tip_elt_n)
      n_add_elt_tip[i] = 0;  // only 1 tip elt (by default) , check now if
                             // fracture has grown in other elements
//      double aux = (tipstruct.st - h_ribbon / 2.) - h_ribbon;

      if (((tipstruct.st - h_ribbon / 2.) > h_ribbon) && v_tip_k[i] > 0.) {
        n_add_elt_tip[i] =
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
      //    prepare opening to impose as volume / mesh size.
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
      double fill_f = (s_o_new[i]-h_ribbon/2.) / h_ribbon;

      if (n_add_elt_tip[i] > 0) {
        // in that case,
        //  we know that the corresponding new tip element is just the last one
        // inserted in mesh_k, i.e. the one that has been just added
        ti_e = mesh_k.tipElts(n_tips - 1);
        ti_b = mesh_k.tipNodes(n_tips - 1);
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
      // get tip nodes b, and the other node a of the tip element (see in mesh)
      // get their coordinates,
      // do   b - (b - a)*(1 - fill_f) to get the tip position.
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
            //region  has changed compared to the previous
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
          sig0[i] = sig0[0];  // constant only -> here would need to call the
                              // initial stress fction
          tau0[i] = tau0[0];
          Pn[i] = 0.;
        }
        // From initial mesh
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
        } else {  // cannot have less element than in the mesh of the solution
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
  std::cout << "end of fracture front loop iterations after " << k
            << " its, errorF=" << errorF << " new number of elts "
            << mesh_k.numberOfElts() << "\n";

  std::cout << " n elts" << Sol_n_k.currentMesh().numberOfElts() << "\n";

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