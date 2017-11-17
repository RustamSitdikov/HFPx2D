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
#include <src/tip/tipAsymptote.h>

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

  hfp2d::Mesh mesh(xy, myconn, p);

  //  hfp2d::Mesh mesh2(xy, myconn, p);

  hfp2d::ElasticProperties myelas(10.e9, 0.3);

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
  il::Array<double> sig_o{Ntot, 0.1e6}, tau_o{Ntot, 0.};

  // initial fluid pressure
  il::Array<double> pf_o{
      Ntot, 0.1e6 + 0.4597e6};  // slightly above sig_o to have initial width

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
  hfp2d::Solution Soln =
      hfp2d::Solution(mesh, 0., width, sheardd, pf_o, sig_o, tau_o);
  il::Array<double> s0{4,3*h/2.};

  il::Array<double> vel0{4,0.};

  Soln.setTipsVelocity(vel0);
  Soln.setRibbonDistances(s0);

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


  hfp2d::SimulationParameters SimulParam;
  SimulParam.Frac_Front_max_its=20;

   double dt = 0.1;
   il::int_t  jt=0;
  il::int_t nsteps=1;

  while (jt<nsteps){
    jt++;

    Solution Soln1=hfp2d::FractureFrontLoop(Soln, K, water, the_rock,the_source,
                                          1000.,dt,  SimulParam,true);

    Soln=Soln1;
    std::cout << " time  " << Soln.time() << "\n";
    std::cout << "size of K: " << K.size(0) << " by " << K.size(1) << "\n";
    std::cout << "n elts " << Soln.CurrentMesh().numberOfElements()<< "\n" ;


  }

//
//  bool imp_tip = false;
//  il::Array<il::int_t> tip_elt{4};
//  il::Array<double> tip_width{4, 0.16};
//  tip_elt = mesh.tip_elts();
//
//  // create a function to get ribbon elements.
//
//  il::Array<il::int_t> ribbon = mesh.getRibbonElements();
//  il::Array<double> ribbon_widths{ribbon.size(), 0.};
//
//  hfp2d::Solution Soln1 =
//      ReynoldsSolverP0(Soln, K, water, the_rock, the_source, dt, imp_tip,
//                       tip_elt, tip_width, SimulParam,true);
//
//  std::cout << " w tip est prop"
//            << sqrt(32. / il::pi) * sqrt(h / 2.) * the_rock.KIc(0) / myelas.Ep()
//            << "\n";
//  std::cout << "ribbon w n " << Soln.openingDD()[ribbon[0]] << "\n";
//  std::cout << "ribbon w  n+1 " << Soln1.openingDD()[ribbon[0]] << "\n";
//
//  // get ribbon cell width
//  for (il::int_t i = 0; i < ribbon.size(); i++) {
//    ribbon_widths[i] = Soln1.openingDD()[ribbon[i]];
//  }
//
//  tip::TipParameters tipstruct;
//  tipstruct.e_p = the_rock.ElasticProperties().Ep();
//  tipstruct.k1c = the_rock.KIc(0);
//  tipstruct.mu = water.fluidViscosity();
//  tipstruct.cl = the_rock.Cl(0);
//  tipstruct.wa = ribbon_widths[0];
//  tipstruct.s0 = h / 2.;
//  tipstruct.dt = dt;
//
//  bool prop = tip::isPropagating(tipstruct);
//  std::cout << " is propagating ? " << prop << "\n";
//
//  tip::tipInversion(tip::res_u_0_m, tipstruct, 100 * h, 1.e-3, 40, false);
//  std::cout << " old distance " << tipstruct.s0 << "\n";
//  std::cout << " new distance " << tipstruct.st << "\n";
//  std::cout << " vel " << tipstruct.vt << "\n";
//
//  hfp2d::Mesh newmesh = mesh;
//
//  newmesh.AddNTipElements(0, 0, 1, 0.);
//
//  std::cout << " new mesh " << newmesh.numberOfElements() << "\n";
//  std::cout << " old mesh " << mesh.numberOfElements() << "\n";
//
//  std::cout << " old tip " << tip_elt[0] << " " << tip_elt[1] << "\n";
//  std::cout << " old tip " << mesh.tip_nodes(0) << " " <<  mesh.tip_nodes(1) << "\n";
//
//  std::cout << " new tip " << newmesh.tip_elts()[0] << " " << newmesh.tip_elts()[1] << " " << "\n";
//  std::cout << " new tip " << newmesh.tip_nodes(0) << " " <<  newmesh.tip_nodes(1) << "\n";

  //  il::Array2D<double> tt=Soln1.TipsLocation();

  // loop on each tips
  // invert tip asymptotics for each tip....
  // check if elements have to be added
  // compute tip width at tn+1 according to asymptote for these tip elements.

  // if yes, new mesh
  // pad with zero the width and fluid pressure array of the previous time
  // solution accordingly
  // adapt elasticity matrix

  //  Soln=std::move(Soln1);
  //  std::cout << "  size tt" << tt.size(1) << "\n";

  std::cout << "now out of  "
            << "\n";

  return 0;
};

////////////////////////////////////////////////////////////////////////////////
// Fracture Front Loop - solve for one time step from tn to tn+timestep
hfp2d::Solution FractureFrontLoop(
    const  hfp2d::Solution &Sol_n, il::Array2D<double> &ElasMat, hfp2d::Fluid &fluid,
    const hfp2d::SolidProperties &rock,const hfp2d::Sources &source, double frac_height,
    double timestep, hfp2d::SimulationParameters &simulParams, bool mute) {
  // INPUTS
  // Sol_n :: solution object containing the solution at time tn
  // ElasMat :: Elasticity matrix on the current mesh  (might be modified)
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

  hfp2d::Mesh mesh_n = Sol_n.CurrentMesh();
  // might evolve during the iterative process.
  hfp2d::Mesh mesh_k = Sol_n.CurrentMesh();

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

  il::Array<il::int_t> ribbon_elt =
      mesh_n.getRibbonElements();  // fix during this time step

  // get initial distance ribbon-cell tip... s_o
  il::Array<double> s_o = Sol_n.RibbonsDistance();
  //
  // array of ribbon tip distance will change during its
  il::Array<double> s_t_k = Sol_n.RibbonsDistance();
  // array of ribbon tip distance during its at previous its.
  il::Array<double> s_t_k_1 = Sol_n.RibbonsDistance();

  // new tip velocity for all tips - // should be initialized at previous timestep
  il::Array<double> v_tip_k = Sol_n.TipsVelocity();
  // the new s_o for the new ribbon cell
  // will be stored for next time step
  il::Array<double> s_o_new = Sol_n.RibbonsDistance();

  il::Array<il::int_t> tip_nodes_n =
      mesh_n.tip_nodes();  // previous time step tip nodes
  il::Array<il::int_t> tip_elt_n = mesh_n.tip_elts();  // and element
  // number of tips that can propagate during the time step
  il::int_t n_tips = tip_elt_n.size();

  // only the tip elements - fixed size, content may evolve during iterations
  il::Array<il::int_t> tips_elt_k = mesh_n.tip_elts();

  il::Array<il::int_t> tip_region_elt_k =
      mesh_n.tip_elts();  // size might evolve during the iterative process.,
  il::int_t ntip_r_elt_k = tip_region_elt_k.size();
  il::int_t ntip_r_elt_k_1 = tip_region_elt_k.size();

  il::Array2D<double> tip_Loc_k{n_tips, 2, 0.};  // new tip location

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

  while ((errorF > simulParams.Frac_Front_tolerance) &&
         (k < simulParams.Frac_Front_max_its)) {
    k++;


    // solution of ELH at fixed fracture front
    Soln_p_1_k = hfp2d::ReynoldsSolverP0(Sol_n_k, ElasMat, fluid, rock, source,
                                  timestep, imp_tip, tip_region_elt_k,
                                  tip_region_width_k, simulParams,mute);

    // Invert tip asymptotes, for all tips ...
    //
    // we always restart from the  mesh at tn
    // tip-regions array, and corresponding tip width array evllces...
    // always restart from previous time step mesh
    mesh_k = mesh_n;
    tip_region_elt_k = mesh_n.tip_elts();
    ntip_r_elt_k = tip_region_elt_k.size();
    tip_region_width_k.resize(ntip_r_elt_k);

    il::int_t n_elt_k = mesh_n.numberOfElements();  // this may increase

    // loop on the different tips / ribbon elts
    for (il::int_t i = 0; i < n_tips; i++) {

      ribbon_width = Soln_p_1_k.openingDD()[ribbon_elt[i]];
      h_ribbon = mesh_n.elt_size(ribbon_elt[i]);
      tipstruct.s0 = s_o[i];
      tipstruct.wa = ribbon_width;
      std::cout << "ribbon opg " << i << " = " << ribbon_width << "is prop ? "<< tip::isPropagating(s_o[i], tipstruct) << "\n";
      // invert tip asymptote for that ribbon elt
      tip::tipInversion(tip::res_u_0_m, tipstruct, 100 * h_ribbon, 1e-4, 100,
                        false);

      s_t_k[i] = tipstruct.st;
      v_tip_k[i] = tipstruct.vt;
      // add element in tip regions if needed (always start from tip_elt_n)
      n_add_elt_tip[i] = 0;  // only 1 tip elt (by default) , check now if
                             // fracture has grown in other elements
      double aux = (tipstruct.st - h_ribbon / 2.) - h_ribbon;

      if (((tipstruct.st - h_ribbon / 2.) > h_ribbon ) && v_tip_k[i]>0.  ) {
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
        mesh_k.AddNTipElements(tip_elt_n[i], tip_nodes_n[i], n_add_elt_tip[i],
                               0.);
        n_elt_k = mesh_k.numberOfElements();//  update the total size
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
      il::int_t ti_e = mesh_n.tip_elts(i);
      il::int_t ti_b = mesh_n.tip_nodes(i);
      s_o_new[i] = s_t_k[i] - h_ribbon / 2. - h_ribbon * n_add_elt_tip[i];
      double fill_f = s_o_new[i] / h_ribbon;

      if (n_add_elt_tip[i] > 0) {
        //  we know that the corresponding new tip element is just the last one
        // inserted in mesh_k, i.e. the one that has been just added
         ti_e = mesh_k.tip_elts(n_tips-1);
          ti_b = mesh_k.tip_nodes(n_tips-1);
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

    }  // end of loops on all tips.

    // ---
    // we now re-compute part of the Elastic stiffness matrix (if necessary)
    // current number of elts in tip regions
    ntip_r_elt_k = tip_region_elt_k.size();
    n_elt_k = mesh_k.numberOfElements(); // update this - effect of last tip reg.

    if (ntip_r_elt_k_1 != ntip_r_elt_k) {  // if the number of element have
                                           // changed compared to the previous
                                           // iteration
      // we resize the vector of the solution at the previous time step ...
      Wn.resize(n_elt_k);
      Pn.resize(n_elt_k);
      Vn.resize(n_elt_k);
      sig0.resize(n_elt_k);
      tau0.resize(n_elt_k);

      // if the number of elements have increased compared to the solution at
      // time tn
      // note that this include also possibly recession between subsequent
      // iterations

      if (mesh_k.numberOfElements() > mesh_n.numberOfElements()) {
        for (il::int_t i = mesh_n.numberOfElements();
             i < mesh_k.numberOfElements(); i++) {
          Wn[i] = 0.;
          Vn[i] = 0.;
          sig0[i] = sig0[0];  // constant only -> here would need to call the
                              // initial stress fction
          tau0[i] = tau0[0];
          Pn[i] = 0.;
        }
        // From initial mesh
        il::int_t ntot_add =
            mesh_k.numberOfElements() - mesh_n.numberOfElements();

        // resize elasticity matrix to previous time-step dofs size
        ElasMat.resize(mesh_n.numberOfDDDofs(), mesh_n.numberOfDDDofs());

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
            ElasMat);  // need to input fracHeight here ...
      } else
      // case where the number of elements recede from the previous iteration
      // to the number of elements at the previous time step
      {
        if (mesh_k.numberOfElements() == mesh_n.numberOfElements()) {
          ElasMat.resize(mesh_n.numberOfDDDofs(), mesh_n.numberOfDDDofs());
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
      tips_elt_k = mesh_k.tip_elts();
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
      errorF += std::abs(1. - s_t_k_1[i] / s_t_k[i]);  // use L1 norm
    }
    s_t_k_1 = s_t_k;  // switch new estimate into previous estimate for the next iterations

  };

  // end of iteration of fracture front position
  std::cout << "end of fracture front loop iterations after " << k
            << " its, errorF=" << errorF <<   " new number of elts " << mesh_k.numberOfElements() << "\n";

  // update last pieces of the solution object
  Soln_p_1_k.setRibbonDistances(s_t_k);  //

  Soln_p_1_k.setTipsVelocity(v_tip_k);

  Soln_p_1_k.setTipsLocation(tip_Loc_k);  //

  Soln_p_1_k.setErrorFront(errorF);
  Soln_p_1_k.setItsFront(k);

  //
  return Soln_p_1_k;
}

}