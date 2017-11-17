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

#include <src/Core/ElasticProperties.h>
#include <src/Core/FluidProperties.h>
#include <src/Core/Mesh.h>

#include <src/Elasticity/AssemblyDDM.h>
#include <src/Elasticity/Simplified3D.h>
#include <src/FractureFluidFlow/ReynoldsP0.h>
#include <src/Core/SimulationParameters.h>
#include <src/Tip/TipAsymptote.h>

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
  hfp2d::Solution  Soln =
      hfp2d::Solution(mesh, 0., width, sheardd, pf_o, sig_o, tau_o);

  // fluid Properties.
  hfp2d::FluidProperties water(1., 0.1, 1.e-10);

  // create source obj. - hardcoded for now....
  il::Array<il::int_t> elt_source{2, 0};
  elt_source[0] = 1;
  elt_source[1] = 4;
  il::Array<double> Qo{2, 0.001};
  hfp2d::Sources the_source = Sources(elt_source, Qo);
  //  std::cout << elt_source[0] << " " << Qo[0] << "\n";

  // create rock Properties obj
  il::Array<double> wh_o{1, 1.e-6}, toughness{1, 1.e6}, Carter{1, 0.};
  hfp2d::SolidProperties the_rock =
      SolidProperties(myelas, toughness, wh_o, Carter);

  // call to Reynolds
  double dt = 0.01;

  hfp2d::SimulationParameters SimulParam;

  bool imp_tip = false;  //
  il::Array<il::int_t> tip_elt{4};
  il::Array<double> tip_width{4, 0.16};
  tip_elt = mesh.tip_elts();

// create a function to get ribbon elements.

  il::Array<il::int_t> ribbon=mesh.getRibbonElements();
  il::Array<double> ribbon_widths{ribbon.size(),0.};

  hfp2d::Solution Soln1 =
      ReynoldsSolverP0(Soln, K, water, the_rock, the_source, dt, imp_tip,
                       tip_elt, tip_width, SimulParam);

  std::cout << " w Tip est prop" <<  sqrt(32./il::pi)*sqrt(h/2.)*the_rock.KIc(0)/myelas.Ep() << "\n";
  std::cout << "ribbon w n " << Soln.openingDD()[ribbon[0]] << "\n";
  std::cout << "ribbon w  n+1 " << Soln1.openingDD()[ribbon[0]] << "\n";


  // get ribbon cell width
   for (il::int_t i=0; i< ribbon.size();i++){
      ribbon_widths[i] = Soln1.openingDD()[ribbon[i]];
   }

//  Tip::TipParameters tipstruct;
//  tipstruct.e_p = the_rock.ElasticProperties().Ep();
//  tipstruct.k1c= the_rock.KIc(0);
//  tipstruct.mu=water.fluidViscosity();
//  tipstruct.cl=the_rock.Cl(0);
//  tipstruct.ta = 1000.;
//  tipstruct.wa= ribbon_widths[0];
//  tipstruct.st=h/2.;

//  bool prop=Tip::isPropagating(tipstruct);
//  std::cout << " is propagating ? "<< prop << "\n";
//
//  Tip::TAInParam tipstruct_a;
//  tipstruct_a.taPrev=tipstruct;
//  tipstruct_a.wa=ribbon_widths[0];
//  tipstruct_a.dt=dt;
//  std::cout << " old distance "<< tipstruct.st <<"\n";
//  Tip::TipParameters tipout = Tip::propagateTip(Tip::res_g_0_s,tipstruct_a,1.e-3,40,false);

//  std::cout << " new distance "<< tipout.st <<"\n";

  hfp2d::Mesh newmesh = mesh;

  newmesh.AddNTipElements(0,0,1,0.);

  std::cout << " new mesh "<<  newmesh.numberOfElements() << "\n";

  std::cout <<" old mesh " << mesh.numberOfElements() << "\n";

  std::cout <<" old Tip " << tip_elt[0] << " " << tip_elt[1] << "\n";

  tip_elt.resize(tip_elt.size()+2);
  std::cout <<" new Tip " << tip_elt[0] << " " << tip_elt[1] << " " << tip_elt[5] << "\n";
  tip_elt.resize(tip_elt.size()-2);
  std::cout <<" back Tip " << tip_elt[0] << " " << tip_elt[1]  << "\n";

//  il::Array2D<double> tt=Soln1.TipsLocation();

  // loop on each tips
  // invert Tip asymptotics for each Tip....
  // check if elements have to be added
  // compute Tip width at tn+1 according to asymptote for these Tip elements.

  // if yes, new mesh
  // pad with zero the width and fluid pressure array of the previous time solution accordingly
  // adapt elasticity matrix


//  Soln=std::move(Soln1);
//  std::cout << "  size tt" << tt.size(1) << "\n";

  std::cout << "now out of reynolds"
            << "\n";

  return 0;
};




hfp2d::Solution FractureFrontLoop(hfp2d::Solution &Sol_n,
                                     il::Array2D<double> &ElasMat,
                                     hfp2d::FluidProperties &fluid,
                                     hfp2d::SolidProperties &rock,
                                     hfp2d::Sources &source, double timestep,
                                     hfp2d::SimulationParameters &simulParams){

// prepare things ....

// notation
// _n : cvged at time tn
// _k : current iteration
// _k_1 : previous iteration

  hfp2d::Mesh mesh_n= Sol_n.CurrentMesh();
  hfp2d::Mesh mesh_k= Sol_n.CurrentMesh();  // might evolve during the iterative process.

  hfp2d::Solution  Sol_n_k = Sol_n; // Solution at time tn
  // might evolve during the iterative process (just because of element number at tn+1 might increase)

  hfp2d::Solution  Soln1_k; // Solution at time tn+1

  // preparation for the LHFM Tip inversion for propagation
  il::Array<il::int_t> ribbon_elt=mesh_n.getRibbonElements(); // fix during this time step
  il::Array<double> s_o{ribbon_elt.size(),0.};
// get initial distance ribbon-cell Tip... s_o (should be stored in solution )
  s_o=Sol_n.RibbonsDistance();

  double ribbon_width,h_ribbon; // for ribbon width and ribbon elt size.

  bool imp_tip = false;
  il::Array<il::int_t>  tip_nodes_n=mesh_n.tip_nodes();
  il::Array<il::int_t>  tip_elt_n = mesh_n.tip_elts();// size might evolve during the iterative process.,

  il::Array<il::int_t>  tip_elt_k = mesh_n.tip_elts();// size might evolve during the iterative process.,
  il::int_t ntip_elt_k = tip_elt_k.size();
  il::Array<double> tip_width_k{tip_elt_k.size(),0.};

//  Tip::TipParameters tipstruct;
//  tipstruct.e_p = rock.ElasticProperties().Ep();
//  tipstruct.k1c= rock.KIc(0);
//  tipstruct.mu=fluid.fluidViscosity();
//  tipstruct.cl=rock.Cl(0);
//  tipstruct.ta = 1000.;

  il::int_t  k=0;
  double errorF=1.;
  while ( (errorF>simulParams.Frac_Front_tolerance) && (k<simulParams.Frac_Front_max_its)){
    k++;

    // solution of Reynolds.
    Soln1_k=ReynoldsSolverP0(Sol_n_k, ElasMat, fluid, rock, source, timestep, imp_tip,
                             tip_elt_k, tip_width_k, simulParams);

    // Invert Tip asymptotes, for all tips ...
    //
    // Tip-regions array, and corresponding Tip width...
    tip_elt_k = mesh_n.tip_elts(); // restart from previous time step value always
    ntip_elt_k = tip_elt_k.size();
    tip_width_k.resize(ntip_elt_k);
    mesh_k =mesh_n; // always restart from the  mesh at tn
    il::int_t n_elt_k=mesh_n.numberOfElements();

    for (il::int_t i=0; i< ribbon_elt.size();i++){

      ribbon_width = Soln1_k.openingDD()[ribbon_elt[i]];
      h_ribbon=mesh_n.elt_size(ribbon_elt[i]);
////      tipstruct.st=s_o[i];
////      tip::TAInParam tipstruct_a;
////      tipstruct_a.taPrev=tipstruct;
////      tipstruct_a.wa=ribbon_width;
////      tipstruct_a.dt=timestep;
////
////      // invert tip asymptote
////      tip::TipParameters tipout = tip::propagateTip(tip::res_g_0_s,tipstruct_a,1.e-5,40,true);
////
////      // add   element in tip regions if needed (always start from tip_elt_n)
//      il::int_t n_add=0; // a priori only 1 Tip elt.
//      if ((tipout.st-h_ribbon/2.)>h_ribbon) {
//        n_add = std::round((tipout.st-h_ribbon/2.)/h_ribbon);
//        ntip_elt_k = tip_elt_k.size(); // store
//        tip_elt_k.resize(ntip_elt_k+n_add); // add space.
//        tip_width_k.resize(ntip_elt_k+n_add);
//        for (il::int_t et=0;et<n_add;et++){
//          tip_elt_k[ntip_elt_k+et]=n_elt_k+et;
//        };
//        // modify mesh accordingly, always restart from previous mesh
//        mesh_k.AddNTipElements(tip_elt_n[i],tip_nodes_n[i],n_add,0.);
//        n_elt_k=mesh_k.numberOfElements(); // don t forget to update.
//      };
//      // compute corresponding tip_volume to impose and thus tip_widths
//      il::Array<double> tipVol{n_add+1,0.},tipW{n_add+1,0.};
//      double sc;
//      for (il::int_t et=0;et<n_add+1;et++){
//         sc= (tipout.st-h_ribbon/2.)+h_ribbon*et;
//         tipVol[et]=Tip::moment0(sc,tipout);
//      };
//      // be careful, do things in reverse from the new Tip location...
//      for (il::int_t et=n_add;et>=0;et--){
//        if (et==n_add){
//          tipW[et]=tipVol[et]/h_ribbon;
//        } else
//        {
//          tipW[et]=(tipVol[et]-tipVol[et+1])/h_ribbon;
//        }
//      };
//      //create the proper tip_width_k contribution for that Tip
//      tip_width_k[i]=tipW[0]; // the pre-existing Tip-elt.
//      for (il::int_t et=0;et<n_add;et++){
//        tip_width_k[ntip_elt_k+et]=tipW[1+et];
//      }
//      // if n_add>=1 - modify Sol_n_k (append elt with zero dds, and pressure)

    } // end of loops on all tips.

    // if number of elements have changed modify the Elasticity_Matrix -> increase etc.
    // always resize to original at tn (in case of element decrease)


    imp_tip = true;

  };

  return Sol_n_k;

}


}