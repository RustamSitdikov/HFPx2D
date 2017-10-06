//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <il/Array.h>

#include <il/Array2D.h>

#include "DevLHFMSolver.h"


#include <src/core/Mesh.h>
#include <src/core/ElasticProperties.h>
#include <src/core/DOF_Handles.h>
#include <src/core/Fluid.h>

#include <src/Elasticity/AssemblyDDM.h>
#include <src/Elasticity/Simplified3D.h>

namespace hfp2d{

int TwoParallelHFs(int nelts, double dist){

  // Test for Reynolds solver
  // 2 fractures parallel separated by dist
  // nelts per frac (nelts should be odd for symmetry)

  // for now for debug we hardcode nelts = 9  so injection is in elt 5 and elt 14

  // step 1 create mesh
  int p = 0;
  double h = 2. / (nelts);  //  element size

  // il::Array<double> x{nelts + 1}; // Not needed
  int Ntot=2*nelts;

  il::Array2D<double> xy{Ntot + 2, 2, 0.0};
  il::Array2D<il::int_t> myconn{Ntot, 2, 0};
  il::Array2D<il::int_t> id_DD{Ntot, 2*(p+1), 0};
  il::Array2D<il::int_t> id_press{Ntot,p+1, 0};
  il::Array<il::int_t> fracID {Ntot,1};
  il::Array<il::int_t> matID {Ntot,1};
  il::Array<il::int_t> condID {Ntot,1}; // not needed

  double Ep = 1.;                    // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh .... first fracture
  for (il::int_t i = 0; i < nelts+1; ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };
  for (il::int_t i = nelts+1; i < Ntot+2; ++i) {
    xy(i, 0) = -1. + (i-(nelts+1)) * h;
    xy(i, 1) = dist;
  };


  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  for (il::int_t i = nelts ; i < Ntot; ++i) {
    myconn(i, 0) = i+1;
    myconn(i, 1) = i + 2;
  };

  for (il::int_t i=0; i < Ntot; i++){
    for (il::int_t j=0; j < 2*(p+1); j++){
      id_DD(i,j)=i * 2 * (p+1) + j;
    }
  }

  for (il::int_t i=0; i < Ntot; i++){
    id_press(i,0)=i;
  }

// second fracture
  for (il::int_t i=nelts; i < Ntot; i++){
    fracID[i]=2;
  }


  hfp2d::Mesh mesh(p,xy,myconn,id_DD,id_press,fracID,matID,condID);

  il::int_t ndof = mesh.numberOfDisplDofs();

  hfp2d::ElasticProperties myelas(1, 0.);

  std::cout << "EP :" << myelas.Ep() << "\n";

  il::Array2D<double> K{ndof, ndof};
  K = hfp2d::basic_assembly(mesh, myelas,
                            hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
                            1000.);  // large pseudo-heigth to reproduce plane-strain kernel

  AddTipCorrectionP0(mesh, myelas,0, K);
  AddTipCorrectionP0(mesh, myelas,nelts-1, K);

  AddTipCorrectionP0(mesh, myelas,nelts, K);
  AddTipCorrectionP0(mesh, myelas,Ntot-1, K);

  // RE_ARRANGE THE DOF FOR DD HERE !!!
  K= ReArrangeKP0(mesh,K);

// initial stress field uniform to test.
  il::Array<double>  sig_o{Ntot,0.1}, tau_o{Ntot,0.};

  // initial fluid pressure
  il::Array<double>  pf_o{Ntot,0.1+1.e-1};// slightly above sig_o to have initial width

  // fluid properties.
  hfp2d::Fluid water(1., 0.1, 1.e-10);

  // injection location
  il::StaticArray<int,1> inj_location;
  inj_location[0]=5   ;inj_location[0]=14;  // hardcocded for nelts=9 for now!

  // solve the initial elastic system
  il::Array<double> fini{2*Ntot,1.};
  il::int_t j=0;

  for (il::int_t i=0; i < Ntot; i=i+1){
    fini[i]=tau_o[i];
  }
  for (il::int_t i=Ntot; i < 2*Ntot; i=i+1){
    fini[i]=-(pf_o[j]-sig_o[j]);
    j++;
  }

//  for (il::int_t i=0; i < Ntot; i=i+1){
//
//    fini[j]=tau_o[i];
//    std::cout << fini[j] <<"\n";
//    j++;
//    fini[j]=-(pf_o[i]-sig_o[i]);
//    std::cout << fini[j] <<"\n";
//    j++;
//
//  }

  std::cout << " size of K" << K.size(0) << " by " << K.size(1) <<"\n";
  std::cout << " size of f" << fini.size() <<"\n";

  il::Status status;
  // use a direct solver
  il::Array<double> dd_ini = il::linearSolve(K, fini, il::io, status);
  status.ok();
  std::cout << "solved \n";
  for (il::int_t i=0; i<2*Ntot;i++){
    std::cout << " dd " << dd_ini[i] <<"\n" ;
  }

return 0;
};
}