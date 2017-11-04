//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cmath>
#include <complex>
#include <iostream>
//#include <string>

#include <il/Array.h>
//#include "il/math.h"
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/Timer.h>
#include <il/linear_algebra.h>
#include <il/norm.h>

#include "SimpleElastic.h"
#include "src/Elasticity/Simplified3D.h"
#include "src/Elasticity/AssemblyDDM.h"
#include "src/Elasticity/PlaneStrainInfinite.h"
#include "src/core/Mesh.h"
#include "src/core/DOF_Handles.h"
#include "src/core/ElasticProperties.h"

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
// analytical solution of the griffith-crack (ct pressure)
il::Array<double> griffithcrack(const il::Array<double> &x, double a, double Ep,
                                double sig) {
  double coef = 4. * sig / (Ep);
  il::Array<double> wsol{x.size(), 0.};

  for (int i = 0; i < x.size(); ++i) {
    if (std::abs(x[i]) < a) {
      wsol[i] = coef * sqrt(pow(a, 2) - pow(x[i], 2));
    }
  }
  return wsol;
}

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleLinearElement(int nelts) {
  int p = 1;
  double h = 2. / (nelts);  //  element size

  // il::Array<double> x{nelts + 1}; // Not needed

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  il::Array2D<il::int_t> id_displ{nelts, 2*(p+1), 0};
  il::Array2D<il::int_t> id_press{nelts, 2, 0};
  il::Array<il::int_t> fracID {nelts,1};
  il::Array<il::int_t> matID {nelts,1};
  il::Array<il::int_t> condID {nelts,1};

  //int ndof = (nelts) * (p + 1) * 2;  // number of dofs
  double Ep = 1.;                    // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  for (il::int_t i=0; i < nelts; i++){
    for (il::int_t j=0; j < 2*(p+1); j++){
      id_displ(i,j)=i * 2 * (p+1) + j;
    }
  }

  for (il::int_t i=0; i < nelts; i++){
    id_press(i,0)=i;
    id_press(i,1)=i+1;
  }


  //il::Array<il::int_t> matid{nelts, 1};
  // create mesh object
  //hfp2d::Mesh mesh;

  //mesh.init1DMesh(xy, myconn, matid);

  hfp2d::Mesh mesh(p,xy,myconn,id_displ,id_press,fracID,matID,condID);

  il::int_t ndof = mesh.numDDDofs();

  hfp2d::ElasticProperties myelas(1, 0.);

  std::cout << "EP :" << myelas.Ep() << "\n";

  il::Array2D<int> id = hfp2d::dofhandle_dp(2, nelts, p, il::io);  // dof handle for DDs

  // some definitions needed for matrix assembly
  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  //  SegmentData mysege,mysegc;


  il::Array2D<double> K{ndof, ndof};

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
            << (nelts) * (p + 1) * 2 << "---" << mesh.numDDDofs() << "\n";
  std::cout << myconn.size(0) << "\n";
;
  std::cout << "------\n";
  std::time_t result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  il::Timer timer{};
  timer.start();

  K = hfp2d::basic_assembly_new(mesh, myelas,
                            hfp2d::normal_shear_stress_kernel_dp1_dd,
                            0.);  // passing p could be avoided here

  timer.stop();

  std::cout << "------ " << timer.elapsed() << "  \n";
  std::cout << "---#---\n";
  result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, -1.};
  // just opening dds - set shear loads to zero
  for (il::int_t i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  //

  ////// Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  il::int_t i = 0;
  // this piece of codes gets 1D mesh of x doubling the nodes of
  // adjacent elements (for comparison with analytical solution)
  for (il::int_t e = 0; e < nelts; ++e) {
    thex[i] = mesh.node(mesh.connectivity(e, 0), 0);
    thex[i + 1] = mesh.node(mesh.connectivity(e, 1), 0);
    i = i + 2;
  }

  wsol = griffithcrack(thex, 1., 1., 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  il::Array<double> rel_err{ndof / 2 - 2};

  double rel_err_norm;
  // do not compute the relative error in x=\pm 1 as it will be infinite
  for (int j = 1; j < ndof / 2 - 1; ++j) {
    rel_err[j - 1] = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    std::cout << "x : " << thex[j] << "..w anal:" << wsol[j]
              << " w num: " << dd[j * 2 + 1] << " rel error: " << rel_err[j - 1]
              << "\n";
  }

  std::cout << " end of Simple Griffith crack example \n";
  status.ok();

  return il::norm(rel_err, il::Norm::L2);
}

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleS3D_P0(int nelts) {
  int p = 0;                // piece wise constant element
  double h = 2. / (nelts);  //  element size

  //il::Array<double> x{nelts + 1}; // Not needed

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  il::Array2D<il::int_t> id_displ{nelts, 2*(p+1), 0};
  il::Array2D<il::int_t> id_press{nelts, 2, 0};
  il::Array<il::int_t> fracID {nelts,1};
  il::Array<il::int_t> matID {nelts,1};
  il::Array<il::int_t> condID {nelts,1};

  //int ndof = (nelts) * (p + 1) * 2;  // total number of dofs
  double Ep = 1.;                    // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  for (il::int_t i=0; i < nelts; i++){
    for (il::int_t j=0; j < 2*(p+1); j++){
      id_displ(i,j)=i*2*(p+1)+j;
    }
  }

  for (il::int_t i=0; i < nelts; i++){
    id_press(i,0)=i;
    id_press(i,1)=i+1;
  }

  //il::Array<il::int_t> matid{nelts, 1};
  // create mesh object
  //hfp2d::Mesh mesh;

  //mesh.init1DMesh(xy, myconn, matid);
  //hfp2d::Mesh mesh(xy,myconn);
  hfp2d::Mesh mesh(p,xy,myconn,id_displ,id_press,fracID,matID,condID);

  il::int_t ndof = mesh.numDDDofs();

  hfp2d::ElasticProperties myelas(1, 0.);
  //  myelas.ElasticProperties(1.,0.);
  std::cout << "EP :" << myelas.Ep() << "\n";

  il::Array2D<int> id = hfp2d::dofhandle_dp(2, nelts, p, il::io);  // dof handle for DDs

  // some definitions needed for matrix assembly
  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  //  SegmentData mysege,mysegc;

  il::Array2D<double> K{ndof, ndof};

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
            << (nelts) * (p + 1) * 2 << "---" << ndof << "\n";
  std::cout << myconn.size(0) << "\n";

  std::cout << "------\n";
  std::time_t result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  il::Timer timer{};
  timer.start();

//  K = hfp2d::basic_assembly(mesh, id, p, myelas,
//                            hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
//                            1000.);  // large pseudo-heigth to reproduce plane-strain kernel
//  // passing p could be avoided here

  K = hfp2d::basic_assembly_new(mesh, myelas,
                            hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
                            1000.);  // large pseudo-heigth to reproduce plane-strain kernel


  timer.stop();
  std::cout << "---end of assembly--- \n";

  std::cout << "------ " << timer.elapsed() << "  \n";
  std::cout << "---#---\n";
  result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, -1.};  // be careful of sign
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  //
  //  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  //SegmentData sege;
  for (int e = 0; e < nelts; ++e) {
    hfp2d::SegmentData sege = mesh.getElementData(e);
    //sege=get_segment_DD_data(mesh,e,p);
    thex[e] = sege.Xmid(0);
  }

  wsol = griffithcrack(thex, 1., 1., 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  il::Array<double> rel_err{ndof / 2 - 2};

  double rel_err_norm;
  // do not compute the relative error in x=\pm 1 as it will be infinite
  for (int j = 1; j < ndof / 2 - 1; ++j) {
    rel_err[j - 1] = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    std::cout << "x : " << thex[j] << "..w anal:" << wsol[j]
              << " w num: " << dd[j * 2 + 1] << " rel error: " << rel_err[j - 1]
              << "\n";
  }

  std::cout << " end of Simple Griffith crack example \n";
  status.ok();

  return il::norm(rel_err, il::Norm::L2);
}
}