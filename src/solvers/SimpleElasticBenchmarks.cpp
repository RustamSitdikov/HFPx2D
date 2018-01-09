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
#include <string>

#include "il/Array.h"
#include "il/math.h"
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/Timer.h>
#include <il/linear_algebra.h>
#include <il/norm.h>

#include "src/Elasticity/AssemblyDDM.h"

#include <src/Elasticity/PlaneStrainInfinite.h>
#include <src/Elasticity/Simplified3D.h>
#include <src/core/Mesh.h>

#include "SimpleElasticBenchmarks.h"
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

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D wellMesh ....
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, p);

  il::int_t ndof = mesh.numberDDDofs();

  hfp2d::ElasticProperties myelas(1, 0.3);

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd,
      0.);  // passing p could be avoided here

  // solve a constant pressurized crack problem...
  il::Array<double> f{
      ndof, 1.};  // here we solve with +, so we get positive DD in opening
  // just opening dds - set shear loads to zero
  for (il::int_t i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  status.ok();

  ////// Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  il::int_t i = 0;
  // this piece of codes gets 1D wellMesh of x doubling the nodes of
  // adjacent elements (for comparison with analytical solution)
  for (il::int_t e = 0; e < nelts; ++e) {
    thex[i] = mesh.coordinates(mesh.connectivity(e, 0), 0);
    thex[i + 1] = mesh.coordinates(mesh.connectivity(e, 1), 0);
    i = i + 2;
  }

  wsol =
      griffithcrack(thex, 1., myelas.Ep(), 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  il::Array<double> rel_err{ndof / 2 - 2};

  double rel_err_norm;
  // do not compute the relative error in x=\pm 1 as it will be infinite
  for (int j = 1; j < ndof / 2 - 1; ++j) {
    rel_err[j - 1] = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    //    std::cout << "x : " << thex[j] << "..w anal:" << wsol[j]
    //              << " w num: " << dd[j * 2 + 1] << " rel error: " <<
    //              rel_err[j - 1]
    //              << "\n";
  }

  std::cout << " end of Simple Griffith crack example P1 \n";

  return il::norm(rel_err, il::Norm::L2);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleLinearElement_byNodes(int nelts) {
  int p = 1;
  double h = 2. / (nelts);  //  element size

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D wellMesh ....
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, p);

  il::int_t ndof = mesh.numberDDDofs();

  hfp2d::ElasticProperties myelas(1, 0.3);

  il::Array2D<double> K = hfp2d::basic_assembly_nodal(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal,
      0.);  // passing p could be avoided here

  // solve a constant pressurized crack problem...
  il::Array<double> f{
      ndof, 1.};  // here we solve with +, so we get positive DD in opening
  // just opening dds - set shear loads to zero
  for (il::int_t i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  status.ok();

  ////// Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  il::int_t i = 0;
  // this piece of codes gets 1D wellMesh of x doubling the nodes of
  // adjacent elements (for comparison with analytical solution)
  for (il::int_t e = 0; e < nelts; ++e) {
    thex[i] = mesh.coordinates(mesh.connectivity(e, 0), 0);
    thex[i + 1] = mesh.coordinates(mesh.connectivity(e, 1), 0);
    i = i + 2;
  }

  wsol =
      griffithcrack(thex, 1., myelas.Ep(), 1.);  // call to analytical solution

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

  std::cout << " end of Simple Griffith crack example P1 \n";

  return il::norm(rel_err, il::Norm::L2);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleLinearElement_AddMesh(int nelts) {
  int p = 1;
  double h = 2. / (nelts);  //  element size

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};

  // create a basic 1D wellMesh ....
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  hfp2d::Mesh mesh(xy, myconn, p);

  il::int_t ndof = mesh.numberDDDofs();

  hfp2d::ElasticProperties myelas(1, 0.);

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd,
      0.);  // passing p could be avoided here

  // add a nelts anpther element to wellMesh.
  il::int_t t1 = mesh.tipElts(0), t2 = mesh.tipElts(1);
  il::int_t n1 = mesh.tipNodes(0), n2 = mesh.tipNodes(1);

  mesh.addNTipElts(t1, n1, nelts / 2, 0.);
  mesh.addNTipElts(t2, n2, nelts / 2, 0.);

  basic_assembly_add_elts(mesh, nelts, myelas,
                          hfp2d::normal_shear_stress_kernel_dp1_dd, 0., il::io,
                          K);

  // check
  std::cout << "Nelts" << mesh.numberOfElts() << "\n";

  // solve a constant pressurized crack problem...
  il::Array<double> f{mesh.numberDDDofs(), 1.};
  // just opening dds - set shear loads to zero
  for (il::int_t i = 0; i < mesh.numberDDDofs() / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  status.ok();

  ////// Analytical solution at nodes
  il::Array<double> thex{mesh.numberDDDofs() / 2, 0},
      wsol{mesh.numberDDDofs() / 2, 0};

  il::int_t i = 0;

  // this piece of codes gets 1D wellMesh of x doubling the nodes of
  // adjacent elements (for comparison with analytical solution)
  for (il::int_t e = 0; e < mesh.numberOfElts(); ++e) {
    thex[i] = mesh.coordinates(mesh.connectivity(e, 0), 0);
    thex[i + 1] = mesh.coordinates(mesh.connectivity(e, 1), 0);
    i = i + 2;
  }

  wsol = griffithcrack(thex, 2., 1., 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  il::Array<double> rel_err{mesh.numberDDDofs() / 2 - 2};

  double rel_err_norm;
  // do not compute the relative error in x=\pm 1 as it will be infinite
  for (int j = 1; j < mesh.numberDDDofs() / 2 - 1; ++j) {
    rel_err[j - 1] = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    std::cout << "x : " << thex[j] << "..w anal:" << wsol[j]
              << " w num: " << dd[j * 2 + 1] << " rel error: " << rel_err[j - 1]
              << "\n";
  }

  std::cout << " end of Simple Griffith crack example P1 "
            << rel_err[mesh.numberDDDofs() / 2 - 4] << " \n";
  // Note that due to the weird ordering... one tip is computed in the rel_err

  return rel_err[mesh.numberDDDofs() / 2 - 4];
  il::norm(rel_err, il::Norm::L2);
}

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleS3D_P0(int nelts) {
  int p = 0;                // piece wise constant element
  double h = 2. / (nelts);  //  element size

  double Ep = 1.;  // Plane strain Young's modulus

  // create a basic 1D wellMesh ....
  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // create wellMesh object
  hfp2d::Mesh mesh(xy, myconn, p);

  hfp2d::ElasticProperties myelas(1, 0.3);

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
      1000.);  // large pseudo-heigth to reproduce plane-strain kernel

  AddTipCorrectionP0(mesh, myelas, 0, K);
  AddTipCorrectionP0(mesh, myelas, nelts - 1, K);

  il::int_t ndof = mesh.numberDDDofs();
  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, 1.};  // be careful of sign
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  //  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  for (int e = 0; e < nelts; ++e) {
    hfp2d::SegmentData sege = mesh.getElementData(e);
    thex[e] = sege.Xmid(0);
  }

  wsol =
      griffithcrack(thex, 1., myelas.Ep(), 1.);  // call to analytical solution

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

  std::cout << " end of Simple Griffith crack example P0 \n";
  status.ok();

  return il::norm(rel_err, il::Norm::L2);
}

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleS3D_P0_AddMesh(int nelts) {
  int p = 0;                // piece wise constant element
  double h = 2. / (nelts);  //  element size

  double Ep = 1.;  // Plane strain Young's modulus

  // create a basic 1D wellMesh ....
  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // create wellMesh object
  hfp2d::Mesh mesh(xy, myconn, p);

  hfp2d::ElasticProperties myelas(1, 0.3);

  il::Array2D<double> K = hfp2d::basic_assembly(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd,
      1000.);  // large pseudo-heigth to reproduce plane-strain kernel

  // add a nelts anpther element to wellMesh.
  il::int_t t1 = mesh.tipElts(0), t2 = mesh.tipElts(1);
  il::int_t n1 = mesh.tipNodes(0), n2 = mesh.tipNodes(1);

  mesh.addNTipElts(t1, n1, nelts / 2, 0.);
  mesh.addNTipElts(t2, n2, nelts / 2, 0.);

  basic_assembly_add_elts(mesh, nelts, myelas,
                          hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, 1000.,
                          il::io, K);

  AddTipCorrectionP0(mesh, myelas, mesh.tipElts(0), K);
  AddTipCorrectionP0(mesh, myelas, mesh.tipElts(1), K);

  il::int_t ndof = mesh.numberDDDofs();
  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, 1.};  // be careful of sign
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  //  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  for (int e = 0; e < mesh.numberOfElts(); ++e) {
    hfp2d::SegmentData sege = mesh.getElementData(e);
    thex[e] = sege.Xmid(0);
  }

  wsol =
      griffithcrack(thex, 2., myelas.Ep(), 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  il::Array<double> rel_err{ndof / 2};

  double rel_err_norm;
  // do not compute the relative error in x=\pm 1 as it will be infinite
  for (int j = 0; j < ndof / 2; ++j) {
    rel_err[j] = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    std::cout << "x : " << thex[j] << "..w anal:" << wsol[j]
              << " w num: " << dd[j * 2 + 1] << " rel error: " << rel_err[j]
              << "\n";
  }

  std::cout << " end of Simple Griffith crack example P0 "
            << il::norm(rel_err, il::Norm::L2) << "\n";
  status.ok();

  return il::norm(rel_err, il::Norm::L2);
}

////////////////////////////////////////////////////////////////////////////////
double SimpleGriffithExampleS3D_P0_byNodes(int nelts) {
  int p = 0;                // piece wise constant element
  double h = 2. / (nelts);  //  element size

  double Ep = 1.;  // Plane strain Young's modulus

  // create a basic 1D wellMesh ....
  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // create wellMesh object
  hfp2d::Mesh mesh(xy, myconn, p);

  hfp2d::ElasticProperties myelas(1, 0.3);

  il::Array2D<double> K = hfp2d::basic_assembly_nodal(
      mesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd_nodal,
      1000.);  // large pseudo-heigth to reproduce plane-strain kernel

  AddTipCorrectionP0(mesh, myelas, 0, K);
  AddTipCorrectionP0(mesh, myelas, nelts - 1, K);

  il::int_t ndof = mesh.numberDDDofs();
  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, 1.};  // be careful of sign
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
  //  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  for (int e = 0; e < nelts; ++e) {
    hfp2d::SegmentData sege = mesh.getElementData(e);
    thex[e] = sege.Xmid(0);
  }

  wsol =
      griffithcrack(thex, 1., myelas.Ep(), 1.);  // call to analytical solution

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

  std::cout << " end of Simple Griffith crack example P0 \n";
  status.ok();

  return il::norm(rel_err, il::Norm::L2);
}

////////////////////////////////////////////////////////////////////////////////
// Simple solver for a json mesh file and constant remote loading for stress
// - need to project on element to get tn ts (normal and shear traction) -

double SimpleCircleCrackExample_P0_byNodes(hfp2d::Mesh &MyMesh) {
  double Ep = 1.;   // Plane strain Young's modulus
  double nu = 0.0;  // Poisson's ratio

  hfp2d::ElasticProperties myelas(Ep, nu);

  //  il::Array2D<double> K = hfp2d::basic_assembly_nodal(
  //      MyMesh, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd_nodal,
  //      1000.);

  //    il::Array2D<double> K = hfp2d::basic_assembly_nodal(
  //        MyMesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal, 0.);

  il::Array2D<double> K = hfp2d::basic_assembly(
      MyMesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

  il::int_t ndof = MyMesh.numberDDDofs();
  il::int_t nelts = MyMesh.numberOfElts();

  il::Array2D<double> insitu_stress_distribution(nelts, 2);
  double far_field_vert_stress = 1.0;
  double far_field_horiz_stress = 0.0;

  for (il::int_t elmt_k = 0; elmt_k < nelts; ++elmt_k) {
    hfp2d::SegmentData mysege = MyMesh.getElementData(elmt_k);
    for (il::int_t coll_k = 0; coll_k < MyMesh.interpolationOrder() + 1;
         ++coll_k) {
      insitu_stress_distribution(elmt_k, 0) =
          -1 * (far_field_vert_stress * sin(mysege.theta())) +
          (far_field_horiz_stress * cos(mysege.theta()));
      insitu_stress_distribution(elmt_k, 1) =
          -1 * (far_field_vert_stress * cos(mysege.theta())) +
          (far_field_horiz_stress * sin(mysege.theta()));
    }
  }

  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, 1.};
  // just opening dds - set shear loads to zero
  //  for (il::int_t i = 0; i < ndof / 2; ++i) {
  //    f[2 * i] = insitu_stress_distribution(i, 0);
  //    f[2 * i + 1] = insitu_stress_distribution(i, 1);
  //  }

  il::Status status;
  // use a direct solver
  il::Array<double> dd = il::linearSolve(K, f, il::io, status);

  double A = 10.;

  return A;
}
}