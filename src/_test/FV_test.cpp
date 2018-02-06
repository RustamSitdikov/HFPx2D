//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 05.02.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>

#include <gtest/gtest.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include <src/core/Mesh.h>
#include <src/devt/FiniteVolumeRoutines.h>

//// TESTS FOR FINITE VOLUME MATRICES ///
// Control volumes centred on mesh nodes
// Spatial integration via Cavalieri-Simpson rule

/// TEST 1
TEST(FV, test_finite_difference_matrix) {
  //  a very simple mesh with 4 elements  (0,1,2,3)
  // ensure the assembling of finite difference matrix

  // create the Mesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  //  // create a basic 1D wellMesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // Interpolation order
  il::int_t p = 1;

  hfp2d::Mesh mesh(xy, myconn, p);

  // Null slip vector -> easier
  il::Array<double> slip{2 * nelts, 0.};

  // Null opening vector -> easier
  il::Array<double> opening{2 * nelts, 0.};

  // Fluid properties
  double density = 1.;
  double viscosity = 0.003;
  double compressibility = 1.;
  hfp2d::Fluid fluidProperties(density, compressibility, viscosity);

  // Time step
  double dt = 0.1;

  // Fracture properties
  il::Array<double> init_permeab{2 * nelts, 0.36};
  il::Array<double> incr_permeab{2 * nelts, 0.};
  il::Array<double> init_hydr_width{2 * nelts, 0.6};
  il::Array<double> incr_hydr_width{2 * nelts, 0.};
  il::Array<double> res_slip{2 * nelts, 0.4};

  hfp2d::FractureEvolution fractureEvolution(
      init_permeab, incr_permeab, res_slip, init_hydr_width, incr_hydr_width);

  auto L = hfp2d::buildLMatrix(mesh, slip, opening, fluidProperties,
                               fractureEvolution, dt);

  il::Array2D<double> Lma{nelts + 1, nelts + 1, 0.};
  Lma(0, 0) = -0.6;
  Lma(0, 1) = 0.6;
  Lma(1, 0) = 0.6;
  Lma(1, 1) = -1.2;
  Lma(1, 2) = 0.6;
  Lma(2, 1) = 0.6;
  Lma(2, 2) = -1.2;
  Lma(2, 3) = 0.6;
  Lma(3, 2) = 0.6;
  Lma(3, 3) = -1.2;
  Lma(3, 4) = 0.6;
  Lma(4, 3) = 0.6;
  Lma(4, 4) = -0.6;

  double my_sum = 0.;

  for (il::int_t j = 0; j < L.size(1); j++) {
    for (il::int_t i = 0; i < L.size(0); i++) {
      my_sum += abs(L(i, j) - Lma(i, j));
    }
  }

  std::cout << my_sum << "\n";

  ASSERT_NEAR(my_sum, 0., 1.e-5);
}

/// TEST 2
TEST(FV, test_mass_matrix) {
  //  a very simple mesh with 4 elements  (0,1,2,3)
  // ensure the assembling of finite difference matrix

  // create the Mesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  //  // create a basic 1D wellMesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // Interpolation order
  il::int_t p = 1;

  hfp2d::Mesh mesh(xy, myconn, p);

  // Null slip vector -> easier
  il::Array<double> slip{2 * nelts, 0.};

  // Fluid properties
  double density = 1.;
  double viscosity = 0.003;
  double compressibility = 1.;
  hfp2d::Fluid fluidProperties(density, compressibility, viscosity);

  // Fracture properties
  il::Array<double> init_permeab{2 * nelts, 0.36};
  il::Array<double> incr_permeab{2 * nelts, 0.};
  il::Array<double> init_hydr_width{2 * nelts, 0.6};
  il::Array<double> incr_hydr_width{2 * nelts, 0.1};
  il::Array<double> res_slip{2 * nelts, 0.4};

  hfp2d::FractureEvolution fractureEvolution(
      init_permeab, incr_permeab, res_slip, init_hydr_width, incr_hydr_width);

  auto Vd =
      hfp2d::buildVdMatrix(mesh, fractureEvolution, fluidProperties, slip);

  il::Array2D<double> Vdma{nelts + 1, 4 * nelts, 0.};
  Vdma(0, 0) = 0.09375;
  Vdma(0, 2) = 0.03125;
  Vdma(1, 0) =  0.03125;
  Vdma(1, 2) = 0.09375;
  Vdma(1, 4) = 0.09375;
  Vdma(1, 6) =  0.03125;
  Vdma(2, 4) =  0.03125;
  Vdma(2, 6) = 0.09375;
  Vdma(2, 8) = 0.09375;
  Vdma(2, 10) =  0.03125;
  Vdma(3, 8) =  0.03125;
  Vdma(3, 10) = 0.09375;
  Vdma(3, 12) = 0.09375;
  Vdma(3, 14) =  0.03125;
  Vdma(4, 12) =  0.03125;
  Vdma(4, 14) = 0.09375;

  double my_sum = 0.;

  for (il::int_t j = 0; j < Vd.size(1); j++) {
    for (il::int_t i = 0; i < Vd.size(0); i++) {
      my_sum += abs(Vd(i, j) - Vdma(i, j));
    }
  }

  std::cout << my_sum << "\n";

  ASSERT_NEAR(my_sum, 0., 1.e-5);
}

/// TEST 3
TEST(FV, test_compressibility_matrix) {
  //  a very simple mesh with 4 elements  (0,1,2,3)
  // ensure the assembling of finite difference matrix

  // create the Mesh.
  il::int_t nelts = 4;

  il::Array2D<double> xy{nelts + 1, 2, 0.};

  //  // create a basic 1D wellMesh .... first fracture
  for (il::int_t i = 0; i < nelts + 1; ++i) {
    xy(i, 0) = -1.0 + 1. * i;
    xy(i, 1) = 0.;
  };

  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  for (il::int_t i = 0; i < nelts; ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // Interpolation order
  il::int_t p = 1;

  hfp2d::Mesh mesh(xy, myconn, p);

  // Null slip vector -> easier
  il::Array<double> slip{2 * nelts, 0.};

  // Fluid properties
  double density = 1.;
  double viscosity = 0.003;
  double compressibility = 1.;
  hfp2d::Fluid fluidProperties(density, compressibility, viscosity);

  // Time step
  double dt = 0.1;

  // Fracture properties
  il::Array<double> init_permeab{2 * nelts, 0.36};
  il::Array<double> incr_permeab{2 * nelts, 0.};
  il::Array<double> init_hydr_width{2 * nelts, 0.6};
  il::Array<double> incr_hydr_width{2 * nelts, 0.};
  il::Array<double> res_slip{2 * nelts, 0.4};

  hfp2d::FractureEvolution fractureEvolution(
      init_permeab, incr_permeab, res_slip, init_hydr_width, incr_hydr_width);

  auto Vp =
      hfp2d::buildVpMatrix(mesh, fractureEvolution, fluidProperties, slip);

  il::Array2D<double> Vpma{nelts + 1, nelts + 1, 0.};
  Vpma(0, 0) = 0.225;
  Vpma(0, 1) = 0.075;
  Vpma(1, 0) = 0.075;
  Vpma(1, 1) = 0.45;
  Vpma(1, 2) = 0.075;
  Vpma(2, 1) = 0.075;
  Vpma(2, 2) = 0.45;
  Vpma(2, 3) = 0.075;
  Vpma(3, 2) = 0.075;
  Vpma(3, 3) = 0.45;
  Vpma(3, 4) = 0.075;
  Vpma(4, 3) = 0.075;
  Vpma(4, 4) = 0.225;

  double my_sum = 0.;

  for (il::int_t j = 0; j < Vp.size(1); j++) {
    for (il::int_t i = 0; i < Vp.size(0); i++) {
      my_sum += abs(Vp(i, j) - Vpma(i, j));
    }
  }

  std::cout << my_sum << "\n";

  ASSERT_NEAR(my_sum, 0., 1.e-5);
}
