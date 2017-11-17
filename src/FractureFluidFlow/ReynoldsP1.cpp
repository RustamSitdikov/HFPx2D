//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 17.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard libtrary
#include <algorithm>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include <src/FractureFluidFlow/ReynoldsP1.h>
#include <src/Devt/FiniteVolumeRoutines.h>

namespace hfp2d {

Solution reynoldsP1(Mesh &theMesh, il::Array2D<double> &elast_submatrix,
                    il::Array2D<double> &fetc_dd,
                    il::Array2D<double> &fetc_press, Solution &SolutionAtTn,
                    SimulationParameters &SimulationParameters,
                    SolidEvolution &SolidEvolution,
                    FractureEvolution &FractureEvolution, Sources &Source,
                    il::Array<int> &dof_active_elmnts, il::Status &status,
                    il::Norm &norm) {
  //// IMPLICIT SOLUTION OF THE COUPLED PROBLEM ////
  // Initialization of the system BigA*BigX = BigB
  il::Array2D<double> BigA{dof_active_elmnts.size() + theMesh.numberOfNodes(),
                           dof_active_elmnts.size() + theMesh.numberOfNodes(),
                           0.};
  il::Array<double> BigB{dof_active_elmnts.size() + theMesh.numberOfNodes(),
                         0.};

  // Assembling elasticity part of BigA
  for (int i = 0; i < elast_submatrix.size(0); ++i) {
    for (int j = 0; j < elast_submatrix.size(1); ++j) {
      BigA(i, j) = elast_submatrix(i, j);
    }
  }

  // Initialization of Finite Volume matrices
  il::Array2D<double> Vd;
  il::Array2D<double> Vp;
  il::Array2D<double> L{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0.};

  // Initialization of the while loop
  int j = 0;

  while (j < SimulationParameters.EHL_max_its &&
         (SolutionAtTn.err_opening() > SimulationParameters.EHL_tolerance ||
          SolutionAtTn.err_shear() > SimulationParameters.EHL_tolerance ||
          SolutionAtTn.err_pressure() > SimulationParameters.EHL_tolerance)) {
      ++j;


  }
};
}