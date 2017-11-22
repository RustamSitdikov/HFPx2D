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
#include <src/Devt/FiniteVolumeRoutines.h>
#include <src/FractureFluidFlow/ReynoldsP1.h>

namespace hfp2d {

Solution reynoldsP1(
    Mesh &theMesh, il::Array2D<double> &elast_submatrix,
    il::Array2D<double> &fetc_dds, il::Array2D<double> &fetc_dd,
    il::Array2D<double> &fetc_press, il::Array2D<double> &Fetc_active_dofs,
    Solution &SolutionAtTn, SimulationParameters &SimulationParameters,
    FluidProperties &FluidProperties, SolidEvolution &SolidEvolution,
    FractureEvolution &FractureEvolution, Sources &Source,
    il::Array<int> &dof_active_elmnts, il::Status &status, il::Norm &norm) {
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

  // Current pore-pressure at collocation point
  il::Array<double> press_coll{2 * theMesh.numberOfElts(), 0};
  auto p_coll = il::dot(fetc_press, SolutionAtTn.pressure());
  for (il::int_t i = 0, k = 1; i < press_coll.size(); ++i, k = k + 2) {
    press_coll[i] = p_coll[k];
  }

  // Initialization of Finite Volume matrices
  il::Array2D<double> Vd;
  il::Array2D<double> Vp;
  il::Array2D<double> L{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0.};

  // Initialization of the while loop
  il::int_t k = 0;
  double err_shearDD = SolutionAtTn.errShear();
  double err_openingDD = SolutionAtTn.errOpening();
  double err_press = SolutionAtTn.errPressure();
  il::Array<double> BigX{BigA.size(0), 0};
  il::Array<double> incrm_shearDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_shearDD_old{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> diff_incrm_shearDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openinigDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openinigDD_old{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> diff_incrm_openingDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_press{theMesh.numberOfNodes(), 0.};
  il::Array<double> incrm_press_old{theMesh.numberOfNodes(), 0.};
  il::Array<double> diff_incrm_press{theMesh.numberOfNodes(), 0.};
  il::Array<double> shearDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> shearDD_coll_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> openingDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> openingDD_coll_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> fric_coeff_k{2 * theMesh.numberOfElts(), 0.};
  il::Array2D<double> Nf{dof_active_elmnts.size(), dof_active_elmnts.size(),
                         0.};

  while (k < SimulationParameters.ehl_max_its &&
         (err_shearDD > SimulationParameters.ehl_tolerance ||
          err_openingDD > SimulationParameters.ehl_tolerance ||
          err_press > SimulationParameters.ehl_tolerance)) {
    ++k;

    // Current nodal slip
    for (il::int_t i = 0; i < incrm_shearDD.size(); ++i) {
      shearDD_k[i] = SolutionAtTn.shearDD(i) + incrm_shearDD[i];
    }
    // Current slip at collocation points
    shearDD_coll_k = il::dot(fetc_dd, shearDD_k);

    // Current opening
    for (il::int_t i = 0; i < incrm_openinigDD.size(); ++i) {
      openingDD_k[i] = SolutionAtTn.openingDD(i) + incrm_openinigDD[i];
    }
    // Current opening at collocation points
    openingDD_coll_k = il::dot(fetc_dd, openingDD_k);

    // Current friction coefficient at collocation points
    fric_coeff_k =
        SolidEvolution.linearFricWeakLaw(shearDD_coll_k, SolidEvolution);
    for (il::int_t j = 0; j < Nf.size(0); j = j + 2) {
      Nf(j, j) = fric_coeff_k[dof_active_elmnts[j] / 2];
    }

    // Finite difference matrix
    L = hfp2d::buildLMatrix(theMesh, shearDD_k, openingDD_k, FluidProperties,
                            FractureEvolution, SolutionAtTn.timestep());

    // Dilatancy matrix
    Vd = hfp2d::buildVdMatrix(theMesh, FractureEvolution, FluidProperties,
                              shearDD_k);

    // Pressure matrix
    Vp = buildVpMatrix(theMesh, FractureEvolution, FluidProperties, shearDD_k);

    /// Assembling the system ///

    // Matrix of coefficient

    auto Npk = il::dot(Nf, Fetc_active_dofs);
    for (il::int_t k2 = 0; k2 < dof_active_elmnts.size(); ++k2) {
      for (il::int_t i = 0; i < theMesh.numberOfNodes(); ++i) {
        BigA(k2, dof_active_elmnts.size() + i) = Npk(k2, i);
      }
    }

    for (il::int_t q = 0; q < Vd.size(0); ++q) {
      for (il::int_t i = 0; i < dof_active_elmnts.size(); ++i) {
        BigA(dof_active_elmnts.size() + q, i) = Vd(q, dof_active_elmnts[i]);
      }
    }

    for (il::int_t m = 0; m < Vp.size(0); ++m) {
      for (il::int_t i = 0; i < Vp.size(1); ++i) {
        BigA(dof_active_elmnts.size() + m, dof_active_elmnts.size() + i) =
            Vp(m, i) - L(m, i);
      }
    }

    // Right hand-side vector

    for (il::int_t n = 0; n < dof_active_elmnts.size(); n = n + 2) {
      BigB[n] = -((fric_coeff_k[dof_active_elmnts[n] / 2] *
                   (SolutionAtTn.sigmaN(dof_active_elmnts[n] / 2) -
                    press_coll[dof_active_elmnts[n] / 2])) -
                  SolutionAtTn.tau(dof_active_elmnts[n] / 2));
    }

    auto Lp = il::dot(L, SolutionAtTn.pressure());
    for (il::int_t i1 = 0; i1 < L.size(0); ++i1) {
      BigB[dof_active_elmnts.size() + i1] =
          Lp[i1];  // No source term!! add here source term for injection rate
                   // for instance
    }

    /// Boundary conditions ///

    for (il::int_t k1 = 0; k1 < BigA.size(1); ++k1) {
      BigA(dof_active_elmnts.size() + Source.getSourcePoint(), k1) = 0;
    }

    for (il::int_t l1 = 0; l1 < BigA.size(0); ++l1) {
      BigA(l1, dof_active_elmnts.size() + Source.getSourcePoint()) = 0;
    }

    BigA(dof_active_elmnts.size() + Source.getSourcePoint(),
         dof_active_elmnts.size() + Source.getSourcePoint()) = 1;

    BigB[dof_active_elmnts.size() + Source.getSourcePoint()] = 0;

    /// Solve the system ///

    BigX = il::linearSolve(BigA, BigB, il::io, status);
    status.abortOnError();

    ///  Under relaxation technique & updating ///

    for (il::int_t q = 0; q < dof_active_elmnts.size(); q = q + 2) {
      incrm_shearDD[dof_active_elmnts[q] / 2] =
          (1 - SimulationParameters.ehl_relaxation) *
              incrm_shearDD_old[dof_active_elmnts[q] / 2] +
          SimulationParameters.ehl_relaxation * BigX[q];
    }

    for (il::int_t q = 1; q <= dof_active_elmnts.size(); q = q + 2) {
      incrm_openinigDD[dof_active_elmnts[q] / 2] =
          (1 - SimulationParameters.ehl_relaxation) *
              incrm_openinigDD_old[dof_active_elmnts[q] / 2] +
          SimulationParameters.ehl_relaxation * BigX[q];
    }

    if (dof_active_elmnts.size() == 0) {
      for (il::int_t i = 0; i < incrm_press.size(); ++i) {
        incrm_press[i] = BigX[dof_active_elmnts.size() + i];
      }
    } else {
      for (il::int_t q = 0; q < incrm_press.size(); ++q) {
        incrm_press[q] =
            (1 - SimulationParameters.ehl_relaxation) * incrm_press_old[q] +
            SimulationParameters.ehl_relaxation *
                BigX[dof_active_elmnts.size() + q];
      }
    }

    // Error on increments
    for (il::int_t i2 = 0; i2 < diff_incrm_shearDD.size(); ++i2) {
      diff_incrm_shearDD[i2] = incrm_shearDD[i2] - incrm_shearDD_old[i2];
    }

    for (il::int_t i2 = 0; i2 < diff_incrm_openingDD.size(); ++i2) {
      diff_incrm_openingDD[i2] =
          incrm_openinigDD[i2] - incrm_openinigDD_old[i2];
    }

    for (il::int_t i2 = 0; i2 < diff_incrm_press.size(); ++i2) {
      diff_incrm_press[i2] = incrm_press[i2] - incrm_press_old[i2];
    }

    if (dof_active_elmnts.size() == 0) {
      err_shearDD = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_shearDD =
          (il::norm(diff_incrm_shearDD, norm) / il::norm(incrm_shearDD, norm));
    }

    if (dof_active_elmnts.size() == 0 ||
        il::norm(incrm_openinigDD, norm) == 0) {
      err_openingDD = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_openingDD = (il::norm(diff_incrm_openingDD, norm) /
                       il::norm(incrm_openinigDD, norm));
    }

    if (dof_active_elmnts.size() == 0) {
      err_press = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_press =
          (il::norm(diff_incrm_press, norm) / il::norm(incrm_press, norm));
    }

    std::cout << "error on shearDD : " << err_shearDD
              << "  error on openingDD: " << err_openingDD
              << "  error on pressure: " << err_press << "\n";

    // Update
    incrm_openinigDD_old = incrm_openinigDD;
    incrm_shearDD_old = incrm_shearDD;
    incrm_press_old = incrm_press;
  }


  return SolutionAtTn;
};
}