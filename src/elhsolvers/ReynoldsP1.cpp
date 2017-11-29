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
#include <src/devt/FiniteVolumeRoutines.h>
#include <src/elhsolvers/ReynoldsP1.h>

namespace hfp2d {

Solution reynoldsP1(
    Mesh &theMesh, il::Array2D<double> &elast_matrix,
    il::Array2D<double> &fetc_dds, il::Array2D<double> &fetc_dd,
    il::Array2D<double> &fetc_press, Solution &SolutionAtTn,
    il::Array<double> &tau_old, il::Array<double> &sigmaN_old,
    il::Array<double> &shearDD_old, il::Array<double> &openingDD_old,
    il::Array<double> &press_old, il::Array<double> &incrm_shearDD,
    il::Array<double> &incrm_openingDD,
    SimulationParameters &SimulationParameters,
    FluidProperties &FluidProperties, SolidEvolution &SolidEvolution,
    FractureEvolution &FractureEvolution, Sources &Source,
    il::Array<int> &dof_active_elmnts, il::Status &status, il::Norm &norm,
    bool damping_term, double damping_coeff, double previous_time) {
  //// IMPLICIT SOLUTION OF THE COUPLED PROBLEM ////
  // Initialization of the system BigA*BigX = BigB
  il::Array2D<double> BigA{dof_active_elmnts.size() + theMesh.numberOfNodes(),
                           dof_active_elmnts.size() + theMesh.numberOfNodes(),
                           0.};
  il::Array<double> BigB{dof_active_elmnts.size() + theMesh.numberOfNodes(),
                         0.};

  // Select the elasticity matrix for shear and opening DDs of slipping DOFs
  il::Array2D<double> elast_submatrix{dof_active_elmnts.size(),
                                      dof_active_elmnts.size(), 0};
  for (il::int_t m = 0; m < elast_submatrix.size(0); ++m) {
    for (il::int_t i = 0; i < elast_submatrix.size(1); ++i) {
      elast_submatrix(m, i) =
          elast_matrix(dof_active_elmnts[m], dof_active_elmnts[i]);
    }
  }

  if (damping_term == true) {
    for (il::int_t i = 0; i < elast_submatrix.size(0); ++i) {
      elast_submatrix(i, i) =
          elast_submatrix(i, i) - (damping_coeff / SolutionAtTn.timestep());
    }
  }

  // Assembling elasticity part of BigA
  for (int i = 0; i < elast_submatrix.size(0); ++i) {
    for (int j = 0; j < elast_submatrix.size(1); ++j) {
      BigA(i, j) = elast_submatrix(i, j);
    }
  }

  // Current pore-pressure at collocation point
  il::Array<double> press_coll{2 * theMesh.numberOfElts(), 0};
  auto p_coll = il::dot(fetc_press, press_old);
  for (il::int_t i = 0, k = 1; i < press_coll.size(); ++i, k = k + 2) {
    press_coll[i] = p_coll[k];
  }

  // Select just the "active" part of the matrix to switch from nodal
  // points to collocation points
  il::Array2D<double> Fetc_active_dofs{dof_active_elmnts.size(),
                                       SolutionAtTn.pressure().size(), 0.};
  for (il::int_t l2 = 0; l2 < Fetc_active_dofs.size(0); ++l2) {
    for (il::int_t i = 0; i < Fetc_active_dofs.size(1); ++i) {
      Fetc_active_dofs(l2, i) = fetc_press(dof_active_elmnts[l2], i);
    }
  }

  // Initialization of Finite Volume matrices
  il::Array2D<double> Vd;
  il::Array2D<double> Vp;
  il::Array2D<double> L{theMesh.numberOfNodes(), theMesh.numberOfNodes(), 0.};

  // Initialization of the while loop
  il::int_t k = 0;
  double err_shearDD = 2.;
  double err_openingDD = 2.;
  double err_press = 2.;
  il::Array<double> BigX{BigA.size(0), 0};
  il::Array<double> incrm_shearDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_shearDD_k_old{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> diff_incrm_shearDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openingDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openingDD_k_old{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> diff_incrm_openingDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_press_k{theMesh.numberOfNodes(), 0.};
  il::Array<double> incrm_press_k_old{theMesh.numberOfNodes(), 0.};
  il::Array<double> diff_incrm_press_k{theMesh.numberOfNodes(), 0.};
  il::Array<double> shearDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> shearDD_coll_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> openingDD_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> openingDD_coll_k{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> fric_coeff_k{2 * theMesh.numberOfElts(), 0.};
  il::Array2D<double> Nf{dof_active_elmnts.size(), dof_active_elmnts.size(),
                         0.};

  while ((k < SimulationParameters.ehl_max_its) &&
         (err_shearDD > SimulationParameters.ehl_tolerance ||
          err_openingDD > SimulationParameters.ehl_tolerance ||
          err_press > SimulationParameters.ehl_tolerance)) {
    ++k;

    // Current nodal slip
    for (il::int_t i = 0; i < incrm_shearDD_k.size(); ++i) {
      shearDD_k[i] = shearDD_old[i] + incrm_shearDD_k[i];
    }
    // Current slip at collocation points
    shearDD_coll_k = il::dot(fetc_dd, shearDD_k);

    // Current opening
    for (il::int_t i = 0; i < incrm_openingDD_k.size(); ++i) {
      openingDD_k[i] = openingDD_old[i] + incrm_openingDD_k[i];
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

    auto Lp = il::dot(L, press_old);
    for (il::int_t i1 = 0; i1 < L.size(0); ++i1) {
      BigB[dof_active_elmnts.size() + i1] =
          Lp[i1];  // No source term!! add here source term for injection rate
                   // for instance
    }

    /// Boundary conditions ///

    for (il::int_t k1 = 0; k1 < BigA.size(1); ++k1) {
      BigA(dof_active_elmnts.size() + Source.getSourcePoint(), k1) = 0;
    }

//    for (il::int_t l1 = 0; l1 < BigA.size(0); ++l1) {
//      BigA(l1, dof_active_elmnts.size() + Source.getSourcePoint()) = 0;
//    }

    BigA(dof_active_elmnts.size() + Source.getSourcePoint(),
         dof_active_elmnts.size() + Source.getSourcePoint()) = 1;

    BigB[dof_active_elmnts.size() + Source.getSourcePoint()] = 0;

    /// Solve the system ///

    BigX = il::linearSolve(BigA, BigB, il::io, status);
    status.abortOnError();

    ///  Under relaxation technique & updating ///

    for (il::int_t q = 0; q < dof_active_elmnts.size(); q = q + 2) {
      incrm_shearDD_k[dof_active_elmnts[q] / 2] =
          (1 - SimulationParameters.ehl_relaxation) *
              incrm_shearDD_k_old[dof_active_elmnts[q] / 2] +
          SimulationParameters.ehl_relaxation * BigX[q];
    }

    for (il::int_t q = 1; q <= dof_active_elmnts.size(); q = q + 2) {
      incrm_openingDD_k[dof_active_elmnts[q] / 2] =
          (1 - SimulationParameters.ehl_relaxation) *
              incrm_openingDD_k_old[dof_active_elmnts[q] / 2] +
          SimulationParameters.ehl_relaxation * BigX[q];
    }

    if (dof_active_elmnts.size() == 0) {
      for (il::int_t i = 0; i < incrm_press_k.size(); ++i) {
        incrm_press_k[i] = BigX[dof_active_elmnts.size() + i];
      }
    } else {
      for (il::int_t q = 0; q < incrm_press_k.size(); ++q) {
        incrm_press_k[q] =
            (1 - SimulationParameters.ehl_relaxation) * incrm_press_k_old[q] +
            SimulationParameters.ehl_relaxation *
                BigX[dof_active_elmnts.size() + q];
      }
    }

    // Error on increments
    for (il::int_t i2 = 0; i2 < diff_incrm_shearDD_k.size(); ++i2) {
      diff_incrm_shearDD_k[i2] = incrm_shearDD_k[i2] - incrm_shearDD_k_old[i2];
    }

    for (il::int_t i2 = 0; i2 < diff_incrm_openingDD_k.size(); ++i2) {
      diff_incrm_openingDD_k[i2] =
          incrm_openingDD_k[i2] - incrm_openingDD_k_old[i2];
    }

    for (il::int_t i2 = 0; i2 < diff_incrm_press_k.size(); ++i2) {
      diff_incrm_press_k[i2] = incrm_press_k[i2] - incrm_press_k_old[i2];
    }

    if (dof_active_elmnts.size() == 0) {
      err_shearDD = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_shearDD = (il::norm(diff_incrm_shearDD_k, norm) /
                     il::norm(incrm_shearDD_k, norm));
    }

    if (dof_active_elmnts.size() == 0 ||
        il::norm(incrm_openingDD_k, norm) == 0) {
      err_openingDD = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_openingDD = (il::norm(diff_incrm_openingDD_k, norm) /
                       il::norm(incrm_openingDD_k, norm));
    }

    if (dof_active_elmnts.size() == 0) {
      err_press = SimulationParameters.ehl_tolerance / 2;
    } else {
      err_press =
          (il::norm(diff_incrm_press_k, norm) / il::norm(incrm_press_k, norm));
    }

//    std::cout << "error on shearDD : " << err_shearDD
//              << "  error on openingDD: " << err_openingDD
//              << "  error on pressure: " << err_press << "\n";

    // Update
    incrm_openingDD_k_old = incrm_openingDD_k;
    incrm_shearDD_k_old = incrm_shearDD_k;
    incrm_press_k_old = incrm_press_k;
  }

  for (il::int_t l = 0; l < incrm_shearDD.size(); ++l) {
    incrm_shearDD[l] = (incrm_shearDD[l] + incrm_shearDD_k[l]);
  }

  for (il::int_t l = 0; l < incrm_openingDD.size(); ++l) {
    incrm_openingDD[l] = (incrm_openingDD[l] + incrm_openingDD_k[l]);
  }

  il::Array2D<double> elast_matrix_shear{2 * theMesh.numberOfElts(),
                                         2 * theMesh.numberOfElts(), 0.};
  for (il::int_t m1 = 0, n1 = 0; m1 < elast_matrix_shear.size(0);
       ++m1, n1 = n1 + 2) {
    for (il::int_t i = 0, j = 0; i < elast_matrix_shear.size(1);
         ++i, j = j + 2) {
      elast_matrix_shear(m1, i) = -1 * elast_matrix(n1, j);
    }
  }

  if (damping_term == true) {
    for (il::int_t i = 0; i < elast_matrix_shear.size(0); ++i) {
      elast_matrix_shear(i, i) =
          elast_matrix_shear(i, i) - (damping_coeff / SolutionAtTn.timestep());
    }
  }

  il::Array2D<double> elast_matrix_sigmaN{2 * theMesh.numberOfElts(),
                                          2 * theMesh.numberOfElts(), 0.};
  for (il::int_t m1 = 0, n1 = 1; m1 < elast_matrix_sigmaN.size(0);
       ++m1, n1 = n1 + 2) {
    for (il::int_t i = 0, j = 1; i < elast_matrix_sigmaN.size(1);
         ++i, j = j + 2) {
      elast_matrix_sigmaN(m1, i) = -1 * elast_matrix(n1, j);
    }
  }

  il::Array<double> incrm_shear_stress{2 * theMesh.numberOfElts(), 0.};
  incrm_shear_stress = il::dot(elast_matrix_shear, incrm_shearDD);

  il::Array<double> incrm_normal_stress{2 * theMesh.numberOfElts(), 0.};
  incrm_normal_stress = il::dot(elast_matrix_sigmaN, incrm_openingDD);

  /// Calculate new stress state, new pore pressure new cumulative DDs and new
  /// friction

  // New shear stress
  il::Array<double> tau_new{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t j2 = 0; j2 < tau_new.size(); ++j2) {
    tau_new[j2] = tau_old[j2] + incrm_shear_stress[j2];
  }

  // Force M-C criterion
  for (il::int_t j3 = 0; j3 < dof_active_elmnts.size(); j3 = j3 + 2) {
    tau_new[dof_active_elmnts[j3] / 2] =
        fric_coeff_k[dof_active_elmnts[j3] / 2] *
        (SolutionAtTn.sigmaN(dof_active_elmnts[j3] / 2) -
         press_coll[dof_active_elmnts[j3] / 2]);
  }

  // New normal stress
  il::Array<double> sigmaN_new{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t j2 = 0; j2 < sigmaN_new.size(); ++j2) {
    sigmaN_new[j2] = sigmaN_old[j2] + incrm_normal_stress[j2];
  }

  // New pore pressure profile
  il::Array<double> pore_press_new{theMesh.numberOfNodes(), 0.};
  for (il::int_t m2 = 0; m2 < pore_press_new.size(); ++m2) {
    pore_press_new[m2] = press_old[m2] + incrm_press_k[m2];
  }

  // New slip vector (shear DD)
  il::Array<double> shearDD_new{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t n2 = 0; n2 < shearDD_new.size(); ++n2) {
    shearDD_new[n2] = shearDD_old[n2] + incrm_shearDD[n2];
  }

  // New opening vector
  il::Array<double> openingDD_new{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t i3 = 0; i3 < openingDD_new.size(); ++i3) {
    openingDD_new[i3] = openingDD_old[i3] + incrm_openingDD[i3];
  }

  // New friction coefficient
  SolidEvolution.setFrictionCoefficient(fric_coeff_k);

  return hfp2d::Solution(theMesh, previous_time,
                         SolutionAtTn.timestep(), openingDD_new, shearDD_new,
                         pore_press_new, sigmaN_new, tau_new,
                         SolutionAtTn.activeElts(), SolutionAtTn.frontIts(),
                         k, SolutionAtTn.errFront(),
                         err_openingDD, err_shearDD, err_press);
};
}