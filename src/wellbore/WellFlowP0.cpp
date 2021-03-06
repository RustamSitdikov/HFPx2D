//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <stdio.h>
#include <algorithm>
#include <cmath>

#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/blas/blas.h>
#include <il/linear_algebra/dense/blas/dot.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include <il/linear_algebra/dense/factorization/linearSolve.h>
#include <il/linear_algebra/dense/norm.h>

#include <src/core/Mesh.h>
#include <src/wellbore/WellFlowP0.h>
#include <src/util/RootFinder.h>

//#include <src/wellbore/WellMesh.h>
//#include <src/wellbore/WellSolution.h>

namespace hfp2d {

// bracketing the root
imf::BracketS bracketingHD(imf::ImplicitFunction fun, IFParametersHD &params,
                           double lw_bound, double up_bound, const int max_iter,
                           bool mute) {
  imf::BracketS ab;

  ab.l_b = lw_bound;
  ab.u_b = up_bound;

  double fa = fun(ab.l_b, params);
  double fb = fun(ab.u_b, params);

  int count = 0;
  double wt = 2.0;
  double mid = ab.u_b;
  while (fa * fb > 0.0 && count < max_iter) {
    // take an intermediate point for lower bound
    mid = (ab.l_b + wt * mid) / (wt + 1.0);
    ab.l_b = mid;
    fa = fun(ab.l_b, params);
    ++count;
  }

  if (count >= max_iter && fa * fb > 0.0) {
    // if root isn't bracketed
    if (!mute) {
      std::cout << "The root of f(s) can't be bracketed in " << count
                << " iterations" << std::endl;
    }
  } else if (!mute) {
    std::cout << "The root of f(s) is bracketed between " << ab.l_b << " and "
              << ab.u_b << " in " << count << " iterations" << std::endl;
  }
  return ab;
}

// overload: bracketing the root (friction factor)
imf::BracketS bracketingHD(imf::ImplicitFunction fun, IFParametersHD &params,
                           const int max_iter, bool mute) {
  const double Re_L = 2300., Re_T = 4000.;

  imf::BracketS ab;
  IFParametersHD par_T;
  par_T.Re_num = Re_T;
  par_T.rough = params.rough;

  // laminar regime
  double ff_lam = 64. / params.Re_num;
  // Maximum drag reduction (MDR) approximation
  double ff_MDR_appr = ffMDRappr(params);
  // Gauckler-Manning-Strickler asymptote (rough turbulence)
  double ff_GMS = 0.143 * std::pow(params.rough, 1. / 3.);
  // Reichert asymptote (rough turbulence)
  double ff_rei = 4. / std::pow(2.28 - 4. * std::log10(params.rough / 2.), 2);
  // ~ max. of Churchill approximation
  double ff_T = ffChurchill(par_T);

  if (params.Re_num < Re_L) {
    ab.l_b = ff_lam / 2.;
    ab.u_b = ff_lam * 2.;
  } else if (params.Re_num < Re_T) {
    ab.l_b = ff_lam / 2.;
    ab.u_b = il::max(ff_GMS, ff_rei, ff_T) * 2.;
  } else {
    // ab.l_b = ff_bla / 4.;
    ab.l_b = ff_MDR_appr / 4.;
    ab.u_b = il::max(ff_GMS, ff_rei, ff_T) * 2.;
  }

  double fb = fun(ab.u_b, params);

  double fa = fun(ab.l_b, params);

  int count = 0;
  double wt = 2.0;
  double mid = ab.u_b;
  while (fa * fb > 0.0 && count < max_iter) {
    // take an intermediate point for lower bound
    mid = (ab.l_b + wt * mid) / (wt + 1.0);
    ab.l_b = mid;
    fa = fun(ab.l_b, params);
    ++count;
  }

  if (count >= max_iter && fa * fb > 0.0) {
    // if root isn't bracketed
    if (!mute) {
      std::cout << "The root of f(s) can't be bracketed in " << count
                << " iterations" << std::endl;
    }
  } else if (!mute) {
    std::cout << "The root of f(s) is bracketed between " << ab.l_b << " and "
              << ab.u_b << " in " << count << " iterations" << std::endl;
  }
  return ab;
}

//////////////////////////////////////////////////////////////////////////
//        FRICTION FACTOR MODELS
//////////////////////////////////////////////////////////////////////////

//--------------------------------------------------
// Yang & Dou model (implicit)
double imYangDou(double s, IFParametersHD &params) {
  // todo: move e, pi, etc to core
  const double e = 2.71828182845904523536;
  const int n_max = 5;

  double v_sc = std::sqrt(8. / s);
  double r = params.Re_num / (2. * v_sc);  // 1/2 HD
  double r_star = r * params.rough;
  double theta = il::pi * std::log(r_star / 1.25) / std::log(100. / 1.25);
  double alpha = (1. - std::cos(theta)) / 2.;
  double beta = 1. - (1. - 0.107) * (alpha + theta / il::pi) / 2.;

  double a_r_s = alpha * r_star;
  double a_b_r_s = beta * a_r_s;
  double c_a = a_r_s / (5. + a_r_s);
  double c_b = a_b_r_s / (5. + a_b_r_s);
  double c_l = (5. + a_r_s) / (5. + a_b_r_s);

  double rt = 1.;
  int n_fact = 1;
  for (int n = 1; n < n_max; ++n) {
    n_fact *= n;
    rt -= 1. / e * ((double)n / (double)n_fact * std::pow(67.8 / r, 2 * n));
  }

  double res =
      v_sc - ((1. - rt) * r / 4. +
              rt * (2.5 * std::log(r) - 66.69 * std::pow(r, -0.72) + 1.8 -
                    (2.5 * std::log(c_l) + (5.8 + 1.25) * std::pow(c_a, 2) +
                     2.5 * c_a - (5.8 + 1.25) * std::pow(c_b, 2) - 2.5 * c_b)));

  return res;
}

//--------------------------------------------------
// Mixed model: Laminar / Yang-Dou / Reichert
double ffYDRmixed(IFParametersHD &params) {
  const bool mute = true;

  const double epsilon = 1E-12;
  const int max_iter = 100;
  const double Re_L = 2300.;

  double ff;
  // laminar regime
  double ff_lam = 64. / params.Re_num;
  // Reichert asymptote (rough turbulence)
  double ff_rei = 4. / std::pow(2.28 - 4. * std::log10(params.rough / 2.), 2);

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else {
    // bracketing
    imf::BracketS ab = bracketingHD((imf::ImplicitFunction)imYangDou, params,
                                    // ff_lam / 2., 3200.,
                                    max_iter, mute);
    // finding the root (friction factor)
    ff = imf::brent((imf::ImplicitFunction)imYangDou, params, ab.l_b, ab.u_b,
                    epsilon, max_iter, mute);

    double v_sc = std::sqrt(8. / ff);
    double r_plus = params.Re_num / (2. * v_sc);
    // use Reichert asymptote for high Re
    if (r_plus > 100. / params.rough) {
      ff = ff_rei;
    }
    if (params.rough > 1. / 32. && ff > ff_rei) {
      ff = ff_rei;
    }
  }
  return ff;
}

//--------------------------------------------------
// Mixed model: Laminar / Yang-Dou / Gauckler-Manning-Strickler
double ffYDSmixed(IFParametersHD &params) {
  const bool mute = true;

  const double epsilon = 1E-12;
  const int max_iter = 100;
  const double Re_L = 2300.;

  double ff;
  // laminar regime
  double ff_lam = 64. / params.Re_num;
  // Gauckler-Manning-Strickler asymptote (rough turbulence)
  double ff_GMS = 4. * 0.143 * std::pow(params.rough, 1. / 3.);

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else {
    // bracketing
    imf::BracketS ab = bracketingHD((imf::ImplicitFunction)imYangDou, params,
                                    // ff_lam / 2., 3200.,
                                    max_iter, mute);
    // finding the root (friction factor)
    ff = imf::brent((imf::ImplicitFunction)imYangDou, params, ab.l_b, ab.u_b,
                    epsilon, max_iter, mute);

    double v_sc = std::sqrt(8. / ff);
    double r_plus = params.Re_num / (2. * v_sc);
    // use Gauckler-Manning-Strickler asymptote for high Re
    if (r_plus > 100. / params.rough) {
      ff = ff_GMS;
    }
    if (1. / params.rough < 32. && ff > ff_GMS) {
      ff = ff_GMS;
    }
  }
  return ff;
}

//--------------------------------------------------
// Churchill model (explicit)
double ffChurchill(IFParametersHD &params) {
  double rey = params.Re_num;  // rescale if necessary

  double rou = params.rough / 2.;

  double f0 = -std::log(std::pow(7. / rey, 0.9) + 0.27 * rou);

  double ff = 8. * std::pow(std::pow(8. / rey, 12) +
                                std::pow(std::pow(2.457 * f0, 16) +
                                             std::pow(37530. / rey, 16),
                                         -1.5),
                            1. / 12.);
  return ff;
}

//--------------------------------------------------
// Haaland model (explicit)
double ffHaaland(IFParametersHD &params) {
  double rey = params.Re_num;  // rescale if necessary

  double rou = params.rough;

  double ff =
      std::pow(1.8 * std::log10(std::pow(rou / 3.7, 1.11) + 6.9 / rey), -2);
  return 4. * ff;
}

// Maximum drag reduction (MDR) asymptote (implicit)
double imMDR(double s, IFParametersHD &params) {
  double sf = std::sqrt(s / 4.);

  // log10(Re_num / sqrt(s)) must be a misprint
  double res = 1. / sf - (19. * std::log10(params.Re_num * sf) - 32.4);

  return res;
}

//--------------------------------------------------
// Maximum drag reduction (MDR) approximation, explicit
double ffMDRappr(IFParametersHD &params) {
  const double Re_L = 1510.08;

  double ff;

  double ff_lam = 64. / params.Re_num;

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else {
    ff = 4. * 1.78 * std::pow(params.Re_num, -0.7);
  }
  return ff;
}

//--------------------------------------------------
// Mixed model: Laminar / MDR
double ffMDRmixed(IFParametersHD &params) {
  const bool mute = true;

  const double Re_L = 1379.33, Re_H = 14022.1;

  double ff;
  // laminar regime
  double ff_lam = 64. / params.Re_num;
  // Blasius approximation
  double ff_bla = 0.316 * std::pow(params.Re_num, -0.25);
  // Maximum drag reduction (MDR) approximation
  double ff_MDR_appr = ffMDRappr(params);

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else {
    imf::BracketS ab;
    if (params.Re_num < Re_H) {
      // bracketing
      ab = bracketingHD((imf::ImplicitFunction)imMDR, params, ff_lam / 2.,
                        ff_MDR_appr * 2., max_iter, mute);
    } else {
      // bracketing
      ab = bracketingHD((imf::ImplicitFunction)imMDR, params, ff_MDR_appr / 4.,
                        ff_bla * 2., max_iter, mute);
    }
    // finding the root (friction factor)
    ff = imf::brent((imf::ImplicitFunction)imMDR, params, ab.l_b, ab.u_b,
                    epsilon, max_iter, mute);
  }
  return ff;
}

//--------------------------------------------------
// Mixed model: Laminar / MDR approximate / MDR
double ffMDRmixed2(IFParametersHD &params) {
  const bool mute = true;

  const double Re_L = 1379.33, Re_H = 1966.03;

  double ff;

  double ff_lam = 64. / params.Re_num;
  double ff_bla = 0.316 * std::pow(params.Re_num, -0.25);
  double ff_MDR_appr = ffMDRappr(params);

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else if (params.Re_num < Re_H) {
    ff = ff_MDR_appr;
  } else {
    // bracketing
    imf::BracketS ab = bracketingHD((imf::ImplicitFunction)imMDR, params,
                                    il::min(ff_lam, ff_MDR_appr) / 4.,
                                    ff_bla * 2., max_iter, mute);
    // finding the root (friction factor)
    ff = imf::brent((imf::ImplicitFunction)imMDR, params, ab.l_b, ab.u_b,
                    epsilon, max_iter, mute);
  }
  return ff;
}

//--------------------------------------------------
// Mixed model: Laminar / MDR approximate / MDR
double ffMDRmixed3(IFParametersHD &params) {
  const bool mute = true;

  const double Re_L = 1510.08, Re_H = 14022.1;

  double ff;

  double ff_lam = 64. / params.Re_num;
  double ff_bla = 0.316 * std::pow(params.Re_num, -0.25);
  double ff_MDR_appr = ffMDRappr(params);

  if (params.Re_num < Re_L) {
    ff = ff_lam;
  } else if (params.Re_num < Re_H) {
    ff = ff_MDR_appr;
  } else {
    // bracketing
    imf::BracketS ab =
        bracketingHD((imf::ImplicitFunction)imMDR, params, ff_MDR_appr / 2.,
                     ff_bla * 2., max_iter, mute);
    // finding the root (friction factor)
    ff = imf::brent((imf::ImplicitFunction)imMDR, params, ab.l_b, ab.u_b,
                    epsilon, max_iter, mute);
  }
  return ff;
}

//////////////////////////////////////////////////////////////////////////
//        FLOW SOLVER ROUTINES
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// nodal (edge) conductivities based on Darcy friction factor Ff(Re_num, rough)
il::Array<double> edgeConductivitiesP0(
    hfp2d::WellMesh &w_mesh, il::Array<double> &velocity, hfp2d::Fluid &fluid,
    double (*ffFunction)(IFParametersHD &params)) {
  // input should be an array of integers (il::int_t),
  // with rows for inner edges, each column for the 2 adjacent elements,
  // and hydraulic_width array with hydraulic width of each element.
  // returns A dv/d(grad P) = A 2 * hd / (density * abs(velocity) * ff(Re, eps))

  IFParametersHD params;

  il::Array2D<il::int_t> edgecommon = w_mesh.edgeCommon();

  il::Array<double> h_width = w_mesh.hd();
  il::Array<double> rough = w_mesh.rough();

  il::Array<double> edge_conduct{edgecommon.size(0), 0.};

  // loop on the inner nodes (edges)
  for (il::int_t i = 0; i < edgecommon.size(0); i++) {
    // harmonic mean of hydraulic width, left & right of the i-th edge
    double hd_mean = 2. *
                     (h_width[edgecommon(i, 1)] * h_width[edgecommon(i, 0)]) /
                     (h_width[edgecommon(i, 1)] + h_width[edgecommon(i, 0)]);

    double A_mean = hd_mean * hd_mean * il::pi / 4;

    // Reynolds number
    params.Re_num = hd_mean * fluid.fluidDensity() * std::fabs(velocity[i]) /
                    fluid.fluidViscosity();

    // dimensionless surface roughness (harmonic mean)
    params.rough = 2. * (rough[edgecommon(i, 0)] * rough[edgecommon(i, 1)]) /
                   (rough[edgecommon(i, 0)] + rough[edgecommon(i, 1)]);

    // Poiseuille (~zero Reynolds) switch
    if (params.Re_num < epsilon) {
      double poiseuille_ct = 1. / (32. * fluid.fluidViscosity());

      // i-th edge conductivity (Poiseuille law)
      edge_conduct[i] = poiseuille_ct * hd_mean * A_mean;  // cube?
    } else {
      // friction factor (Darcy)
      double f_f = ffFunction(params);

      // i-th edge conductivity
      edge_conduct[i] = 2. * hd_mean * A_mean / f_f / fluid.fluidDensity() /
                        std::fabs(velocity[i]);
    }
  }

  return edge_conduct;
}

////////////////////////////////////////////////////////////////////////////////
// Finite Volume (Finite Difference) Matrix L built from edge conductivities
// should return a sparse matrix....
// for now, for quick debug assemble a dense one
il::Array2D<double> buildWellFiniteDiffP0(hfp2d::WellMesh &w_mesh,
                                          il::Array<double> &edg_cond,
                                          double coef) {
  // Inputs:
  // wellMesh     :: Mesh object
  // edg_cond :: node (edge) Darcy conductivities array
  // coef     :: the factor of all entries (typically, the time-step)

  // checks here...

  il::Array2D<il::int_t> edgecommon = w_mesh.edgeCommon();

  // todo: L would better be a sparse matrix.....
  il::Array2D<double> L{w_mesh.numberOfElts(), w_mesh.numberOfElts(), 0.};

  il::int_t er, el, ml = edgecommon.size(0);

  double hi, dv;

  // loop on edges (connecting 2 elements each)
  for (il::int_t i = 0; i < ml; i++) {
    er = edgecommon(i, 0);
    el = edgecommon(i, 1);
    hi = (w_mesh.eltSize(er) + w_mesh.eltSize(el)) / 2.;
    // cross-sections, right & left of the edge (node)
    dv = coef * edg_cond[i] / hi;  //
    // diagonal entries
    L(er, er) += dv;
    L(el, el) += dv;
    // other entries
    L(er, el) -= dv;
    L(el, er) -= dv;
  }
  return L;  // such that L.P is Div.q in FV sense
}

////////////////////////////////////////////////////////////////////////////////
// function for current cell (element) volume -> returning a vector
il::Array<double> wellVolumeCompressibilityP0(
    hfp2d::WellMesh &mesh, hfp2d::Fluid &fluid,
    il::Array<double> &hydraulic_diameter) {
  il::Array<double> volume{mesh.numberOfElts(), 0.};
  // Area \times c_f \times h_i
  for (il::int_t e = 0; e < mesh.numberOfElts(); e++) {
    volume[e] = mesh.eltSize(e) * il::pi / 4. * hydraulic_diameter[e] *
                hydraulic_diameter[e] * fluid.fluidCompressibility();
  }
  return volume;
}

////////////////////////////////////////////////////////////////////////////////
// Solver of the well Hydrodynamics over a time-step with given in and out flow
hfp2d::WellSolution wellFlowSolverP0(
    hfp2d::WellSolution &well_soln, hfp2d::WellMesh &w_mesh,
    hfp2d::WellInjection &w_inj, hfp2d::Sources &out_flow,
    double (*ffFunction)(IFParametersHD &), double timestep,
    hfp2d::SimulationParameters &simul_params, bool mute, hfp2d::Fluid &fluid) {
  // Solver for wellbore Hydrodynamics over a time step.
  // (solves for fluid pressure over a time step)
  // PICARD / Fixed Pt Iterations SCHEME
  // Wellbore: fixed wellMesh, P0 elements
  //
  // soln         :: solution object at time tn
  //                 (converged solution)
  // w_mesh :: well mesh
  // w_inj :: well parameters description for in and outflow
  // fluid        :: Fluid properties object
  // ffFunction   :: Friction factor dependency on Reynolds number etc.
  // timestep     :: double, current size of the time step
  //
  // returns a Solution at t_n+1 object
  //

  il::int_t n_elts = w_mesh.numberOfElts();

  // call it here just in case
  // w_mesh.setNodeAdjElts();

  il::Array2D<il::int_t> edgecommon = w_mesh.edgeCommon();
  il::int_t n_s_e = edgecommon.size(0);

  //  WellInjection w_inj = well_soln.wellInjection();

  // Well plug location: element
  // il::int_t plug_location_e = w_inj.plugLocation();
  // element (farthest pressurized element)
  // todo: review the well plug handling -> most probably not needed.

  // remember, it's a P0 scheme
  il::int_t tot_dofs = n_elts;  // il::min(n_elts, plug_location_e + 1);

  // pressure @ prev. time step
  il::Array<double> p_n = well_soln.pressure();

  il::Array<double> p_hs = well_soln.hydrostaticPressure();

  il::Array<double> p_tilde_n = p_n;
  for (il::int_t i = 0; i < p_n.size(); i++) {
    p_tilde_n[i] -= p_hs[i];
  }
  // we need a copy to iterate as we use both
  il::Array<double> p_tilde_iter = p_tilde_n;

  // INITIALIZE the tangent matrix
  // todo: it would better be a sparse matrix.....
  il::Array2D<double> matrix_Xi{tot_dofs, tot_dofs, 0.};

  // INITIALIZE the RHS vector
  il::Array<double> vector_RS{tot_dofs, 0.};

  // velocity
  il::Array<double> v_k = well_soln.velocity();

  // INITIALIZE hydraulic width, sin(inclination), vectors of rel-errors
  // todo change to hyd_radius
  il::Array<double> hyd_diam = w_mesh.hd(), stat_p{n_elts, 0.},
                    err_Dv{n_s_e, 1.}, err_Dp{n_elts, 1.};

  // INITIALIZE the part of the tangent matrix matrix_Xi that won't change
  // during iterations
  il::Array<double> cfV = wellVolumeCompressibilityP0(w_mesh, fluid, hyd_diam);

  il::Array<double> DX_k{tot_dofs, 0.}, DX_k_1{tot_dofs, 0.};

  il::Array<double> Residuals{tot_dofs, 1.};
  double res_norm = 0.;

  double beta_relax = simul_params.ehl_relaxation;
  il::int_t k = 0;

  //-------------------------------------
  // Fixed Point Iteration Solver.
  while ((k < simul_params.ehl_max_its) &&
         (il::norm(err_Dp, il::Norm::L2) > simul_params.ehl_tolerance) &&
         (il::norm(err_Dv, il::Norm::L2) > simul_params.ehl_tolerance)) {
    k++;

    // calculate node (edge) conductivities with current values of velocities.
    il::Array<double> edg_cond =
        edgeConductivitiesP0(w_mesh, v_k, fluid, ffFunction);

    // calculate the updated Finite Diff Matrix
    // todo: it would better be a sparse matrix
    il::Array2D<double> L = buildWellFiniteDiffP0(w_mesh, edg_cond, timestep);

    il::Array<double> LdotPn = il::dot(L, p_tilde_n);

    // construct the tangent system (matrix + RHS)
    // with the effect of L (conductivity),
    // the rest being constant (compressibility + hydrostatic P)
    for (il::int_t j = 0; j < tot_dofs; j++) {
      for (il::int_t i = 0; i < tot_dofs; i++) {
        if (i == j) {
          matrix_Xi(i, j) = L(i, j) + cfV[i];
        } else {
          matrix_Xi(i, j) = L(i, j);
        }
      }
      vector_RS[j] = -LdotPn[j];
    }

    // don't forget to add source/sinks/plug to the RHS
    // source -> 0-th cell
    vector_RS[0] += w_inj.wellInjRate() * timestep;
    // sinks -> all HFs
    for (il::int_t i = 0; i < out_flow.SourceElt().size(); i++) {
      //      // check if the HF is above the plug
      //      if (w_inj.hfLocation(i) <= plug_location_e) {
      // add outflow to the RHS
      vector_RS[out_flow.SourceElt(i)] -= out_flow.InjectionRate(i) * timestep;
      //      }
    }
    // plug at an edge (node) -> zero flux through the edge
    // todo: review the well plug handling, plug not really needed, mesh can be
    // adjusted in inputs.

    // compute current residuals w.r.to previous iteration solution
    // (before solution of the new iteration)
    Residuals = vector_RS;
    il::blas(1.0, matrix_Xi, DX_k, -1.0, il::io, Residuals);
    res_norm = il::norm(Residuals, il::Norm::L2);

    // Solve the tangent system - solving for pressure increment
    il::Status status{};
    // use a direct solver
    DX_k = il::linearSolve(matrix_Xi, vector_RS, il::io, status);
    // il::LU<il::Array2D<double>> lu_decomposition(matrix_Xi, il::io, status);
    // if (!status.ok()) {
    //     // The matrix is singular to the machine precision.
    //     // You should deal with the error.
    // }
    // double cnd = lu_decomposition.conditionNumber(il::Norm::L2, );
    // std::cout << cnd << std::endl;
    // status.abortOnError();
    // DX_k = lu_decomposition.solve(vector_RS);
    status.ok();

    // under-relaxation of the solution
    il::blas((1. - beta_relax), DX_k_1, beta_relax, il::io, DX_k);

    // compute relative difference between sub-its
    // (relative error in pressure)
    for (il::int_t i = 0; i < n_elts; i++) {
      err_Dp[i] = std::fabs((DX_k[i] - DX_k_1[i]) / DX_k[i]);
    }

    // update pressure - p_hydro at current time step
    for (il::int_t i = 0; i < n_elts; i++) {
      p_tilde_iter[i] = p_tilde_n[i] + DX_k[i];
    }

    // update velocity according to pressure-p_hydro
    // loop through inner edges (shared by cells)
    for (il::int_t i = 0; i < n_s_e; i++) {
      // harmonic mean of hydraulic width, left & right of the i-th edge
      // this should be done in WellMesh
      double hd_mean =
          2. * (hyd_diam[edgecommon(i, 1)] * hyd_diam[edgecommon(i, 0)]) /
          (hyd_diam[edgecommon(i, 1)] + hyd_diam[edgecommon(i, 0)]);

      double hi = (w_mesh.eltSize(edgecommon(i, 0)) +
                   w_mesh.eltSize(edgecommon(i, 1))) /
                  2.;

      // flux through i-th edge (node)
      double v_new = edg_cond[i] * (p_tilde_iter[edgecommon(i, 0)] -
                                    p_tilde_iter[edgecommon(i, 1)]) /
                     hi;

      // divide by cross-sections to get av. velocity
      v_new /= il::pi / 4. * hd_mean * hd_mean;

      // compute relative difference between sub-its
      err_Dv[i] = std::fabs((v_k[i] - v_new) / v_new);

      // update velocities
      v_k[i] = v_new;
    }

    // old is new
    DX_k_1 = DX_k;

    if (!mute) {
      std::cout << " its. " << k
                << " rel. err. dp: " << il::norm(err_Dp, il::Norm::L2)
                << std::endl;
    };
  }  // end of iterative loop

  if (!mute) {
    std::cout << " end of Picard Scheme for Well flow, after " << k << " its. "
              << " rel. err. dp: " << il::norm(err_Dp, il::Norm::L2)
              << " rel. err. dv: " << il::norm(err_Dv, il::Norm::L2)
              << " norm of residuals: " << res_norm << std::endl;
  };

  // update total pressure
  for (il::int_t i = 0; i < p_n.size(); i++) {
    p_n[i] = p_tilde_iter[i] + p_hs[i];
  }

  // compute Reynolds number
  il::Array<double> Re(v_k.size(), 0.);
  for (il::int_t i = 0; i < n_s_e; i++) {
    // harmonic mean of hydraulic width, left & right of the i-th edge
    // this should be done in WellMesh
    double hd_mean = 2. *
                     (hyd_diam[edgecommon(i, 1)] * hyd_diam[edgecommon(i, 0)]) /
                     (hyd_diam[edgecommon(i, 1)] + hyd_diam[edgecommon(i, 0)]);

    Re[i] = hd_mean * fluid.fluidDensity() * std::fabs(v_k[i]) /
            fluid.fluidViscosity();
  }

  return hfp2d::WellSolution(
      well_soln.time() + timestep, timestep, p_hs, p_n, v_k, Re, out_flow, k,
      il::norm(err_Dp, il::Norm::L2), il::norm(err_Dv, il::Norm::L2));
}

////////////////////////////////////////////////////////////////////////////////
}