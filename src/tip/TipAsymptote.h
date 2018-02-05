//
// This file is part of HFPx2D
//
// Created by D.Nikolski on 10/2/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_TIP_ASYMPTOTE_H
#define HFPX2D_TIP_ASYMPTOTE_H

#include <cmath>
// #include <il/math.h>
#include <il/StaticArray.h>

// fracture tip asymptote(s) inversion
namespace tip {
// constants
// const double beta_m = std::pow(2.0, 1.0 / 3.0) * std::pow(3.0, 5.0 / 6.0);
// = 3.1473451902649443;
const double beta_m_3 = 2.0 * std::pow(3.0, 2.5);
// = 31.176914536239792; // = std::pow(beta_m, 3);
// const double beta_c = 4.0 / std::pow(15.0 * (std::pow(2.0, 0.5) - 1.0),
// 0.25);
// = 2.5335594408265694;
const double beta_c_4 = 256.0 / 15.0 / (std::pow(2.0, 0.5) - 1.0);
// = 41.202578131167478; // = std::pow(beta_c, 4);
const double con_m = beta_m_3 / 3.0;
// = 10.392304845413264; // = c_3(1.0 / 3.0);
// const double con_c = beta_c_4 / 4.0;
// = ; // c_2(0.25);
const double con_mc = 0.75 * beta_c_4 / beta_m_3;
// = 0.99117998230567206; // = con_c / con_m;

// v tolerance (> machine precision)
const double m_tol = 2.221e-016;

// critical Chi (misbehaving residual functions)
const double chi_c = 3000.0;

// fracture tip parameters (input & output)
struct TipParameters {
  double k1c;  // SIF // K' = 4.0 * std::pow(2.0 / il::pi, 0.5) * K1c
  double e_p;  // Plane strain modulus = youngPS_ = E / (1.0 - nu*nu)
  double cl;   // Carter leak-off coefficient; C' = 2.0 * Cl is used
  double mu;   // Fluid viscosity
  double wa;   // Opening at the ribbon cell
  double s0;   // Previous distance from the ribbon cell center to the tip
  double dt;   // Time step
  double st;   // Distance from the ribbon cell center to the tip
  double vt;   // Tip velocity
};

// scaling (see Dontsov & Peirce 2015, 2017)
double k_P(double k1c);
double k_H(double k1c, double e_p, double w, double s);
double c_H(double cl, double v, double w, double s);
double s_H(double mu, double e_p, double v, double w, double s);

// auxiliary functions
double c_1(double d);
double c_2(double d);
double b_12(double d);

// approximate solution of integral eqn (see Dontsov & Peirce 2015, 2017)
double effe(double k_h, double b_h, double c1_h);

// various scaled tip asymptote approximations
// approximations for k-m edge
double g_km_0(double s_t);  // zero-order approximation
double g_km_1(double s_t);  // 1st order delta-correction

// approximations for k-m~ edge
double g_kc_0(double k_h, double c_h);  // zero-order approximation
double g_kc_1(double k_h, double c_h);  // 1st order delta-correction

// approximations for universal tip asymptote (functions to minimize)
double g_un_0(double k_h, double c_h);  // zero-order approximation
double g_un_1(double k_h, double c_h);  // 1st order delta-correction

// criteria for unstable residual function behavior
bool isMisbehaving(double s, TipParameters &taParam);

// checking propagation criterion
bool isPropagating(TipParameters &taParam);
// overload just in case...
bool isPropagating(double s, TipParameters &taParam);

// (virtual) residual function of distance to minimize (set to zero)
typedef double (*ResidualFunction)(double s, TipParameters &taParam);

// particular residual functions to find the root (distance to the tip)
// (modified to overcome misbehavior at high chi values)
// zero-order approximation
double res_u_0_m(double s, TipParameters &taParam);
// 1st order delta-correction
double res_u_1_m(double s, TipParameters &taParam);

// bracketing the tip
// (works for modified residual function for high chi)
il::StaticArray<double, 2> bracket(ResidualFunction resF,
                                   TipParameters &taParam, double up_bound,
                                   int maxIter, bool mute);

// Brent root finder
double brent(tip::ResidualFunction fun, tip::TipParameters &params, double a0,
             double b0, double epsilon, int maxIter);

// finding the distance to the tip & velocity after a time step dt
// bl - I don t really understand why dt and wa are inputs if they already are
// stored in taIn
void tipInversion(ResidualFunction resF, TipParameters &tipS, double up_bound,
                  double epsilon, int maxIter, bool mute);

// moments (volume of the tip)
double deltaP(double k_h, double c_h, double p);
double moment0(TipParameters &taParam);
double moment1(TipParameters &taParam);
double moment0(double s, TipParameters &taParam);
double moment1(double s, TipParameters &taParam);

}

#endif  // HFPX2D_TIP_ASYMPTOTE_H
