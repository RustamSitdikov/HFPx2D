//
// This file is part of HFPx2D
//
// Created by nikolski on 10/2/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_TIP_ASYMPTOTE_H
#define HFPX2D_TIP_ASYMPTOTE_H

#include <cmath>
// #include <il/math.h>
#include <il/StaticArray.h>

// fracture tip asymptote inversion
namespace tip {
    // constants
    const double beta_m = std::pow(2.0, 1.0 / 3.0) * std::pow(3.0, 5.0 / 6.0);
    //
    const double beta_m_3 = std::pow(beta_m, 3);
    //
    const double beta_c = 4.0 / std::pow(15.0 * (std::pow(2.0, 0.5) - 1.0), 0.25);
    //
    const double beta_c_4 = std::pow(beta_c, 4);
    //
    const double con_m = beta_m_3 / 3.0;
    // 10.392304845
    const double con_mc = 0.75 * beta_c_4 / beta_m_3;
    // 0.9911799823

    // tolerance, close to machine precision
    const double m_tol = 2.221e-016;

    // fracture tip parameters
    struct TAParam {
        double k1c; // SIF // K' = 4.0 * std::pow(2.0 / il::pi, 0.5) * K1c
        double e_p; // Plane strain modulus = youngPS_ = E / (1.0 - nu*nu)
        double cl; // Carter leak-off coefficient; C' = 2.0 * Cl is used
        double mu; // Fluid viscosity
        double wa; // Opening at the ribbon cell
        double st; // Distance from the ribbon cell center to the tip
        double vt; // Tip velocity
        double m0; // 0th moment (volume of the tip)
        double m1; // 1st moment (volume of the tip)
        double ta; // Triggering time for tip leak-off
        double cv; // Tip leak-off volume
    };

    // input
    struct TAInParam {
        TAParam taPrev; // Previous tip state
        double wa; // New opening at the ribbon cell
        double dt; // Time step
    };

    // scaling
    double k_H(double k1c, double e_p, double w, double s);
    double c_H(double cl, double v, double w, double s);
    double s_H(double mu, double e_p, double v, double w, double s);

    // old scaling
    double s_T(double mu, double k1c, double e_p, double v, double s);

    // auxiliary functions
    double c_1(double d);
    double c_2(double d);
    double b_12(double d);

    // approximate solution of integral eqn (see Dontsov & Peirce 2015, 2017)
    double effe(double k_h, double b_h, double c1_h);

    // approximations for k-m edge
    double g_km_0(double s_t); // zero-order approximation
    double g_km_1(double s_t); // 1st order delta-correction

    // approximations for k-m~ edge
    double g_kc_0(double chi, double s_t); // zero-order approximation
    double g_kc_1(double chi, double s_t); // 1st order delta-correction

    // approximations for universal tip asymptote (functions to minimize)
    double g_0(double k_h, double c_h); // zero-order approximation
    double g_1(double k_h, double c_h); // 1st order delta-correction

    // moments
    double deltaP(double k_h, double c_h, double p);
    double moment0(TAParam &taParam);
    double moment1(TAParam &taParam);
    double moment0(double s, TAParam &taParam);
    double moment1(double s, TAParam &taParam);

    // (virtual) residual function of distance to minimize (set to zero)
    typedef double (*ResFun)(double s, TAParam &taParam);

    // possibly, viscosity & leak-off asymptotes to be added
    // zero-order approximation
    double res_g_0(double s, TAParam &taParam);
    double res_g_0_s(double s, TAParam &taParam);
    // 1st order delta-correction
    double res_g_1(double s, TAParam &taParam);
    double res_g_1_s(double s, TAParam &taParam);

    // checking propagation criterion
    bool isPropagating(TAParam &taParam);
    bool isPropagating(double s, TAParam &taParam); // overload just in case...

    // bracketing the tip
    il::StaticArray<double, 2> bracket(ResFun resF,
                                       TAParam &taParam,
                                       int maxIter);

    // Brent root finder
    double brent(tip::ResFun fun,
                 tip::TAParam &params,
                 double a0, double b0,
                 double epsilon, int maxIter);

    // finding the distance to the tip, moments, etc.
    TAParam propagateTip(ResFun resF,
                         TAInParam &taIn,
                         double epsilonS, int maxIterS,
                         double epsilonV, int maxIterV,
                         bool mute);

    TAParam propagateTip(ResFun resF,
                         TAInParam &taIn,
                         double epsilon, int maxIter,
                         bool mute);

}

#endif //HFPX2D_TIP_ASYMPTOTE_H
