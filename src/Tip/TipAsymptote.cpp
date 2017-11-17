//
// This file is part of HFPx2D
//
// Created by D.Nikolski on 10/2/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cmath>
#include <il/math.h>
#include <il/StaticArray.h>
#include "TipAsymptote.h"
#include <iostream>

// fracture Tip asymptote(s) inversion
namespace tip {

////////////////////////////////////////////////////////////////////////////////
// scaling (see Dontsov & Peirce 2015, 2017)
    double k_P(double k1c) {
      return 4.0 * std::pow(2.0 / il::pi, 0.5) * k1c;
    }

    double k_H(double k1c, double e_p, double w, double s) {
      double k_p = k_P(k1c);
      return (k_p * std::pow(s, 0.5)) / (e_p * w);
    }

    double c_H(double cl, double v, double w, double s) {
      double c_p = 2.0 * cl;
      return (2.0 * c_p * std::pow(s, 0.5)) / (std::pow(v, 0.5) * w);
    }

    double s_H(double mu, double e_p, double v, double w, double s) {
      double mu_p = 12.0 * mu;
      return (mu_p * v * s * s) / (e_p * std::pow(w, 3));
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// auxiliary functions
    double c_1(double d) {
      return 4.0 * (1.0 - 2.0 * d) * std::tan(il::pi * d) / d / (1.0 - d);
    }

    double c_2(double d) {
      return 16.0 * (1.0 - 3.0 * d) * std::tan(1.5 * il::pi * d) /
             3.0 / d / (2.0 - 3.0 * d);
    }

    double b_12(double d) {
      return c_2(d) / c_1(d);
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// approximate solution of integral eqn (see Dontsov & Peirce 2015, 2017)
    double effe(double k_h, double b_h, double c1) {
      double k_h_2 = k_h * k_h, k_h_3 = k_h * k_h_2;
      double b_h_2 = b_h * b_h, b_h_3 = b_h * b_h_2;
      double l_h = std::log((1.0 + b_h) / (k_h + b_h));
      return (1.0 - k_h_3 - 1.5 * b_h * (1.0 - k_h_2) +
              3.0 * b_h_2 * (1.0 - k_h) - 3.0 * b_h_3 * l_h) / (3.0 * c1);
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// various scaled Tip asymptote approximations
// zero-order approximation for universal Tip asymptote
    double g_un_0(double k_h, double c_h) {
      return effe(k_h, con_mc * c_h, con_m);
    }

// 1st order delta-correction for universal Tip asymptote
    double g_un_1(double k_h, double c_h) {
      double delta = con_m * (1.0 + con_mc * c_h) * g_un_0(k_h, c_h);
      double c_1_d = c_1(delta);
      double b_12_d = b_12(delta);
      return effe(k_h, c_h * b_12_d, c_1_d);
    }

// zero-order approximation for k-m edge solution (zero leak-off)
    double g_km_0(double k_h) {
      return (1.0 - std::pow(k_h, 3.0)) / beta_m_3;
    }

// 1st order delta-correction for k-m edge solution (zero leak-off)
    double g_km_1(double k_h) {
      double delta = con_m * g_km_0(k_h);
      double c_1_d = c_1(delta);
      return (1.0 - std::pow(k_h, 3.0)) / 3.0 / c_1_d;
    }

// zero-order approximation for k-m~ edge solution
    double g_kc_0(double k_h, double c_h) {
      return (1.0 - std::pow(k_h, 4.0)) / beta_c_4 / c_h;
    }

// 1st order delta-correction for k-m~ edge solution
    double g_kc_1(double k_h, double c_h) {
      double delta = con_m * (1.0 + con_mc * c_h) * g_kc_0(k_h, c_h);
      double c_2_d = c_2(delta);
      return (1.0 - std::pow(k_h, 4.0)) / 4.0 / c_2_d / c_h;
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// residual functions to find the root (distance to the Tip)
// (modified to overcome misbehavior at high chi values)
// zero-order approximation
    double res_u_0_m(double s, TipParameters &taParam) {
      // assume fully coupled iteration
      double v = (s - taParam.s0) / taParam.dt;
      // cutting out negative values (remove zero if necessary)
      v = std::max(v, 0.0 * m_tol);
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
      double s_h = s_H(taParam.mu, taParam.e_p, v, taParam.wa, s);
      if (taParam.cl == 0.0) {
        // use k-m approximate solution for zero leak-off
        double s_m = g_km_0(k_h);
        return s_m - s_h;
      } else {
        double c_h = c_H(taParam.cl, std::max(v, m_tol), taParam.wa,
                         s);
        if (isMisbehaving(s, taParam)) {
          // use k-m~ asymptote
          double s_l = g_kc_0(k_h, c_h);
          return s_l - s_h;
        } else {
          // use universal Tip asymptote
          double s_u = g_un_0(k_h, c_h);
          return s_u - s_h;
        }
      }
    }

// 1st order delta-correction
    double res_u_1_m(double s, TipParameters &taParam) {
      // assume fully coupled iteration
      double v = (s - taParam.s0) / taParam.dt;
      // cutting out negative values (remove zero if necessary)
      v = std::max(v, m_tol);
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
      double s_h = s_H(taParam.mu, taParam.e_p, v, taParam.wa, s);
      if (taParam.cl == 0.0) {
        // use k-m approximate solution for zero leak-off
        double s_m = g_km_1(k_h);
        return s_m - s_h;
      } else {
        double c_h = c_H(taParam.cl, v, taParam.wa, s);
        if (isMisbehaving(s, taParam)) {
          // use k-m~ asymptote
          double s_l = g_kc_1(k_h, c_h);
          return s_l - s_h;
        } else {
          // use universal Tip asymptote
          double s_u = g_un_1(k_h, c_h);
          return s_u - s_h;
        }
      }
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// criterion for unstable residual function behavior (chi > chi_c)
    bool isMisbehaving(double s, TipParameters &taParam) {
      double k_p = k_P(taParam.k1c);
      return (std::pow(s - taParam.s0, 0.5) * k_p * chi_c <
              4.0 * std::pow(taParam.dt, 0.5) * taParam.cl * taParam.e_p);
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// checking propagation criterion
    bool isPropagating(TipParameters &taParam) {
      double k_p = k_P(taParam.k1c);
      return (k_p * std::pow(taParam.s0, 0.5)
              <= taParam.e_p * taParam.wa);
    }

// overload for s as an independent parameter
    bool isPropagating(double s, TipParameters &taParam) {
      double k_p = k_P(taParam.k1c);
      return (k_p * std::pow(s, 0.5) <= taParam.e_p * taParam.wa);
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// bracketing the root (Tip position)
// simple way (works for modified residual function for high chi)
// todo: bracketing for zero toughness - now fixed upper bound is used
    il::StaticArray<double, 2> bracket
            (ResidualFunction resF,
             TipParameters &taParam,
             double up_bound,
             int maxIter,
             bool mute) {
      il::StaticArray<double, 2> ab;

      if (taParam.k1c == 0.0) {
        ab[1] = up_bound;
      } else {
        // inverse K-asymptote
        double k_p = k_P(taParam.k1c);
        double s_k = std::pow((taParam.e_p * taParam.wa) / k_p, 2);
        // take K-solution as upper bound
        ab[1] = s_k;
      }

      double fb = resF(ab[1], taParam);

      // take previous Tip position as lower bound
      ab[0] = taParam.s0; // + m_tol;
      double fa = resF(ab[0], taParam);

      // search for a better lower bound if signs of fa and fb are equal
      int count = 0;
      double wt1 = 2.0, wt2 = 1.0;
      double mid = ab[1];
      while (fa * fb > 0.0 && count < maxIter) {
        // take an intermediate point for lower bound
        mid = (wt1 * ab[0] + wt2 * mid) / (wt1 + wt2);
        ab[0] = mid;
        fa = resF(ab[0], taParam);
        ++count;
      }

      if (count >= maxIter && fa * fb > 0.0) {
        // if root isn't bracketed
        if (!mute) {
          std::cout << "The root of f(s) can't be bracketed in "
                    << count << " iterations" << std::endl;
        }
      }

      if (!mute){
        std::cout << "The root of f(s) is bracketed between "
                  << ab[0] << " and " << ab[1] << " in "
                  << count << " iterations" << std::endl;
      }
      return ab;
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Brent's root finding method
//
// Brent R.P. (1973),
// "Chapter 4: An Algorithm with Guaranteed Convergence
// for Finding a Zero of a Function",
// Algorithms for Minimization without Derivatives,
// Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2
//
    double brent
            (tip::ResidualFunction fun,
             tip::TipParameters &params,
             double a0, double b0,
             double epsilon, int maxIter) {
  //
  double a = a0;
  double b = b0;
  double fa = fun(a, params);  // calculated now to save function calls
  double fb = fun(b, params);  // calculated now to save function calls
  double fs = 0;               // initialize

  if (std::isnan(fa)) {
    std::cout << "f(a) isn't finite; a=" << a << std::endl;
    return -111;
  }

  if (std::isnan(fb)) {
    std::cout << "f(b) isn't finite; b=" << b << std::endl;
    return -111;
  }

  if (fa * fb > 0.0) {
    // throw exception if root isn't bracketed
    std::cout << "Signs of f(a) and f(b) must be opposite" << std::endl;
    return -11;
  }

  // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  if (std::abs(fa) < std::abs(b)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  // c now equals the largest magnitude of the lower and upper bounds
  double c = a;
  double fc = fa;  // precompute function evalutation for point c by assigning
                   // it the same value as fa
  bool mflag = true;  // boolean flag used to evaluate if statement later on
  double s = 0;       // Our Root that will be returned
  double d = 0;       // Only used if mflag is unset (mflag == false)

  for (unsigned int iter = 1; iter < maxIter; ++iter) {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b - a) < epsilon) {
      //                std::cout << "After " << iter
      //                          << " iterations the root is: "
      //                          << s << std::endl;
      return s;
    }  // end if

    if (fa != fc && fb != fc) {
      // use inverse quadratic interopolation
      s = (a * fb * fc / ((fa - fb) * (fa - fc))) +
          (b * fa * fc / ((fb - fa) * (fb - fc))) +
          (c * fa * fb / ((fc - fa) * (fc - fb)));
    } else {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }

    // check to see whether we can use the faster converging quadratic
    // && secant methods or if we need to use bisection
    if (((s < (3 * a + b) * 0.25) || (s > b)) ||
        (mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
        (!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
        (mflag && (std::abs(b - c) < epsilon)) ||
        (!mflag && (std::abs(c - d) < epsilon))) {
      // bisection method
      s = (a + b) * 0.5;

      mflag = true;
    } else {
      mflag = false;
    }

    fs = fun(s, params);  // calculate fs
    d = c;                // first time d is being used
    // (wasnt used on first iteration because mflag was set)
    c = b;    // set c equal to upper bound
    fc = fb;  // set f(c) = f(b)

    // fa and fs have opposite signs
    if (fa * fs < 0) {
      b = s;
      fb = fs;  // set f(b) = f(s)
    } else {
      a = s;
      fa = fs;  // set f(a) = f(s)
    }

    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);    // swap a and b
      std::swap(fa, fb);  // make sure f(a) and f(b) are correct after swap
    }

  }  // end of loop
  // error message
  std::cout << "The solution did not converge after " << maxIter
            << " iterations" << std::endl;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// finding the distance to the Tip & velocity after a time step dt
    TipParameters tipStep
            (ResidualFunction resF,
             TipParameters &taPrev,
             double dt, double wa,
             double up_bound,
             double epsilon, int maxIter,
             bool mute) {
      // Substitute the corresponding residual function
      // resF(x, params)
      // defining the Tip asymptote.

      TipParameters taNew;
      taNew.wa = wa; // new ribbon cell opening
      taNew.s0 = taPrev.st; // previous distance to the Tip
      taNew.vt = std::max(taPrev.vt, m_tol); // previous Tip velocity
      taNew.dt = dt; // taPrev.dt; // time step
      // copying other parameters
      taNew.k1c = taPrev.k1c;
      taNew.e_p = taPrev.e_p;
      taNew.cl = taPrev.cl;
      taNew.mu = taPrev.mu;
      // todo: triggering time for Carter leak-off

      // check K vs K1c (taParam.k1c)
      bool isProp = isPropagating(taNew);

      if (!isProp) {
        // not propagating
        taNew.vt = 0.0;
        taNew.st = taNew.s0;
      } else {
        // fully coupled iteration for v and s by default
        // bracketing the Tip position for fine root finding
        il::StaticArray<double, 2> ab =
                bracket(resF, taNew, up_bound, 30, mute);
        // fine root finding (adjusting Tip position)
        double s_a = brent(resF, taNew, ab[0], ab[1],
                           epsilon, maxIter);
        // adjusting Tip velocity
        double v_a = (s_a - taNew.s0) / taNew.dt;
        taNew.vt = v_a;
        taNew.st = s_a;
      }
      // todo: Carter Tip leak-off
      return taNew;
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// moments (Tip volume)
    double deltaP(double k_h, double c_h, double p) {
      double delta = con_m * (1.0 + con_mc * c_h) * g_un_0(k_h, c_h);
      return (1.0 - p + p * g_un_0(k_h, c_h)) * delta;
    }

    double moment0(TipParameters &taParam) {
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, taParam.st);
      double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, taParam.st);
      double delta_p = deltaP(k_h, c_h, 0.377);
      return (2.0 * taParam.wa * taParam.st)/(3.0 + delta_p);
    }

    double moment1(TipParameters &taParam) {
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, taParam.st);
      double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, taParam.st);
      double delta_p = deltaP(k_h, c_h, 0.26);
      return (2.0 * taParam.wa * taParam.st * taParam.st)/(5.0 + delta_p);
    }

// overload for s as an independent parameter
    double moment0(double s, TipParameters &taParam) {
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
      double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
      double delta_p = deltaP(k_h, c_h, 0.377);
      return (2.0 * taParam.wa * s)/(3.0 + delta_p);
    }

    double moment1(double s, TipParameters &taParam) {
      double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
      double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
      double delta_p = deltaP(k_h, c_h, 0.26);
      return (2.0 * taParam.wa * s * s)/(5.0 + delta_p);
    }
////////////////////////////////////////////////////////////////////////////////

// todo: Carter Tip leak-off

}
