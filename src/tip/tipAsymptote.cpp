//
// This file is part of HFPx2D
//
// Created by D.Nikolski on 10/2/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <iostream>
//#include <cmath>
#include <il/math.h>
//#include <il/StaticArray.h>
#include "tipAsymptote.h"

// fracture tip asymptote(s) inversion
namespace tip {

////////////////////////////////////////////////////////////////////////////////
// scaling (see Dontsov & Peirce 2015, 2017)
double k_P(double k1c) { return 4.0 * std::pow(2.0 / il::pi, 0.5) * k1c; }

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
  return 16.0 * (1.0 - 3.0 * d) * std::tan(1.5 * il::pi * d) / 3.0 / d /
         (2.0 - 3.0 * d);
}

double b_12(double d) { return c_2(d) / c_1(d); }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// approximate solution of integral eqn (see Dontsov & Peirce 2015, 2017)
double effe(double k_h, double b_h, double c1) {
  double l_h = std::log((1.0 + b_h) / (k_h + b_h));
  double k_h_2 = k_h * k_h;
  double b_h_2 = b_h * b_h;
//  double cf1 = (1.0 - k_h) / b_h;
//  double cf2 = (1.0 + k_h) / b_h;
//  double cf3 = (1.0 + k_h + k_h_2) / b_h_2;
//  return (cf1 * (cf3 / 3.0 - cf2 / 2.0 + 1.0) - l_h) / c1;
  double k_h_3 = k_h * k_h_2;
  double b_h_3 = b_h * b_h_2;
  return (1.0 - k_h_3 - 1.5 * b_h * (1.0 - k_h_2)
          + 3.0 * b_h_2 * (1.0 - k_h)
          - 3.0 * b_h_3 * l_h) /
         (3.0 * c1);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// various scaled tip asymptote approximations
// zero-order approximation for universal tip asymptote
double g_un_0(double k_h, double c_h) {
  double b_h = con_mc * c_h;
  return effe(k_h, b_h, con_m);
//  return std::pow(effe(k_h, b_h, con_m), al_exp)
//         * std::pow(b_h, al_exp * be_exp);
}

// 1st order delta-correction for universal tip asymptote
double delta_c(double k_h, double c_h) {
  double b_h = con_mc * c_h;
  return con_m * (1.0 + b_h) * std::pow(b_h, 3.0 - be_exp)
                     * effe(k_h, b_h, con_m);
}

double g_un_1(double k_h, double c_h) {
  double delta = con_m * (1.0 + con_mc * c_h) * g_un_0(k_h, c_h);
//  double delta = delta_c(k_h, c_h);
  double c_1_d = c_1(delta);
  double b_12_d = b_12(delta);
  double b_h_d = c_h * b_12_d;
  return effe(k_h, b_h_d, c_1_d);
//  return std::pow(effe(k_h, b_h_d, c_1_d), al_exp)
//         * std::pow(b_h_d, al_exp * be_exp);
}

// zero-order approximation for k-m edge solution (zero leak-off)
double g_km_0(double k_h) { return (1.0 - std::pow(k_h, 3.0)) / beta_m_3; }

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
// residual functions to find the root (distance to the tip)
// zero-order approximation
double res_u_0(double s, TipParameters &taParam) {
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
    // use universal tip asymptote
    double c_h = c_H(taParam.cl, std::max(v, m_tol), taParam.wa, s);
    double s_u = g_un_0(k_h, c_h);
//    return s_u - s_h;
    double b_h = con_mc * c_h;
    return s_u - std::pow(s_h / std::pow(b_h, 3.0 - be_exp), al_exp);
  }
}

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
        double c_h = c_H(taParam.cl, std::max(v, m_tol), taParam.wa, s);
        if (isMisbehaving(s, taParam)) {
          // use k-m~ asymptote
          double s_l = g_kc_0(k_h, c_h);
          return s_l - s_h;
        } else {
          // use universal tip asymptote
          double s_u = g_un_0(k_h, c_h);
          return s_u - s_h;
//      double b_h = con_mc * c_h;
//      return s_u - std::pow(s_h / std::pow(b_h, 3.0 - be_exp), al_exp);
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
      // use universal tip asymptote
      double s_u = g_un_1(k_h, c_h);
      return s_u - s_h;
//      double delta = delta_c(k_h, c_h);;
//      double b_12_d = b_12(delta);
//      double b_h_d = c_h * b_12_d;
//      return s_u - std::pow(s_h / std::pow(b_h_d, 3.0 - be_exp), al_exp);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// criterion for unstable residual function behavior (chi > chi_c)
// todo: with modified residual function, it might be unnecessary
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
  return (k_p * std::pow(taParam.s0, 0.5) <= taParam.e_p * taParam.wa);
}

// overload for s as an independent parameter
bool isPropagating(double s, TipParameters &taParam) {
  double k_p = k_P(taParam.k1c);
  return (k_p * std::pow(s, 0.5) <= taParam.e_p * taParam.wa);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// bracketing the root (tip position)
// simple way (works for modified residual function for high chi)
// todo: bracketing for zero toughness - now fixed upper bound is used
imf::BracketS bracketingTip(imf::ImplicitFunction resF,
                            TipParameters &taParam,
                            double up_bound,
                            int maxIter, bool mute) {
  // resF :: is the tip asymptote function to bracket
  // taParam :: is the tipParameters structure containing all necessary
  // parameters
  // up_bound :: is the initial upper bound (for the zero toughness case)
  // maxIter :: is the number of maximum number of iteration for the root
  // bracketing scheme
  // mute :: boolean - true for muting output.

  imf::BracketS ab;

  if (taParam.k1c == 0.0) {
    ab.u_b = up_bound;
  } else {
    // inverse K-asymptote as upper bound
    double k_p = k_P(taParam.k1c);
    double s_k = std::pow((taParam.e_p * taParam.wa) / k_p, 2);
    // take K-solution as upper bound
    ab.u_b = s_k;
  }

  double fb = resF(ab.u_b, taParam);

  // take previous tip position + inverse asymptote at sqrt(machine precision)
  // as lower bound
  double ds = std::sqrt(m_tol)
            * 27.0 * (4.0 * taParam.cl * taParam.cl)
            * taParam.s0 * taParam.dt
            / (taParam.wa * taParam.wa);
  ab.l_b = taParam.s0 + ds;
  double fa = resF(ab.l_b, taParam);

  // search for a better lower bound if signs of fa and fb are equal
  int count = 0;
  double wt1 = 2.0, wt2 = 1.0;
  double mid = ab.u_b;
  while (fa * fb > 0.0 && count < maxIter) {
    // take an intermediate point for lower bound
    mid = (wt1 * ab.l_b + wt2 * mid) / (wt1 + wt2);
    ab.l_b = mid;
    fa = resF(ab.l_b, taParam);
    ++count;
  }

  if (count >= maxIter && fa * fb > 0.0) {
    // if root isn't bracketed
    if (!mute) {
      std::cout << "The root of f(s) can't be bracketed in " << count
                << " iterations" << std::endl;
    }
  }

  if (!mute) {
    std::cout << "The root of f(s) is bracketed between " << ab.l_b << " and "
              << ab.u_b << " in " << count << " iterations" << std::endl;
  }
  return ab;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// for Brent's root finding method, see src/util/RootFinder.cpp
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// finding the distance to the tip & velocity after a time step dt
void tipInversion(imf::ImplicitFunction resF,
                  TipParameters &tipS,
                  double up_bound,
                  double epsilon, int maxIter, bool mute) {
  // resF ::  the corresponding residual function
  //          resF(x, params)  defining the tip asymptote.
  // tipS :: the tip parameters containing all the necessacry parameters, this
  // structure
  //          is UPDATED HERE, i.e. tipS.st contains the new distance ribbon -
  //          tip
  //                           and tipS.vt the corresponding tip velocity
  // up_bound :: upper bound for the new distance ribbon - tip
  // epsilon :: tolerance for finding the zero of the tip inversion w=resF(x,
  // params)
  // maxIter :: maximum number of Iterations of the Brent scheme
  // mute :: boolean for muting outputs

  tipS.vt = std::max(tipS.vt, m_tol);  // previous tip velocity

  // todo: triggering time for Carter leak-off

  // check K vs K1c (taParam.k1c)
  bool isProp = isPropagating(tipS);

  if (!isProp) {
    // not propagating
    tipS.vt = 0.0;
    tipS.st = tipS.s0;
  } else {
    // fully coupled iteration for v and s by default
    // bracketing the tip position for fine root finding
    imf::BracketS ab = bracketingTip((imf::ImplicitFunction)resF,
                                     tipS, up_bound, 30, mute);
    // fine root finding (adjusting tip position)
    double s_a = imf::brent((imf::ImplicitFunction)resF, tipS,
                            ab.l_b, ab.u_b,
                            epsilon, maxIter, mute);
    // adjusting tip velocity
    double v_a = (s_a - tipS.s0) / tipS.dt;
    tipS.vt = v_a;
    tipS.st = s_a;
  }
  // todo: Carter tip leak-off related
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// moments (tip volume)
double deltaP(double k_h, double c_h, double p) {
  double delta;
  if (c_h == 0) {
    delta = con_m * g_km_0(k_h);
    return (1.0 - p + p * g_km_0(k_h)) * delta;
  } else {
    delta = con_m * (1.0 + con_mc * c_h) * g_un_0(k_h, c_h);
//    delta = delta_c(k_h, c_h);
    return (1.0 - p + p * g_un_0(k_h, c_h)) * delta;
  }
}


double moment0(TipParameters &taParam) {
  double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, taParam.st);
  double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, taParam.st);
  double delta_p = deltaP(k_h, c_h, 0.377);
  return (2.0 * taParam.wa * taParam.st) / (3.0 + delta_p);
}

double moment1(TipParameters &taParam) {
  double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, taParam.st);
  double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, taParam.st);
  double delta_p = deltaP(k_h, c_h, 0.26);
  return (2.0 * taParam.wa * taParam.st * taParam.st) / (5.0 + delta_p);
}


// overload for s as an independent parameter
double moment0(double s, TipParameters &taParam) {
  if (isPropagating(taParam)) {
    double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
    double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
    double delta_p = deltaP(k_h, c_h, 0.377);
    return (2.0 * taParam.wa * s) / (3.0 + delta_p);
  } else {
    double K1prime = taParam.e_p * taParam.wa / std::pow(taParam.s0, 0.5);
    return (2./3.) * K1prime / taParam.e_p * std::pow(s, 1.5);
  }
}

double moment1(double s, TipParameters &taParam) {
  if (isPropagating(taParam)) {
    double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
    double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
    double delta_p = deltaP(k_h, c_h, 0.26);
    return (2.0 * taParam.wa * s * s) / (5.0 + delta_p);
  } else
  {
    double K1prime = taParam.e_p * taParam.wa / std::pow(taParam.s0, 0.5);
    return (2./5.) * K1prime / taParam.e_p * std::pow(s, 2.5);
  }
}

////////////////////////////////////////////////////////////////////////////////

// todo: Carter tip leak-off volume

}
