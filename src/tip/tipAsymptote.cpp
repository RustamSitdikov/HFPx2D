//
// This file is part of HFPx2D
//
// Created by nikolski on 10/2/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "tipAsymptote.h"
#include <il/StaticArray.h>
#include <il/math.h>
#include <cmath>
#include <iostream>

// fracture tip asymptote(s)
namespace tip {

// scaling
    double k_H(double k1c, double e_p, double w, double s) {
        double k_p = 4.0 * std::pow(2.0 / il::pi, 0.5) * k1c;
        return (k_p * std::pow(s, 0.5)) / (e_p * w);
    }

    double c_H(double cl, double v, double w, double s) {
        double c_p = 2.0 * cl;
        return (2.0 * c_p * std::pow(s, 0.5)) / (std::pow(v, 0.5) * w);
    }

    double s_H(double mu, double e_p, double v, double w, double s) {
        return (12.0 * mu * v * s * s) / (e_p * std::pow(w, 3));
    }

// old scaling
    double s_T(double mu, double k1c, double e_p, double v, double s) {
        double k_p = 4.0 * std::pow(2.0 / il::pi, 0.5) * k1c;
        return std::pow(s, 0.5) * (12.0 * mu * v * e_p * e_p) / std::pow(k_p, 3);
    }

// auxiliary functions
    double c_1(double d) {
        return 4.0 * (1.0 - 2.0 * d) * std::tan(il::pi * d) / d / (1.0 - d);
    }

    double c_2(double d) {
        return 16.0 * (1.0 - 3.0 * d) * std::tan(1.5 * il::pi * d) / 3.0 / d /
               (2.0 - 3.0 * d);
    }

    double b_12(double d) { return c_2(d) / c_1(d); }

// approximate solution of integral eqn (see Dontsov & Peirce 2015, 2017)
    double effe(double k_h, double b_h, double c1_h) {
        double k_h_2 = k_h * k_h, k_h_3 = k_h * k_h_2;
        double b_h_2 = b_h * b_h, b_h_3 = b_h * b_h_2;
        double l_h = std::log((1.0 + b_h) / (k_h + b_h));
        return (1.0 - k_h_3 - 1.5 * b_h * (1.0 - k_h_2) + 3.0 * b_h_2 * (1.0 - k_h) -
                3.0 * b_h_3 * l_h) /
               3.0 / c1_h;
    }

// zero-order approximation for k-m edge solution
    double g_km_0(double s_t) {
        return std::pow((1.0 + beta_m_3 * s_t), 1.0 / 3.0);
    }

// 1st order delta-correction for k-m edge solution
    double g_km_1(double s_t) {
        double a = beta_m_3 * s_t;
        double delta = a / (1.0 + a) / 3.0;
        return std::pow((1.0 + 3.0 * c_1(delta) * s_t), 1.0 / 3.0);
    }

// zero-order approximation for k-m~ edge solution
    double g_kc_0(double chi, double s_t) {
        return std::pow((1.0 + beta_c_4 * chi * s_t), 0.25);
    }

// 1st order delta-correction for k-m~ edge solution
    double g_kc_1(double chi, double s_t) {
        double a = beta_c_4 * chi * s_t;
        double delta = 0.25 * a / (1.0 + a);
        return std::pow((1.0 + 4.0 * chi * c_2(delta) * s_t), 0.25);
    }

// zero-order approximation for universal tip asymptote
    double g_0(double k_h, double c_h) { return effe(k_h, con_mc * c_h, con_m); }

// 1st order delta-correction for universal tip asymptote
    double g_1(double k_h, double c_h) {
        double delta = con_m * (1.0 + con_mc * c_h) * g_0(k_h, c_h);
        return effe(k_h, c_h * b_12(delta), c_1(delta));
    }

////////////////////////////////////////////////////////////////////////////////
// residual functions to find the root
// zero-order
    double res_g_0(double s, TipParameters &taParam) {
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double s_t = s_T(taParam.mu, taParam.k1c, taParam.e_p, taParam.vt, s);
        // use k-m approximate solution for zero leak-off
        if (taParam.cl == 0.0) {
            return g_km_0(s_t) * k_h - 1.0;
        } else {
            double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
            double s_h = s_H(taParam.mu, taParam.e_p, taParam.vt, taParam.wa, s);
            // use universal tip asymptote
            double s_u = g_0(k_h, c_h);
            return s_u - s_h;
        }
    }

////////////////////////////////////////////////////////////////////////////////
    double res_g_0_s(double s, TipParameters &taParam) {
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double s_t = s_T(taParam.mu, taParam.k1c, taParam.e_p, taParam.vt, s);
        // use k-m approximate solution for zero leak-off
        if (taParam.cl == 0.0) {
            return g_km_0(s_t) * k_h - 1.0;
        } else {
            double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
            double s_h = s_H(taParam.mu, taParam.e_p, taParam.vt, taParam.wa, s);
            // use k-m~ approximate solution for large leak-off - small velocity
            if (c_h > 5000.0 * k_h && s_t < 0.001 && k_h > m_tol) {
                double chi = c_h / k_h;
                double g_kc = g_kc_0(chi, s_t);
                double s_u = s_t / g_kc;
                return s_u - s_h;
            } else {
                // use universal tip asymptote
                double s_u = g_0(k_h, c_h);
                return s_u - s_h;
            }
        }
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// 1st order delta-correction
    double res_g_1(double s, TipParameters &taParam) {
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double s_t = s_T(taParam.mu, taParam.k1c, taParam.e_p, taParam.vt, s);
        // use k-m approximate solution for zero leak-off
        if (taParam.cl == 0.0) {
            return g_km_1(s_t) * k_h - 1.0;
        } else {
            double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
            double s_h = s_H(taParam.mu, taParam.e_p, taParam.vt, taParam.wa, s);
            // use universal tip asymptote
            double s_u = g_1(k_h, c_h);
            return s_u - s_h;
        }
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
    double res_g_1_s(double s, TipParameters &taParam) {
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double s_t = s_T(taParam.mu, taParam.k1c, taParam.e_p, taParam.vt, s);
        // use k-m approximate solution for zero leak-off
        if (taParam.cl == 0.0) {
            return g_km_1(s_t) * k_h - 1.0;
        } else {
            double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
            double s_h = s_H(taParam.mu, taParam.e_p, taParam.vt, taParam.wa, s);
            // use k-m~ approximate solution for large leak-off - small velocity
            if (c_h > 5000.0 * k_h && s_t < 0.001 && k_h > m_tol) {
                double chi = c_h / k_h;
                double g_kc = g_kc_1(chi, s_t);
                double s_u = s_t / g_kc;
                return s_u - s_h;
            } else {
                // use universal tip asymptote
                double s_u = g_1(k_h, c_h);
                return s_u - s_h;
            }
        }
    }
////////////////////////////////////////////////////////////////////////////////

// moments (tip volume)
    double deltaP(double k_h, double c_h, double p) {
        double delta = con_m * (1.0 + con_mc * c_h) * g_0(k_h, c_h);
        return (1.0 - p + p * g_0(k_h, c_h)) * delta;
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
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
        double delta_p = deltaP(k_h, c_h, 0.377);
        return (2.0 * taParam.wa * s) / (3.0 + delta_p);
    }

    double moment1(double s, TipParameters &taParam) {
        double k_h = k_H(taParam.k1c, taParam.e_p, taParam.wa, s);
        double c_h = c_H(taParam.cl, taParam.vt, taParam.wa, s);
        double delta_p = deltaP(k_h, c_h, 0.26);
        return (2.0 * taParam.wa * s * s) / (5.0 + delta_p);
    }

// todo: Carter tip leak-off

// checking propagation criterion
    bool isPropagating(TipParameters &taParam) {
        // return (k_H(taParam.k_p, taParam.e_p,
        // taParam.wa, taParam.s_prev) <= 1.0);
        return (taParam.k1c * std::pow(taParam.st, 0.5) <= taParam.e_p * taParam.wa);
    }

// overload for s as an independent parameter
    bool isPropagating(double s, TipParameters &taParam) {
        // return (k_H(taParam.k_p, taParam.e_p, taParam.wa, s) <= 1.0);
        return (taParam.k1c * std::pow(s, 0.5) <= taParam.e_p * taParam.wa);
    }

////////////////////////////////////////////////////////////////////////////////
    il::StaticArray<double, 2> bracket(ResFun resF, TipParameters &taParam, int maxIter) {
        //
        double k_p = 4.0 * std::pow(2.0 / il::pi, 0.5) * taParam.k1c;
        double s_k = std::pow((taParam.e_p * taParam.wa) / k_p, 2);
        il::StaticArray<double, 2> ab;
        // try previous tip position for lower bound
        ab[0] = taParam.st + m_tol;
        // take K-asymptote as upper bound
        ab[1] = s_k;
        double fa = resF(ab[0], taParam), fb = resF(ab[1], taParam);
        if (fa * fb > 0.0) {
            // step back to ~zero
            ab[0] = m_tol;
            fa = resF(ab[0], taParam);
        }
        double mid = ab[1];
        int count = 0;
        while (fa * fb > 0.0 && count < maxIter) {
            // take an intermediate point for lower bound
            mid = (ab[0] + 2.0 * mid) / 3.0;
            ab[0] = mid;
            fa = resF(ab[0], taParam);
            ++count;
        }
        if (count >= maxIter) {
            // if root isn't bracketed
            std::cout << "The root of f(s) can't be bracketed in " << count
                      << " iterations" << std::endl;
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
    double brent(tip::ResFun fun, tip::TipParameters &params, double a0, double b0,
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
// finding the distance to the tip, moments, etc.
    TipParameters propagateTip(ResFun resF, TAInParam &taIn, double epsilonS,
                               int maxIterS, double epsilonV, int maxIterV, bool mute) {
        // Substitute the corresponding residual function
        // resF(x, params)
        // defining the tip asymptote.

        // under-relaxation parameter
        const double beta = 0.5;

        TipParameters taNew;
        taNew.wa = taIn.wa;                          // new ribbon cell opening
        taNew.st = taIn.taPrev.st;                   // previous distance to the tip
        taNew.vt = std::max(taIn.taPrev.vt, m_tol);  // previous tip velocity
        // copying other parameters
        taNew.k1c = taIn.taPrev.k1c;
        taNew.e_p = taIn.taPrev.e_p;
        taNew.cl = taIn.taPrev.cl;
        taNew.mu = taIn.taPrev.mu;
        // initialize volumes
        taNew.m0 = 0.0;
        taNew.m1 = 0.0;
        taNew.cv = 0.0;
        // todo: triggering time
        taNew.ta = taIn.taPrev.ta;

        // check K vs K1c (taParam.k_p)
        bool isProp = isPropagating(taNew);

        if (!isProp) {
            // not propagating
            taNew.vt = 0.0;
            taNew.st = taIn.taPrev.st;
        } else {
            double s_0, s_a = taIn.taPrev.st + m_tol, v_a;
            int iter_count = 0;
            // iterating for tip velocity
            do {
                s_0 = s_a;
                // bracketing the tip position for fine root finding
                il::StaticArray<double, 2> ab = bracket(resF, taNew, 30);
                // fine root finding (adjusting tip position)
                s_a = brent(resF, taNew, ab[0], ab[1], epsilonS, maxIterS);
                //                s_a = root::brent(resF, taNew, ab[0], ab[1],
                //                                  epsilonS, maxIterS);
                // adjusting tip velocity
                v_a = (s_a - taIn.taPrev.st) / taIn.dt;
                if (v_a < m_tol) {
                    taNew.vt = m_tol;
                    // taNew.vt = (1.0 - beta) * taNew.vt;
                    s_a = taIn.taPrev.st + taNew.vt * taIn.dt;
                } else {
                    // under-relaxation
                    taNew.vt = beta * v_a + (1.0 - beta) * taNew.vt;
                };
                // s_a = taIn.taPrev.st + taNew.vt * taIn.dt;
                // k1c and e_p can be updated here in case of heterogeneity
                ++iter_count;
                if (s_a < 0) {
                    // other errors
                    std::cout << "Iteration " << iter_count << "; Failed to converge on s"
                              << std::endl;
                } else {
                    if (!mute) {
                        std::cout << "Iteration " << iter_count << "; s=" << s_a
                                  << "; v=" << taNew.vt << std::endl;
                    }
                }
            } while (std::abs(s_a - s_0) > epsilonV * s_0 && iter_count < maxIterV &&
                     s_a >= 0);
            if (iter_count >= maxIterV) {
                // if root isn't bracketed
                std::cout << "Failed to converge on v" << std::endl;
            }
            // new tip position
            taNew.st = s_a;
        }
        // moments
        taNew.m0 = moment0(taNew);
        taNew.m1 = moment1(taNew);
        // todo: Carter tip leak-off
        return taNew;
    }
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
    TipParameters propagateTip(ResFun resF, TAInParam &taIn, double epsilon, int maxIter,
                               bool mute) {
        return propagateTip(resF, taIn, epsilon, maxIter, epsilon, maxIter, mute);
    }
////////////////////////////////////////////////////////////////////////////////
}
