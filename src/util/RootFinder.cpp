//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 18.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <algorithm>
#include <cmath>
#include <iostream>

#include <src/util/RootFinder.h>




namespace imf {

// Brent's root finding method
//
// Brent R.P. (1973),
// "Chapter 4: An Algorithm with Guaranteed Convergence
// for Finding a Zero of a Function",
// Algorithms for Minimization without Derivatives,
// Englewood Cliffs, NJ: Prentice-Hall, ISBN 0-13-022335-2
//
double brent(ImplicitFunction fun, IFParameters &params, double a0, double b0,
             const double epsilon, const int max_iter, bool mute) {
  //
  double a = a0;
  double b = b0;
  double fa = fun(a, params);  // calculated now to save function calls
  double fb = fun(b, params);  // calculated now to save function calls
  double fs = 0;               // initialize

  if (std::isnan(fa)) {
    if (!mute) {
      std::cout << "f(a) isn't finite; a=" << a << std::endl;
    }
    return -111;
  }

  if (std::isnan(fb)) {
    if (!mute) {
      std::cout << "f(b) isn't finite; b=" << b << std::endl;
    }
    return -111;
  }

  if (fa * fb > 0.0) {
    // throw exception if root isn't bracketed
    if (!mute) {
      std::cout << "Signs of f(a) and f(b) must be opposite" << std::endl;
    }
    return -11;
  }

  // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  if (std::abs(fa) < std::abs(b)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  double c =
      a;  // c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;  // precompute function evalutation for point c by assigning
  // it the same value as fa
  bool mflag = true;  // boolean flag used to evaluate if statement later on
  double s = 0;       // Our Root that will be returned
  double d = 0;       // Only used if mflag is unset (mflag == false)

  for (unsigned int iter = 1; iter < max_iter; ++iter) {
    // stop if converged on root or error is less than tolerance
    if (std::abs(b - a) < epsilon) {
      if (!mute) {
        std::cout << "After " << iter << " iterations the root is: " << s
                  << std::endl;
      }
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
  if (!mute) {
    std::cout << "The solution did not converge after " << max_iter
              << " iterations" << std::endl;
  }

  return -1;
}

}
