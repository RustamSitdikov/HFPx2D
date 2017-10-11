//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 07.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_SIMULATIONPARAMETERS_H
#define HFPX2D_SIMULATIONPARAMETERS_H

namespace hfp2d {

class SimulationParameters {
  // placeholder for simulation parameters

 private:
  il::int_t ehl_max_its_ =
      100;  // max number of iterations for elastic - fluid coupled problen
  il::int_t fracfront_max_its_ =
      30;  // max number of its for fracture front loop.

  double ehl_betarela_ =
      0.5;  //  ehl under-relaxation parameters. must be in [0,1]
  double ehl_tolerance_ = 1.e-6;  // tolerance for stopping criteria for elastic
                                  // - fluid coupled problem
  double fracfront_tolerance_ =
      1.e-3;  // tolerance for stopping criteria on fracture front position

  // a string for an option to adapt the under-relaxation

 public:
  // constructor:
  SimulationParameters(){};

  SimulationParameters(il::int_t ehl_max_its) { ehl_max_its_ = ehl_max_its; }
  SimulationParameters(const il::int_t ehl_max_its,
                       const il::int_t fracfront_max_its) {
    ehl_max_its_ = ehl_max_its;
    fracfront_max_its_ = fracfront_max_its;
  }
  SimulationParameters(const il::int_t ehl_max_its,
                       const il::int_t fracfront_max_its, const double beta) {
    ehl_max_its_ = ehl_max_its;
    fracfront_max_its_ = fracfront_max_its;
    if (beta > 0 && beta <= 1.) {
      ehl_betarela_ = beta;
    }
  }

  SimulationParameters(const il::int_t ehl_max_its,
                       const il::int_t fracfront_max_its, const double beta,
                       const double ehl_tolerance,
                       const double fracfront_tolerance) {
    ehl_max_its_ = ehl_max_its;
    fracfront_max_its_ = fracfront_max_its;
    if (beta > 0 && beta <= 1.) {
      ehl_betarela_ = beta;
    };
    ehl_tolerance_ =
        ehl_tolerance;  // todo : put some checks here to avoid some weird stuff
    fracfront_tolerance_ = fracfront_tolerance;  // todo : put some checks here
                                                 // to avoid some weird stuff
  }

  // Get functions
  il::int_t EHL_max_its() const { return ehl_max_its_; };
  il::int_t Front_max_its() const { return fracfront_max_its_; };

  double EHL_tol() const { return ehl_tolerance_; };
  double Front_tol() const { return fracfront_tolerance_; };
  double EHL_relaxation() const { return ehl_betarela_; };
};
}
#endif  // HFPX2_SIMULATIONPARAMETERS_H
