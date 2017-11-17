//
// Created by Federico Ciardo on 10/11/17.
//
//
// This file is part of HFPx2DUnitTest.
//
// Created by lorenzo on 9/18/17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2DUNITTEST_SOLIDEVOLUTION_H
#define HFPX2DUNITTEST_SOLIDEVOLUTION_H

// Inclusion from standard library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>
#include <il/base.h>

namespace hfp2d {

class SolidEvolution {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
  il::String type_;
  il::Array<double> friction_coefficients_;
  il::Array<double> peak_friction_coefficients_;
  il::Array<double> residual_friction_coefficients_;
  il::Array<double> residual_slips_;
  il::Array<double> failure_stresses_;
  il::Array<double> maximum_openings_;

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////
 public:
  SolidEvolution() = default;

  SolidEvolution(il::Array<double> &vector_peak_fric_coeff,
                 il::Array<double> &vector_residual_fric_coeff,
                 il::Array<double> &vector_residual_slip) {
    friction_coefficients_ = vector_peak_fric_coeff;
    peak_friction_coefficients_ = vector_peak_fric_coeff;
    residual_friction_coefficients_ = vector_residual_fric_coeff;
    residual_slips_ = vector_residual_slip;
  };

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  inline il::Array<double> getFricCoeff() { return friction_coefficients_; }

  inline double getFricCoeff(il::int_t i) { return friction_coefficients_[i]; }

  inline il::Array<double> getPeakFricCoeff() {
    return peak_friction_coefficients_;
  }
  inline double getPeakFricCoeff(il::int_t i) {
    return peak_friction_coefficients_[i];
  }
  inline il::Array<double> getResidFricCoeff() {
    return residual_friction_coefficients_;
  }
  inline double getResidFricCoeff(il::int_t i) {
    return residual_friction_coefficients_[i];
  }
  inline il::Array<double> getResidSlip() { return residual_slips_; }
  inline double getResidSlip(il::int_t i) { return residual_slips_[i]; }
  inline il::String getType() { return type_; }

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  void setSolidEvolution(il::Array<double> &current_fric_coeff,
                         il::Array<double> &peak_fric_coeff,
                         il::Array<double> &res_fric_coeff,
                         il::Array<double> &res_slip,
                         il::Array<double> &fail_stress,
                         il::Array<double> &max_openings) {
    friction_coefficients_ = current_fric_coeff;
    peak_friction_coefficients_ = peak_fric_coeff;
    residual_friction_coefficients_ = res_fric_coeff;
    residual_slips_ = res_slip;
    failure_stresses_ = fail_stress;
    maximum_openings_ = max_openings;
  }

  il::Array<double> linearFricWeakLaw(
      il::Array<double> &slip, const SolidEvolution &InitSolidEvolution) {
    IL_EXPECT_FAST(InitSolidEvolution.peak_friction_coefficients_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitSolidEvolution.residual_friction_coefficients_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitSolidEvolution.residual_slips_.size() == slip.size());

    for (il::int_t i = 0; i < slip.size(); ++i) {
      if (slip[i] < InitSolidEvolution.residual_slips_[i]) {
        friction_coefficients_[i] =
            InitSolidEvolution.peak_friction_coefficients_[i] -
            (((InitSolidEvolution.peak_friction_coefficients_[i] -
               InitSolidEvolution.residual_friction_coefficients_[i]) /
              InitSolidEvolution.residual_slips_[i]) *
             slip[i]);
      } else {
        friction_coefficients_[i] =
            InitSolidEvolution.residual_friction_coefficients_[i];
      }
    }
    return friction_coefficients_;
  };

  il::Array<double> exponentialFricWeakLaw(
      il::Array<double> &slip, const SolidEvolution &InitSolidEvolution) {
    IL_EXPECT_FAST(InitSolidEvolution.peak_friction_coefficients_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitSolidEvolution.residual_friction_coefficients_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitSolidEvolution.residual_slips_.size() == slip.size());

    for (il::int_t i = 0; i < slip.size(); ++i) {
      friction_coefficients_[i] =
          InitSolidEvolution.peak_friction_coefficients_[i] -
          ((InitSolidEvolution.peak_friction_coefficients_[i] -
            InitSolidEvolution.residual_friction_coefficients_[i]) *
           (1 - exp(-slip[i] / InitSolidEvolution.residual_slips_[i])));
    }
    return friction_coefficients_;
  }
};
}

#endif  // HFPX2DUNITTEST_SOLIDEVOLUTION_H
