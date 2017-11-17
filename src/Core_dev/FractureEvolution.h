//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 10.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FRACTUREEVOLUTION_H
#define HFPX2D_FRACTUREEVOLUTION_H

// Inclusion from standard library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>
#include <src/Core_dev/SolidEvolution.h>

namespace hfp2d {

class FractureEvolution {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
  il::String type_;
  il::Array<double> initial_permeabilitys_;
  il::Array<double> increment_permeabilitys_;
  il::Array<double> residual_slips_;
  il::Array<double> initial_hydraulic_widths_;
  il::Array<double> increment_hydraulic_widths_;

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////

 public:
  FractureEvolution() = default;

  FractureEvolution(il::Array<double> &vector_init_permeab,
                    il::Array<double> &vector_incr_permeab,
                    il::Array<double> &vector_resid_slip,
                    il::Array<double> &vector_init_hydr_width,
                    il::Array<double> &vector_incr_hydr_width) {
    initial_permeabilitys_ = vector_init_permeab;
    increment_permeabilitys_ = vector_incr_permeab;
    residual_slips_ = vector_resid_slip;
    initial_hydraulic_widths_ = vector_init_hydr_width;
    increment_hydraulic_widths_ = vector_incr_hydr_width;
  }

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  il::Array<double> getInitPermab() { return initial_permeabilitys_; }
  il::Array<double> getIncrPermeab() { return increment_permeabilitys_; }
  il::Array<double> getResidSlip() { return residual_slips_; }
  il::Array<double> getInitHydrWidth() { return initial_hydraulic_widths_; }
  il::Array<double> getIncrHydrWidth() { return increment_hydraulic_widths_; }
  il::String getType() { return type_; }

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  void setFractureEvolution(il::Array<double> &init_permeab,
                            il::Array<double> &incr_permeab,
                            il::Array<double> &res_slip,
                            il::Array<double> &init_hydr_width,
                            il::Array<double> &incr_hydr_width) {
    initial_permeabilitys_ = init_permeab;
    increment_permeabilitys_ = incr_permeab;
    residual_slips_ = res_slip;
    initial_hydraulic_widths_ = init_hydr_width;
    increment_hydraulic_widths_ = incr_hydr_width;
  };

  il::Array<double> linearDilatancy(il::Array<double> &slip,
                                    FractureEvolution &InitFractureEvolution) {
    il::Array<double> new_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.initial_permeabilitys_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.increment_permeabilitys_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.initial_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_hydr_width.size(); ++i) {
      if (slip[i] < InitFractureEvolution.residual_slips_[i]) {
        new_hydr_width[i] =
            InitFractureEvolution.initial_hydraulic_widths_[i] +
            (slip[i] * (InitFractureEvolution.increment_hydraulic_widths_[i] /
                        InitFractureEvolution.residual_slips_[i]));

      } else {
        new_hydr_width[i] =
            InitFractureEvolution.initial_hydraulic_widths_[i] +
            InitFractureEvolution.increment_hydraulic_widths_[i];
      }
    }

    return new_hydr_width;
  };

  il::Array<double> derivativeLinearDilatancy(
      il::Array<double> &slip, FractureEvolution &InitFractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_deriv_hydr_width.size(); ++i) {
      if (slip[i] < InitFractureEvolution.residual_slips_[i]) {
        new_deriv_hydr_width[i] =
            InitFractureEvolution.increment_hydraulic_widths_[i] /
            InitFractureEvolution.residual_slips_[i];

      } else {
        new_deriv_hydr_width[i] = 0;
      }
    }

    return new_deriv_hydr_width;
  };

  il::Array<double> exponentialDilatancy(
      il::Array<double> &slip, FractureEvolution &InitFractureEvolution) {
    il::Array<double> new_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_hydr_width.size(); ++i) {
      new_hydr_width[i] =
          InitFractureEvolution.increment_hydraulic_widths_[i] +
          (InitFractureEvolution.increment_hydraulic_widths_[i] *
           (1 - exp(-slip[i] / InitFractureEvolution.residual_slips_[i])));
    }

    return new_hydr_width;
  };

  il::Array<double> derivativeExponentialDilatancy(
      il::Array<double> &slip, FractureEvolution &InitFractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(InitFractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_deriv_hydr_width.size(); ++i) {
      new_deriv_hydr_width[i] =
          (InitFractureEvolution.increment_hydraulic_widths_[i] /
           InitFractureEvolution.residual_slips_[i]) *
          (exp(-slip[i] / InitFractureEvolution.residual_slips_[i]));
    }

    return new_deriv_hydr_width;
  };
};
}

#endif  // HFPX2D_FRACTUREEVOLUTION_H
