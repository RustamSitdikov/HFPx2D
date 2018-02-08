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

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

// Inclusion from the project
#include <src/core_dev/SolidEvolution.h>

namespace hfp2d {

class FractureEvolution {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
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

  il::Array<double> getInitPermeab() { return initial_permeabilitys_; }
  double getInitPermeab(il::int_t i) {
    return initial_permeabilitys_[i];
  };

  il::Array<double> getIncrPermeab() { return increment_permeabilitys_; }
  double getIncrPermeab(il::int_t i) {
    return increment_permeabilitys_[i];
  };

  il::Array<double> getResidSlip() { return residual_slips_; }
  double getResidSlip(il::int_t i) { return residual_slips_[i]; };

  il::Array<double> getInitHydrWidth() {
    return initial_hydraulic_widths_;}
  double getInitHydrWidth(il::int_t i) {
    return initial_hydraulic_widths_[i]; };

  il::Array<double> getIncrHydrWidth() {
    return increment_hydraulic_widths_; }
  double getIncrHydrWidth(il::int_t i) {
    return increment_hydraulic_widths_[i];
  };

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  il::Array<double> linearDilatancy(const il::Array<double> &slip,
                                    FractureEvolution &FractureEvolution) {

    il::Array<double> new_hydr_width{slip.size(), 0.};

    IL_EXPECT_FAST(new_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.initial_permeabilitys_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_permeabilitys_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.initial_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_hydr_width.size(); ++i) {
      if (slip[i] < FractureEvolution.residual_slips_[i]) {
        new_hydr_width[i] =
            FractureEvolution.initial_hydraulic_widths_[i] +
            (slip[i] * (FractureEvolution.increment_hydraulic_widths_[i] /
                        FractureEvolution.residual_slips_[i]));

      } else {
        new_hydr_width[i] = FractureEvolution.initial_hydraulic_widths_[i] +
                            FractureEvolution.increment_hydraulic_widths_[i];
      }
    }

    return new_hydr_width;
  };

  il::Array<double> derivativeLinearDilatancy(
     const il::Array<double> &slip, FractureEvolution &FractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_deriv_hydr_width.size(); ++i) {
      if (slip[i] < FractureEvolution.residual_slips_[i]) {
        new_deriv_hydr_width[i] =
            FractureEvolution.increment_hydraulic_widths_[i] /
            FractureEvolution.residual_slips_[i];

      } else {
        new_deriv_hydr_width[i] = 0;
      }
    }

    return new_deriv_hydr_width;
  };

  il::Array<double> exponentialDilatancy(const il::Array<double> &slip,
                                         FractureEvolution &FractureEvolution) {
    il::Array<double> new_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_hydr_width.size(); ++i) {
      new_hydr_width[i] =
          FractureEvolution.increment_hydraulic_widths_[i] +
          (FractureEvolution.increment_hydraulic_widths_[i] *
           (1 - exp(-slip[i] / FractureEvolution.residual_slips_[i])));
    }

    return new_hydr_width;
  };

  il::Array<double> derivativeExponentialDilatancy(
      const il::Array<double> &slip, FractureEvolution &FractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_deriv_hydr_width.size(); ++i) {
      new_deriv_hydr_width[i] =
          (FractureEvolution.increment_hydraulic_widths_[i] /
           FractureEvolution.residual_slips_[i]) *
          (exp(-slip[i] / FractureEvolution.residual_slips_[i]));
    }

    return new_deriv_hydr_width;
  };

  il::Array<double> linearPermeability(const il::Array<double> &slip,
                                       FractureEvolution &FractureEvolution) {
    il::Array<double> new_permeab{slip.size(), 0};

    IL_EXPECT_FAST(new_permeab.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_permeab.size(); ++i) {
      if (slip[i] < FractureEvolution.getResidSlip(i)) {
        new_permeab[i] = FractureEvolution.getInitPermeab(i) +
                         (slip[i] * (FractureEvolution.getIncrPermeab(i) /
                                     FractureEvolution.getResidSlip(i)));

      } else {
        new_permeab[i] = FractureEvolution.getInitPermeab(i) +
                         FractureEvolution.getIncrPermeab(i);
      }
    }

    return new_permeab;
  };


  il::Array<double> exponentialPermeability(
      const il::Array<double> &slip, FractureEvolution &FractureEvolution) {
    il::Array<double> new_permeab{slip.size(), 0};

    IL_EXPECT_FAST(new_permeab.size() == slip.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() == slip.size());

    for (int i = 0; i < new_permeab.size(); ++i) {
      new_permeab[i] =
          FractureEvolution.getInitPermeab(i) +
          (FractureEvolution.getIncrPermeab(i) *
           (1 - exp(-slip[i] / FractureEvolution.getResidSlip(i))));
    }

    return new_permeab;
  };

};
}

#endif  // HFPX2D_FRACTUREEVOLUTION_H
