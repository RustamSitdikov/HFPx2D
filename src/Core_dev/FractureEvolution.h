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

  inline il::Array<double> getInitPermeab() { return initial_permeabilitys_; }
  inline double getInitPermeab(il::int_t i) {
    return initial_permeabilitys_[i];
  };

  inline il::Array<double> getIncrPermeab() { return increment_permeabilitys_; }

  inline double getIncrPermeab(il::int_t i) {
    return increment_permeabilitys_[i];
  };
  inline il::Array<double> getResidSlip() { return residual_slips_; }
  inline double getResidSlip(il::int_t i) { return residual_slips_[i]; };
  inline il::Array<double> getInitHydrWidth() {
    return initial_hydraulic_widths_;
  }
  inline double getInitHydrWidth(il::int_t i) {
    return initial_hydraulic_widths_[i];
  };
  inline il::Array<double> getIncrHydrWidth() {
    return increment_hydraulic_widths_;
  }
  inline double getIncrHydrWidth(il::int_t i) {
    return increment_hydraulic_widths_[i];
  };
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
                                    FractureEvolution &FractureEvolution) {
    il::Array<double> new_hydr_width{slip.size(), 0};

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

  il::Array<double> linearDilatancyMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_hydr_width{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_hydr_width.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.initial_permeabilitys_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_permeabilitys_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.initial_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() / 2 ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_hydr_width.size(); ++i, k = k + 2) {
      if (slip_middle[i] < ((FractureEvolution.residual_slips_[k] +
                             FractureEvolution.residual_slips_[k + 1]) /
                            2)) {
        new_hydr_width[i] =
            ((FractureEvolution.initial_hydraulic_widths_[k] +
              FractureEvolution.initial_hydraulic_widths_[k + 1]) /
             2) +
            (slip_middle[i] *
             (((FractureEvolution.increment_hydraulic_widths_[k] +
                FractureEvolution.increment_hydraulic_widths_[k + 1]) /
               2) /
              ((FractureEvolution.residual_slips_[k] +
                FractureEvolution.residual_slips_[k + 1]) /
               2)));

      } else {
        new_hydr_width[i] =
            ((FractureEvolution.initial_hydraulic_widths_[k] +
              FractureEvolution.initial_hydraulic_widths_[k + 1]) /
             2) +
            ((FractureEvolution.increment_hydraulic_widths_[k] +
              FractureEvolution.increment_hydraulic_widths_[k + 1]) /
             2);
      }
    }

    return new_hydr_width;
  };

  il::Array<double> derivativeLinearDilatancy(
      il::Array<double> &slip, FractureEvolution &FractureEvolution) {
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

  il::Array<double> derivativeLinearDilatancyMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() / 2 ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_deriv_hydr_width.size(); ++i, k = k + 2) {
      if (slip_middle[i] < ((FractureEvolution.residual_slips_[k] +
                             FractureEvolution.residual_slips_[k + 1]) /
                            2)) {
        new_deriv_hydr_width[i] =
            ((FractureEvolution.increment_hydraulic_widths_[k] +
              FractureEvolution.increment_hydraulic_widths_[k + 1]) /
             2) /
            ((FractureEvolution.residual_slips_[k] +
              FractureEvolution.residual_slips_[k + 1]) /
             2);

      } else {
        new_deriv_hydr_width[i] = 0;
      }
    }

    return new_deriv_hydr_width;
  };

  il::Array<double> exponentialDilatancy(il::Array<double> &slip,
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

  il::Array<double> exponentialDilatancyMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_hydr_width{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_hydr_width.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_hydr_width.size(); ++i, k = k + 2) {
      new_hydr_width[i] =
          ((FractureEvolution.increment_hydraulic_widths_[k] +
            FractureEvolution.increment_hydraulic_widths_[k + 1]) /
           2) +
          (((FractureEvolution.increment_hydraulic_widths_[k] +
             FractureEvolution.increment_hydraulic_widths_[k + 1]) /
            2) *
           (1 -
            exp(-slip_middle[i] / ((FractureEvolution.residual_slips_[k] +
                                    FractureEvolution.residual_slips_[k + 1]) /
                                   2))));
    }
    return new_hydr_width;
  };

  il::Array<double> derivativeExponentialDilatancy(
      il::Array<double> &slip, FractureEvolution &FractureEvolution) {
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

  il::Array<double> derivativeExponentialDilatancyMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_deriv_hydr_width{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_deriv_hydr_width.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() / 2 ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_deriv_hydr_width.size(); ++i, k = k + 2) {
      new_deriv_hydr_width[i] =
          (((FractureEvolution.increment_hydraulic_widths_[k] +
             FractureEvolution.increment_hydraulic_widths_[k + 1]) /
            2) /
           ((FractureEvolution.residual_slips_[k] +
             FractureEvolution.residual_slips_[k + 1]) /
            2)) *
          (exp(-slip_middle[i] / ((FractureEvolution.residual_slips_[k] +
                                   FractureEvolution.residual_slips_[k + 1]) /
                                  2)));
    }

    return new_deriv_hydr_width;
  };

  il::Array<double> linearPermeability(il::Array<double> &slip,
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

  il::Array<double> linearPermeabilityMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_permeab{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_permeab.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() / 2 ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_permeab.size(); ++i, k = k + 2) {
      if (slip_middle[i] < ((FractureEvolution.getResidSlip(k) +
                             FractureEvolution.getResidSlip(k + 1)) /
                            2)) {
        new_permeab[i] =
            ((FractureEvolution.getInitPermeab(k) +
              FractureEvolution.getInitPermeab(k + 1)) /
             2) +
            (slip_middle[i] * (((FractureEvolution.getIncrPermeab(k) +
                                 FractureEvolution.getIncrPermeab(k + 1)) /
                                2) /
                               ((FractureEvolution.getResidSlip(k) +
                                 FractureEvolution.getResidSlip(k + 1)) /
                                2)));

      } else {
        new_permeab[i] = ((FractureEvolution.getInitPermeab(k) +
                           FractureEvolution.getInitPermeab(k + 1)) /
                          2) +
                         ((FractureEvolution.getIncrPermeab(k) +
                           FractureEvolution.getIncrPermeab(k + 1)) /
                          2);
      }
    }
    return new_permeab;
  };

  il::Array<double> exponentialPermeability(
      il::Array<double> &slip, FractureEvolution &FractureEvolution) {
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

  il::Array<double> exponentialPermeabilityMiddle(
      il::Array<double> &slip_middle, FractureEvolution &FractureEvolution) {
    il::Array<double> new_permeab{slip_middle.size(), 0};

    IL_EXPECT_FAST(new_permeab.size() == slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.increment_hydraulic_widths_.size() / 2 ==
                   slip_middle.size());
    IL_EXPECT_FAST(FractureEvolution.residual_slips_.size() / 2 ==
                   slip_middle.size());

    for (int i = 0, k = 0; i < new_permeab.size(); ++i, k = k + 2) {
      new_permeab[i] =
          ((FractureEvolution.getInitPermeab(k) +
            FractureEvolution.getInitPermeab(k + 1)) /
           2) +
          (((FractureEvolution.getIncrPermeab(k) +
             FractureEvolution.getIncrPermeab(k + 1)) /
            2) *
           (1 - exp(-slip_middle[i] / ((FractureEvolution.getResidSlip(k) +
                                        FractureEvolution.getResidSlip(k + 1)) /
                                       2))));
    }

    return new_permeab;
  };
};
}

#endif  // HFPX2D_FRACTUREEVOLUTION_H
