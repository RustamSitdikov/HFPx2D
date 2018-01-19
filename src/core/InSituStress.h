//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 13.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_INSITUSTRESS_H
#define HFPX2D_INSITUSTRESS_H

namespace hfp2d {
class InSituStress {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////
 private:
  il::Array2D<double> local_insitu_stresses_;

  il::Array<double> ambient_pore_pressure_;

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////

 public:
  InSituStress() = default;

  InSituStress(il::Array2D<double> &local_in_situ_stress_state,
               il::Array<double> &ambient_pore_pressure) {
    local_insitu_stresses_ = local_in_situ_stress_state;
    ambient_pore_pressure_ = ambient_pore_pressure;
  }

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  il::Array2D<double> getLocalInSituStresses() {
    return local_insitu_stresses_;
  }

  double getLocalInSituStresses(il::int_t i, il::int_t j) {
    return local_insitu_stresses_(i, j);
  }

  il::Array<double> getAmbientPorePressure() { return ambient_pore_pressure_; }

  double getAmbientPorePressure(il::int_t i) {
    return ambient_pore_pressure_[i];
  }

  il::Array<double> getBackgroundShearStress() {
    il::Array<double> tau_0{this->local_insitu_stresses_.size(0)};

    for (int i = 0; i < tau_0.size(); ++i) {
      tau_0[i] = this->local_insitu_stresses_(i, 0);
    }

    return tau_0;
  };

  double getBackgroundShearStress(il::int_t i) {
    return this->local_insitu_stresses_(i, 0);
  };

  il::Array<double> getBackgroundNormalStress() {
    il::Array<double> sigma_n{this->local_insitu_stresses_.size(0)};

    for (int i = 0; i < sigma_n.size(); ++i) {
      sigma_n[i] = this->local_insitu_stresses_(i, 1);
    }

    return sigma_n;
  };
  double getBackgroundNormalStress(il::int_t i) {
    return this->local_insitu_stresses_(i, 1);
  };
};
}

#endif  // HFPX2D_INSITUSTRESS_H
