//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 03.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <il/norm.h>
#include <src/wellbore/WellInjection.h>


namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
//        Injection class methods
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// adjusting the outfluxes (volume rates) from the wellbore into the HFs
// according to the pressure drop(s) at the HF locations
// see Lecampion & Desroches 2015
// Technically, it's a fixed-point iteration step to be used
// as a feedback between the HF and the wellbore solvers
//void WellInjection::adjHFVolRate() {
//  il::int_t n_hf = hf_rates_.size();
//  il::Array<double> new_vol_rate{n_hf};
//  for (il::int_t i = 0; i < n_hf; ++i) {
//    new_vol_rate[i] =
//        hf_p_drop_[i] /
//        (coef_perf_[i] * hf_rates_[i] +
//         coef_tort_[i] * std::pow(hf_rates_[i], beta_tort_[i] - 1.0));
//  }
//  setHFVolRate(new_vol_rate);
//}
//
//
//// overload - calculate the residue norm (error) ????
//void WellInjection::adjHFVolRate(il::io_t, double &err) {
//  il::int_t n_hf = hf_rates_.size();
//  il::Array<double> new_vol_rate{n_hf};
//  for (il::int_t i = 0; i < n_hf; ++i) {
//    new_vol_rate[i] =
//        hf_p_drop_[i] /
//        (coef_perf_[i] * hf_rates_[i] +
//         coef_tort_[i] * std::pow(hf_rates_[i], beta_tort_[i] - 1.0));
//  }
//  setHFVolRate(new_vol_rate, il::io, err);
//}
//
//// overload - calculate but not update - is this needed ?
//il::Array<double> WellInjection::adjHFVolRate(WellInjection &well_inj) {
//  il::Array<double> &hf_vol_rate = well_inj.hf_rates_;
//  il::Array<double> &hf_p_drop = well_inj.hf_p_drop_;
//  il::int_t n_hf = hf_vol_rate.size();
//  il::Array<double> new_vol_rate{n_hf};
//  for (il::int_t i = 0; i < n_hf; ++i) {
//    new_vol_rate[i] =
//        hf_p_drop[i] /
//        (coef_perf_[i] * hf_vol_rate[i] +
//         coef_tort_[i] * std::pow(hf_vol_rate[i], beta_tort_[i] - 1.0));
//  }
//  well_inj.hf_rates_ = new_vol_rate;
//  return new_vol_rate;
//}

// setting volume rate(s) into the HF(s)
void WellInjection::setHFVolRate(il::Array<double> &new_vol_rate) {
  hf_rates_ = new_vol_rate;
}

// overload - calculate the residue norm (error)
void WellInjection::setHFVolRate(il::Array<double> &new_vol_rate, il::io_t,
                                 double &err) {
  il::Array<double> diff{new_vol_rate.size()};
  for (il::int_t i = 0; i < new_vol_rate.size(); ++i) {
    if (!isnan(hf_rates_[i])) {
      diff[i] = new_vol_rate[i] - hf_rates_[i];
    } else {
      diff[i] = new_vol_rate[i];
    }
  }
  err = il::norm(diff, il::Norm::L2);
  hf_rates_ = new_vol_rate;
}


}
