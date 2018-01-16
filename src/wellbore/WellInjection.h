//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 03.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_WELLINJECTION_H
#define HFPX2D_WELLINJECTION_H

#include <il/math.h>
#include <il/Array.h>
#include <il/Array2D.h>

#include <src/core/Sources.h>

namespace hfp2d {

//////////////////////////////////////////////////////////////////////////
// injection class / outflow to HFs
//////////////////////////////////////////////////////////////////////////

class WellInjection {

 private:
  // injection rate at the well source (1st element for the time being)
  double well_inj_rate_;
  // we are far away of injecting simultaneously in 2 boreholes!
  // then when rate is variable  will have input table of rates and time etc.

  // outflux (HFs) locations on the wellbore
  // see "source_location_" for the corresponding
  // locations on the HF wellMesh
//  il::Array<il::int_t> hf_location_;
  // rate of outflow in hydraulic fractures.
//   il::Array<double> hf_rates_;

  // nonlinear HF entry friction model parameters
  // (see Lecampion & Desroches 2015)
  // perforation resistance (delta_p = coef_perf_ * Q^2)
  il::Array<double> coef_perf_;
  // tortuosity resistance (delta_p = coef_tort_ * Q^beta_tort_)
  il::Array<double> coef_tort_;
  il::Array<double> beta_tort_;

 public:
  ////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS
  ////////////////////////////////////////////////////////////////////////

  // default constructor
  WellInjection(){};

  // normal constructor
  WellInjection(double well_inj_rate,
                //il::Array<il::int_t> &hf_location,
                //il::int_t plug_location,
               // il::Array<double> &hf_vol_rate,
                //il::Array<double> &hf_p_drop,
                il::Array<double> &coef_perf,
                il::Array<double> &coef_tort, il::Array<double> &beta_tort) {
    well_inj_rate_ = well_inj_rate;
    //plug_location_ = plug_location;

 //   IL_EXPECT_FAST(hf_location.size()==hf_vol_rate.size());
    IL_EXPECT_FAST(coef_perf.size()==coef_tort.size());
    IL_EXPECT_FAST(coef_tort.size()==beta_tort.size());
    IL_EXPECT_FAST(beta_tort.size()==coef_perf.size());

 //   hf_location_ = hf_location;
  //  hf_rates_ = hf_vol_rate;
   // hf_p_drop_ = hf_p_drop;

    coef_perf_ = coef_perf;
    coef_tort_ = coef_tort;
    beta_tort_ = beta_tort;
  };

  ////////////////////////////////////////////////////////////////////////
  //        public interfaces
  ////////////////////////////////////////////////////////////////////////

  double wellInjRate() { return well_inj_rate_; }
  //il::int_t numberOfHFs() { return hf_location_.size(); }

//  il::Array<il::int_t> hfLocation() { return hf_location_; }
//  il::int_t hfLocation(il::int_t nf) { return hf_location_[nf]; }
//
//  il::Array<double> hfRate() { return hf_rates_; }
//  double hfRate(il::int_t nf) { return hf_rates_[nf]; }

 // inline il::Array<double> hfPressDrop() { return hf_p_drop_; }
 // inline double hfPressDrop(il::int_t nf) { return hf_p_drop_[nf]; }

  //inline il::int_t plugLocation() { return plug_location_; }

  il::Array<double> coefPerf() { return coef_perf_; }
  double coefPerf(il::int_t nf) { return coef_perf_[nf]; }
  il::Array<double> coefTort() { return coef_tort_; }

  double coefTort(il::int_t nf) { return coef_tort_[nf]; }
  il::Array<double> betaTort() { return beta_tort_; }
  double betaTort(il::int_t nf) { return beta_tort_[nf]; }

  ////////////////////////////////////////////////////////////////////////
  //        set functions
  ////////////////////////////////////////////////////////////////////////

  // setting the well injection rate
  void setWellInjRate(double well_inj_rate) { well_inj_rate_ = well_inj_rate; };

  // setting volume rate(s) into the HF(s)
 // void setHFVolRate(il::Array<double> &newVolRate);

  // return the residue norm (error)
//  void setHFVolRate(il::Array<double> &newVolRate, il::io_t, double &err);

  ////////////////////////////////////////////////////////////////////////
  //        methods
  ////////////////////////////////////////////////////////////////////////

  // adjustment of volume rate(s) into the HF(s)
  // technically, a fixed-point iteration step
//  void adjHFVolRate();
//  // overload - return the residue norm (error)
//  void adjHFVolRate(il::io_t, double &err);
//  // overload - calculate but not update
//  il::Array<double> adjHFVolRate(WellInjection &well_inj);
};
};

#endif  // HFPX2D_WELLINJECTION_H
