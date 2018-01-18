//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 03.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_WELLSOLUTION_H
#define HFPX2D_WELLSOLUTION_H

#include <il/Array.h>
#include <il/Array2D.h>

#include <src/wellbore/WellInjection.h>
#include <src/wellbore/WellMesh.h>
#include <src/util/json.hpp>

namespace hfp2d {

//////////////////////////////////////////////////////////////////////////
// solution for the wellbore flow
// for the time being, a staggered scheme (w. HF EHD solver) is planned
//////////////////////////////////////////////////////////////////////////
class WellSolution {
 private:
  // current time (tn+1 = tn + timestep)
  double time_;
  // time step
  double timestep_;

  // fluid pressure (at cell center)
  il::Array<double> pressure_;
  // hydrostatic fluid pressure (at cell center)
  il::Array<double> hydrostatic_pressure_;
  //   fluid velocity (at cell edges)
  il::Array<double> velocity_;

  // Reynolds number (at cell edges)
  il::Array<double> Re_;

  // max relative difference on unknowns (between successive iteration),
  // ie. (u_(k+1)-u_k)/u_(k+1)
  double err_p_;  // on fluid pressure
  double err_v_;  // on av. velocity

  // number of iterations for well flow step convergence
  il::int_t ehlIts_ = 0;

  // out flow sources
  hfp2d::Sources clusters_sources_;


 public:
  ////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS
  ////////////////////////////////////////////////////////////////////////

  // default constructor
  WellSolution(){};

  //----------------------------------------------------------------------
  // normal constructor
  WellSolution(double time, double dt, il::Array<double> &p_hs,
               il::Array<double> &p, il::Array<double> &v,
               il::Array<double> &Re,
               hfp2d::Sources &outflow,
               il::int_t its, double err_p, double err_v) {
    IL_EXPECT_FAST(p.size() == p_hs.size());
    IL_EXPECT_FAST(v.size() == Re.size());
    IL_EXPECT_FAST(v.size() == p.size() - 1);

    time_ = time;
    timestep_ = dt;
    pressure_ = p;
    velocity_ = v;
    err_p_ = err_p;
    err_v_ = err_v;
    ehlIts_ = its;

    hydrostatic_pressure_ = p_hs;  // we always store the hydrostatic... we do
                                   // not recompute here... unsafe

    Re_ = Re;

   // clusters_sources_ = outflow;
  };
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  // hydrostatic constructor for time_==0
  WellSolution(WellMesh &well_mesh,
               WellInjection &well_inj,
               Fluid &fluid) {
    double g = 9.811;  // m/s^2

    // well_mesh_ = well_mesh; // Not needed to be stored here as it does not
    // change.
    // well_inj_ = well_inj;   // may be not needed to store it here.

    time_ = 0.;
    timestep_ = 0.;

    il::int_t n_el = well_mesh.numberOfElts();
    // need here to cut at plug location.....

    hydrostatic_pressure_ = il::Array<double>{n_el};

    for (il::int_t i = 0; i < n_el; ++i) {
      hydrostatic_pressure_[i] = g * fluid.fluidDensity() / 2. *
                                 (well_mesh.tvd(well_mesh.connectivity(i, 0)) +
                                  well_mesh.tvd(well_mesh.connectivity(i, 1)));
    }

    pressure_ = hydrostatic_pressure_;

    velocity_ = il::Array<double>{n_el - 1, 0.};
    Re_ = il::Array<double>{n_el - 1, 0.};

    err_p_ = 0.;
    err_v_ = 0.;

  //  clusters_sources_ = outflow;
  };
  //----------------------------------------------------------------------

  ////////////////////////////////////////////////////////////////////////
  //        public interfaces
  ////////////////////////////////////////////////////////////////////////

  double time() const { return time_; }
  double timestep() const { return timestep_; }

  il::Array<double> pressure() const { return pressure_; }
  double pressure(il::int_t e) const { return pressure_[e]; }

  il::Array<double> hydrostaticPressure() const {
    return hydrostatic_pressure_;
  }
  double hydrostaticPressure(il::int_t e) const {
    return hydrostatic_pressure_[e];
  }

  il::Array<double> velocity() const { return velocity_; }
  double velocity(il::int_t e) const { return velocity_[e]; }

  double errP() const { return err_p_; }
  double errV() const { return err_v_; }

//  il::Array<il::int_t> clusterElts() const {
//    return clusters_sources_.SourceElt();
//  }
//  il::Array<double> clusterFluxes() const {
//    return clusters_sources_.InjectionRate();
//  };

 // hfp2d::Sources clusterSources() const {return clusters_sources_;};


  ////////////////////////////////////////////////////////////////////////
  //        set functions
  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  //        methods
  ////////////////////////////////////////////////////////////////////////

  // get fluid pressure at a number of element
  il::Array<double> pressureAtElts(il::Array<il::int_t> &elt_list) {
    il::int_t nc = elt_list.size();
    IL_EXPECT_FAST(nc<=pressure_.size());

    il::Array<double> pc{nc, 0.};
    for (il::int_t i = 0; i < nc; i++) {
      pc[i] = pressure_[elt_list[i]];
    }
    return pc;
  }

  //  create json output object

  using json = nlohmann::json;

  //----------------------------------------------------------------------
  json createJsonObject() {
    json json_pressure = json::array();
    json json_velocity = json::array();
    json json_hydrostatic = json::array();
    json json_Reynolds = json::array();

    for (il::int_t m = 0; m < pressure_.size(); ++m) {
      json_pressure[m] = pressure_[m];
      json_hydrostatic[m] = hydrostatic_pressure_[m];
    }

    for (il::int_t m = 0; m < velocity_.size(); ++m) {
      json_velocity[m] = velocity_[m];
      json_Reynolds[m] = Re_[m];
    }

//    json json_hf_loc = json::array();
//    json json_hf_rate = json::array();


    json j_obj = {{"Time", time_},
                  {"Time step", timestep_},
                  {"Pressure", json_pressure},
                  {"Hydrostatic", json_hydrostatic},
                  {"Velocity", json_velocity},
                  {"Reynolds", json_Reynolds},
//                  {"Outflow element", json_hf_loc},
//                  {"Outflow", json_hf_rate},
                  {"Error pressure", err_p_},
                  {"Error velocity", err_v_},
                  {"Its well-flow", ehlIts_}};

    return j_obj;
  };

  //----------------------------------------------------------------------
  void writeToFile(std::string &filename) {
    json j_obj = createJsonObject();
    // write prettified JSON to file
    std::ofstream output(filename);
    output << std::setw(4) << j_obj << std::endl;
  };
  //----------------------------------------------------------------------
};
}

#endif  // HFPX2D_WELLSOLUTION_H
