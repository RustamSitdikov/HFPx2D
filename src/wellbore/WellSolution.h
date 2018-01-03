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

namespace hfp2d {

//////////////////////////////////////////////////////////////////////////
  // solution for the wellbore flow
  // for the time being, a staggered scheme (w. HF EHD solver) is planned
  //////////////////////////////////////////////////////////////////////////
  class WellSolution : Solution {
    friend class Solution;

   protected:

    const double g_ = 9.811;  // m/s^2

    // current time (tn+1 = tn + timestep)
    double time_;
    // time step
    double timestep_;
    // wellbore mesh
    WellMesh well_mesh_;
    // injection parameters
    WellInjection well_inj_;
    // fluid pressure (at nodes)
    il::Array<double> pressure_;
    // av. fluid velocity (at edges)
    il::Array<double> velocity_;

    // "well_mesh_.source_location_" and "well_mesh_.injection_rate_"
    // are inherited and refer to the HF mesh.

    // pressures at the HF injection locations ("source_location_")
    // on the fracture side, so
    il::Array<double> hf_pressure_;
    // max relative difference on unknowns (between successive iteration),
    // ie. (u_(k+1)-u_k)/u_(k+1)

    double err_p_;  // on fluid pressure
    double err_v_;  // on av. velocity

    // "pressure_" is inherited from Solution class
    // and means the pressure on the wellbore side

   public:
    ////////////////////////////////////////////////////////////////////////
    //        CONSTRUCTORS
    ////////////////////////////////////////////////////////////////////////

    // default constructor
    WellSolution(){};

    // normal constructor
    WellSolution(WellMesh &well_mesh, WellInjection &well_inj, double time,
                 double dt, il::Array<double> &p, il::Array<double> &v,
                 double err_p, double err_v) {
      well_mesh_ = well_mesh;
      well_inj_ = well_inj;
      time_ = time;
      timestep_ = dt;
      pressure_ = p;
      velocity_ = v;
      err_p_ = err_p;
      err_p_ = err_v;
    };

    // hydrostatic constructor for time_==0
    WellSolution(WellMesh &well_mesh, il::Array<il::int_t> &hf_location,
                 il::int_t plug_location, il::Array<double> &coef_perf,
                 il::Array<double> &coef_tort, il::Array<double> &beta_tort,
                 Fluid &fluid, double well_inj_rate, double dt) {
      well_mesh_ = well_mesh;
      time_ = 0.;
      timestep_ = dt;

      il::int_t n_el = well_mesh_.numberOfElts();

      pressure_ = il::Array<double>{n_el};
      for (il::int_t i = 0; i < n_el; ++i) {
        pressure_[i] = g_ * fluid.fluidDensity() / 2. *
                       (well_mesh_.tvd(well_mesh_.connectivity(i, 0)) +
                        well_mesh_.tvd(well_mesh_.connectivity(i, 1)));
      }

      velocity_ = il::Array<double>{n_el - 1, 0.};

      il::int_t n_hf = hf_location.size();
      il::Array<double> hf_vol_rate{n_hf, 0.};
      il::Array<double> hf_p_drop{n_hf, 0.};

      well_inj_ =
          WellInjection(well_inj_rate, hf_location, plug_location, hf_vol_rate,
                  hf_p_drop, coef_perf, coef_tort, beta_tort);

      err_p_ = 1.;
      err_p_ = 1.;
    };

    ////////////////////////////////////////////////////////////////////////
    //        public interfaces
    ////////////////////////////////////////////////////////////////////////

    // todo: move g_ to the special file for constants???
    const double g() { return g_; }
    double time() { return time_; }
    double timestep() { return timestep_; }

    WellMesh mesh() const { return well_mesh_; }

    WellInjection wellInjection() const { return well_inj_; }
    il::Array<double> pressure() { return pressure_; }
    il::Array<double> hfPressure() { return hf_pressure_; }
    il::Array<double> velocity() { return velocity_; }
    double errP() { return err_p_; }
    double errV() { return err_v_; }

    ////////////////////////////////////////////////////////////////////////
    //        set functions
    ////////////////////////////////////////////////////////////////////////

    void setMesh(WellMesh &well_mesh) { well_mesh_ = well_mesh; };
    void setInj(WellInjection &well_inj) { well_inj_ = well_inj; };

    ////////////////////////////////////////////////////////////////////////
    //        methods
    ////////////////////////////////////////////////////////////////////////

    // solver: matching "velocity_" to the "pressure_"
    // solver: matching "injection_rate_" to the "well_inj.hf_p_drop_"
    // solver: matching "well_inj.hf_p_drop_" to "pressure_"-"hf_pressure_"


  };

}

#endif  // HFPX2D_WELLSOLUTION_H
