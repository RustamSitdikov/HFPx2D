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
  class WellSolution   {

   private:

    // current time (tn+1 = tn + timestep)
    double time_;
    // time step
    double timestep_;
    // wellbore wellMesh
  //  WellMesh well_mesh_;
    // injection parameters
  //  WellInjection well_inj_;

    // fluid pressure (at cell center)
    il::Array<double> pressure_;
    // hydrostatic fluid pressure (at cell center)
    il::Array<double> hydrostatic_pressure_;
    // av. fluid velocity (at cell edges)
    il::Array<double> velocity_;


    // pressures at the HF injection locations ("source_location_")
    // on the fracture side, so - is this needed here ?
    il::Array<double> hf_pressure_;

    // max relative difference on unknowns (between successive iteration),
    // ie. (u_(k+1)-u_k)/u_(k+1)
    double err_p_;  // on fluid pressure
    double err_v_;  // on av. velocity


   public:

    ////////////////////////////////////////////////////////////////////////
    //        CONSTRUCTORS
    ////////////////////////////////////////////////////////////////////////

    // default constructor
    WellSolution(){}; // is this needed?

    // normal constructor
    WellSolution(//WellMesh &well_mesh,
                 //WellInjection &well_inj,
                 double time,
                 double dt, il::Array<double> &p_hs, il::Array<double> &p,
                                      il::Array<double> &v,
                 double err_p, double err_v) {

//      well_mesh_ = well_mesh; // Not needed to be stored here as it does not change.

    //  well_inj_ = well_inj;

      time_ = time;
      timestep_ = dt;
      pressure_ = p;
      velocity_ = v;
      err_p_ = err_p;
      err_v_ = err_v;

      hydrostatic_pressure_ = p_hs; // we always store the hydrostatic... we do not recompute here... unsafe

    };

    // hydrostatic constructor for time_==0
    WellSolution(WellMesh &well_mesh,
                 WellInjection &well_inj,
                 Fluid &fluid,  double dt) {

      double g = 9.811;  // m/s^2

     // well_mesh_ = well_mesh; // Not needed to be stored here as it does not change.
     // well_inj_ = well_inj;   // may be not needed to store it here.

      time_ = 0.;
      timestep_ = dt;

      il::int_t n_el = well_mesh.numberOfElts();

      hydrostatic_pressure_ = il::Array<double>{n_el};

      for (il::int_t i = 0; i < n_el; ++i) {
        hydrostatic_pressure_[i] = g * fluid.fluidDensity() / 2. *
                       (well_mesh.tvd(well_mesh.connectivity(i, 0)) +
                           well_mesh.tvd(well_mesh.connectivity(i, 1)));
      }

      pressure_ = hydrostatic_pressure_;

      velocity_ = il::Array<double>{n_el - 1, 0.};

      err_p_ = 0.;
      err_v_ = 0.;

    };

    ////////////////////////////////////////////////////////////////////////
    //        public interfaces
    ////////////////////////////////////////////////////////////////////////

    inline double time() { return time_; }
    inline double timestep() { return timestep_; }

   // WellMesh wellMesh() const { return well_mesh_; }
 //   WellInjection wellInjection() const { return well_inj_; }

    inline il::Array<double> pressure() { return pressure_; }
    inline il::Array<double> hydrostaticPressure() { return hydrostatic_pressure_; }

    inline il::Array<double> hfPressure() { return hf_pressure_; }
    inline il::Array<double> velocity() { return velocity_; }
    inline double errP() { return err_p_; }
    inline double errV() { return err_v_; }

    ////////////////////////////////////////////////////////////////////////
    //        set functions
    ////////////////////////////////////////////////////////////////////////

    //  void setMesh(WellMesh &well_mesh) { well_mesh_ = well_mesh; };  // ??? why
  //    void setInj(WellInjection &well_inj) { well_inj_ = well_inj; };

    ////////////////////////////////////////////////////////////////////////
    //        methods
    ////////////////////////////////////////////////////////////////////////

    // solver: matching "velocity_" to the "pressure_"
    // solver: matching "injection_rate_" to the "well_inj.hf_p_drop_"
    // solver: matching "well_inj.hf_p_drop_" to "pressure_"-"hf_pressure_"

    // write solution to Json file


  };

}

#endif  // HFPX2D_WELLSOLUTION_H
