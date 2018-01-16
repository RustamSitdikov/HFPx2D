//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 10.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#ifndef HFPX2DUNITTEST_SOLUTIONATT_H
#define HFPX2DUNITTEST_SOLUTIONATT_H

#include <fstream>

#include <il/Array.h>
#include <il/container/string/String.h>

#include <src/core/Mesh.h>
#include <src/util/json.hpp>

namespace hfp2d {

class Solution {
  //  class for solution of coupled fluid driven fracture problem
 private:

  double time_;  // current time (tn+1 = tn + timestep) at which the solution
                 // refers to
  double
      timestep_;  // time step  (of the last time step taken to arrive at time_)

  hfp2d::Mesh currentmesh_;  // the associated wellMesh  // should be a reference

  il::Array<double> openingDD_;  // opg DD (at nodes if P1)

  il::Array<double> shearDD_;  // shear DD (at nodes if P1)

  il::Array<double> pressure_;  // fluid pressure (at nodes)

  il::Array<double> sigma_n_o_;  // in-situ normal stress at collocation points
  il::Array<double> tau_o_;      // in-situ  shear stress at collocation points

  // total normal traction component to element (at collocation points)
  il::Array<double> sigma_n_;

  //  total shear traction component (at collocation points)
  il::Array<double> tau_;

  // DD  vector
  il::Array<double> dd_values_;

  // traction  vector
  il::Array<double> traction_values_;

  // state variable
  il::Array<double> state_1_;

  // active set of elements -> e.g. satisfying a yield criteria
  il::Array<il::int_t> active_set_elements_;

  il::Array2D<double> tipsLocation_;  // 2D coordinate of the location of the
                                      // tips (i.e. may be inside one element in
                                      // the case of an ILSA scheme )

  il::Array<double> tips_velocity_;

  il::Array<double> ribbon_tips_s_;

  il::int_t frontIts_;  // number of fracture front iterations
  il::int_t ehlIts_;    // number of ElastoHydrodynamics iterations

  // relative difference on unknowns (between successive iteration) , ie.
  // (u^k+1-u^k)/u^k+1

  double err_front_;  // on fracture front location -> humm this could be an
                      // array if there is several fracture (to detect the one
                      // converging in some cases ?)
  double err_openingDD_;  // on opening dds
  double err_shearDD_;    //  on shear dds
  double err_P_;          // on fluid pressure

  // residuals of ehl system (last front iteration)
  double residuals_ehl_;

 public:
  // constructor. 3 cases are only needed really
  // 1. quick construction with wellMesh obj., width, shear and fluid pressure,
  // sigma0 and tau0
  // 2. complete construction with also errors etc.
  // 3. complete construction from an input file (i.e. needed for a restart
  // simulation)

  // for now, we implement #1 and #2
  Solution(){};

  // todo : do also the cae where current stress are also stored....
  //#1. constructor w.o active set and without current stress at collocation
  //points
  Solution(hfp2d::Mesh &mesh, double t, const il::Array<double> &width,
           const il::Array<double> &sheardd, const il::Array<double> &pressure,
           const il::Array<double> &sigma0, const il::Array<double> &tau0) {
    // todo should have checks here on dimensions with wellMesh etc.
    time_ = t;
    currentmesh_ = mesh;
    openingDD_ = width;
    shearDD_ = sheardd;
    //    sigma_n_ = sigma0;
    sigma_n_o_ = sigma0;

    //    tau_ = tau0;
    tau_o_ = tau0;
    pressure_ = pressure;
  };

  //#2.    constructor w.o active set and without current stress at collocation
  //points
  Solution(hfp2d::Mesh &mesh, double t, double dt,
           const il::Array<double> &width, const il::Array<double> &sheardd,
           const il::Array<double> &pressure, const il::Array<double> &sigma0,
           const il::Array<double> &tau0, il::int_t itsFront, il::int_t itsEHL,
           double err_front, double err_width, double err_shear, double err_p) {
    // have checks here on dimensions with wellMesh....

    time_ = t;
    timestep_ = dt;
    currentmesh_ = mesh;

    openingDD_ = width;
    shearDD_ = sheardd;
    //    sigma_n_ = sigma0;
    //    tau_ = tau0;
    sigma_n_o_ = sigma0;
    tau_o_ = tau0;

    pressure_ = pressure;

    frontIts_ = itsFront;
    ehlIts_ = itsEHL;
    err_front_ = err_front;
    err_openingDD_ = err_width;
    err_shearDD_ = err_shear;
    err_P_ = err_p;
  };

  //#3.    constructor w.  active set
  Solution(hfp2d::Mesh &mesh, double t, double dt,
           const il::Array<double> &width, const il::Array<double> &sheardd,
           const il::Array<double> &pressure, const il::Array<double> &sigmaC,
           const il::Array<double> &tauC, const il::Array<double> &sigma0,
           const il::Array<double> &tau0, il::Array<il::int_t> &act_set_elmnts,
           il::int_t itsFront, il::int_t itsEHL, double err_front,
           double err_width, double err_shear, double err_p) {
    // have checks here on dimensions with wellMesh....

    time_ = t;
    timestep_ = dt;
    currentmesh_ = mesh;

    openingDD_ = width;
    shearDD_ = sheardd;
    sigma_n_ = sigmaC;
    tau_ = tauC;

    sigma_n_o_ = sigma0;
    tau_o_ = tau0;

    pressure_ = pressure;

    active_set_elements_ = act_set_elmnts;
    frontIts_ = itsFront;
    ehlIts_ = itsEHL;
    err_front_ = err_front;
    err_openingDD_ = err_width;
    err_shearDD_ = err_shear;
    err_P_ = err_p;
  };

  // some set functions
   void setRibbonDistances(const il::Array<double> &srt) {
    ribbon_tips_s_ = srt;
  };

   void setTipsLocation(const il::Array2D<double> &tips_xy) {
    tipsLocation_ = tips_xy;
  };

   void setErrorFront(const double errF) { err_front_ = errF; };

   void setItsFront(const il::int_t its) { frontIts_ = its; };

   void setTipsVelocity(const il::Array<double> &tips_vel) {
    tips_velocity_ = tips_vel;
  };

   void setTimeStep(const double dt) { timestep_ = dt; };

   void setActiveElts(const il::Array<il::int_t> &act_set_elmnts) {
    active_set_elements_ = act_set_elmnts;
  }

  /////////////////////////////////////////////////////////////////////////////
  // get functions

   il::Array<double> openingDD() const { return openingDD_; };
   il::Array<double> shearDD() const { return shearDD_; };
   il::Array<double> pressure() const { return pressure_; };
   il::Array<double> sigma0() const { return sigma_n_o_; };
   il::Array<double> tau0() const { return tau_o_; };

   double pressure(il::int_t & k) const { return pressure_[k];};

   il::Array<il::int_t> activeElts() const {
    return active_set_elements_;
  };

   il::Array2D<double> tipsLocation() const { return tipsLocation_; };
   il::Array<double> ribbonsDistance() const { return ribbon_tips_s_; };
   il::Array<double> tipsVelocity() const { return tips_velocity_; };

   hfp2d::Mesh currentMesh() const { return currentmesh_; };

   double time() const { return time_; };
   double timestep() const { return timestep_; }

   double errFront() const { return err_front_; }
   double errOpening() const { return err_openingDD_; }
   double errShear() const { return err_shearDD_; }
   double errPressure() const { return err_P_; }

   il::int_t frontIts() const { return frontIts_; }
   il::int_t ehlIts() const { return ehlIts_; }

  //////////////////////////////////////////////////////////////////////////////
  //  METHODS :
  //////////////////////////////////////////////////////////////////////////////

  // write solution to file  -> json format
  // for convenience
  using json = nlohmann::json;

  //----------------------------------------------------------------------
  json createJsonObject(){
    // we output the wellMesh
    json json_coord = json::array();
    for (il::int_t m = 0; m < currentmesh_.coordinates().size(0); ++m) {
      json_coord[m] = {currentmesh_.coordinates(m, 0),
                       currentmesh_.coordinates(m, 1)};
    }
    //  connectivity, and dofs array
    json json_connectivity = json::array();
    json json_dof_handle_dd = json::array();
    json json_dof_handle_pres = json::array();
    json json_mat_id = json::array();

    for (il::int_t m = 0; m < currentmesh_.numberOfElts(); ++m) {
      json_connectivity[m] = {currentmesh_.connectivity(m, 0),
                              currentmesh_.connectivity(m, 1)};
      json_mat_id[m] = currentmesh_.matid(m);
      if (currentmesh_.interpolationOrder() == 0) {
        json_dof_handle_dd[m] = {currentmesh_.dofDD(m, 0),
                                 currentmesh_.dofDD(m, 1)};
        json_dof_handle_pres[m] = currentmesh_.dofPress(m, 0);
      } else if (currentmesh_.interpolationOrder() == 1) {
        json_dof_handle_dd[m] = {
            currentmesh_.dofDD(m, 0), currentmesh_.dofDD(m, 1),
            currentmesh_.dofDD(m, 2), currentmesh_.dofDD(m, 3)};
        json_dof_handle_pres[m] = {currentmesh_.dofPress(m, 0),
                                   currentmesh_.dofPress(m, 1)};
      }
    }

    // note all the loop below should be collapsed into a single one
    json json_shearDD = json::array();
    json json_openingDD = json::array();
    for (il::int_t m = 0; m < shearDD_.size(); ++m) {
      json_shearDD[m] = shearDD_[m];
      json_openingDD[m] = openingDD_[m];
    }

    json json_pressure = json::array();
    for (il::int_t m = 0; m < pressure_.size(); ++m) {
      json_pressure[m] = pressure_[m];
    }

    json json_shear_stress = json::array();
    json json_normal_stress = json::array();
    for (il::int_t m = 0; m < tau_.size(); ++m) {
      json_shear_stress[m] = tau_[m];
      json_normal_stress[m] = sigma_n_[m];
    }

    json json_shear_stress_o = json::array();
    json json_normal_stress_o = json::array();
    for (il::int_t m = 0; m < tau_o_.size(); ++m) {
      json_shear_stress_o[m] = tau_o_[m];
      json_normal_stress_o[m] = sigma_n_o_[m];
    }

    // tips structure etc.
    json json_tip_pos = json::array();
    json json_tip_vel = json::array();
    json json_ribbon_s = json::array();
    json json_tip_elt = json::array();

    for (il::int_t m = 0; m < tips_velocity_.size(); ++m) {
      json_tip_pos[m] = {tipsLocation_(m, 0), tipsLocation_(m, 1)};
      json_tip_vel[m] = tips_velocity_[m];
      json_ribbon_s[m] = ribbon_tips_s_[m];
      json_tip_elt[m] = currentmesh_.tipElts(m);
    }

    json j_mesh = {{"Interpolation order", currentmesh_.interpolationOrder()},
                   {"Node coordinates", json_coord},
                   {"Connectivity", json_connectivity},
                   {"Material ID", json_mat_id},
                   {"Dof handle DD", json_dof_handle_dd},
                   {"Dof handle P", json_dof_handle_pres}};

    json j_tips = {{"Tip coordinates", json_tip_pos},
                   {"Tip velocity", json_tip_vel},
                   {"Ribbon-tip distance", json_ribbon_s},
                   {"Tip elts", json_tip_elt}};

    json j_obj = {{"Time", time_},
                  {"Time step", timestep_},
                  {"Its frac. front", frontIts_},
                  {"Error Fracture front", err_front_},
                  {"Its EHL", ehlIts_},
                  {"Error EHL pressure", err_P_},
                  {"Error EHL opening", err_openingDD_},
                  {"Error EHL opening", err_shearDD_},
                  {"Mesh", j_mesh},
                  {"Tips", j_tips},
                  {"Shear DD", json_shearDD},
                  {"Opening DD", json_openingDD},
                  {"Fluid pressure", json_pressure},
                  {"Initial shear traction", json_shear_stress_o},
                  {"Initial normal traction", json_normal_stress_o},
                  {"Shear traction", json_shear_stress},
                  {"Normal traction", json_normal_stress}};

    return j_obj;
  }
  //----------------------------------------------------------------------


  //----------------------------------------------------------------------
  void writeToFile(std::string &filename) {

    json j_obj=createJsonObject();
    // write prettified JSON to file
    std::ofstream output(filename);
    output << std::setw(4) << j_obj << std::endl;

  };
  //----------------------------------------------------------------------

};

// todo read json from file for restart

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
}

#endif  // HFPX2DUNITTEST_SOLUTIONATT_H
