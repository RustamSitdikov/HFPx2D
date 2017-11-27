//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 10.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2DUNITTEST_SOLUTIONATT_H
#define HFPX2DUNITTEST_SOLUTIONATT_H

// Inclusion from standard library
#include <fstream>

// Inclusion from Inside Loop library
#include <il/Array.h>

// Inclusion from the project
#include <src/Core/Mesh.h>
#include <src/Core/Utilities.h>
#include <src/Core_dev/SolidEvolution.h>
#include <json.hpp>


namespace hfp2d {

class Solution {
  // Base  class for solution of coupled fluid driven fracture problem

  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
  double time_;  // current time (tn+1 = tn + timestep) at which the solution
                 // refers to

  double
      timestep_;  // time step  (of the last time step taken to arrive at time_)

  hfp2d::Mesh currentmesh_;  // the associated mesh  // should be a reference

  il::Array<double> openingDD_;  // opg DD (at nodes if P1)

  il::Array<double> shearDD_;  // shear DD (at nodes if P1)

  il::Array<double> pressure_;  // fluid pressure (at nodes)

  il::Array<double> sigma_n_;  // total normal stress component to element (at
                               // collocation points)

  il::Array<double> tau_;  // shear stress component (at collocation points)

  // DD  vector
  il::Array<double> dd_values_;

  // traction  vector
  il::Array<double> traction_values_;

  // state variable
  il::Array<double> state_1_;

  // TODO: change type to il::int_t
  il::Array<int> active_set_elements_;  // active set of elements ->
                                        // Remember: an element is
                                        // 'active' if and only if both
                                        // two collocation points fail the
                                        // M-C criterion

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

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////

 public:
  // constructor. 3 cases are only needed really
  // 1. quick construction with mesh obj., width, shear and fluid pressure,
  // sigma0 and tau0
  // 2. complete construction with also errors etc.
  // 3. complete construction from an Input file (i.e. needed for a restart
  // simulation)

  // for now, we implement #1 ans #2
  Solution() = default;

  //#1.
  Solution(hfp2d::Mesh &mesh, double t, const il::Array<double> &width,
           const il::Array<double> &sheardd, const il::Array<double> &pressure,
           const il::Array<double> &sigma0, const il::Array<double> &tau0) {
    // todo should have checks here on dimensions with mesh etc.

    time_ = t;
    currentmesh_ = mesh;
    openingDD_ = width;
    shearDD_ = sheardd;
    sigma_n_ = sigma0;
    tau_ = tau0;
    pressure_ = pressure;
  };

  //#2.    constructor w.o active set
  Solution(hfp2d::Mesh &mesh, double t, double dt,
           const il::Array<double> &width, const il::Array<double> &sheardd,
           const il::Array<double> &pressure, const il::Array<double> &sigma0,
           const il::Array<double> &tau0, il::int_t itsFront, il::int_t itsEHL,
           double err_front, double err_width, double err_shear, double err_p) {
    // have checks here on dimensions with mesh....

    time_ = t;
    timestep_ = dt;
    currentmesh_ = mesh;

    openingDD_ = width;
    shearDD_ = sheardd;
    sigma_n_ = sigma0;
    tau_ = tau0;
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
           const il::Array<double> &pressure, const il::Array<double> &sigma0,
           const il::Array<double> &tau0, const il::Array<int> &act_set_elmnts,
           il::int_t itsFront, il::int_t itsEHL, double err_front,
           double err_width, double err_shear, double err_p) {
    // have checks here on dimensions with mesh....

    time_ = t;
    timestep_ = dt;
    currentmesh_ = mesh;

    openingDD_ = width;
    shearDD_ = sheardd;
    sigma_n_ = sigma0;
    tau_ = tau0;
    pressure_ = pressure;

    active_set_elements_ = act_set_elmnts;
    frontIts_ = itsFront;
    ehlIts_ = itsEHL;
    err_front_ = err_front;
    err_openingDD_ = err_width;
    err_shearDD_ = err_shear;
    err_P_ = err_p;
  };

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  // get functions
  inline il::Array<double> openingDD() const { return openingDD_; };
  inline il::Array<double> shearDD() const { return shearDD_; };
  inline il::Array<double> pressure() const { return pressure_; };
  inline il::Array<double> sigmaN() const { return sigma_n_; };
  inline il::Array<double> tau() const { return tau_; };

  inline double openingDD(il::int_t i) { return openingDD_[i]; };
  inline double shearDD(il::int_t i) { return shearDD_[i]; };
  inline double pressure(il::int_t i) { return pressure_[i]; };
  inline double sigmaN(il::int_t i) { return sigma_n_[i]; };
  inline double tau(il::int_t i) { return tau_[i]; };

  inline il::Array<int> activeElts() const { return active_set_elements_; };
  inline il::int_t activeElts(il::int_t i) const {
    return active_set_elements_[i];
  };

  inline il::Array2D<double> tipsLocation() const { return tipsLocation_; };
  inline il::Array<double> ribbonsDistance() const { return ribbon_tips_s_; };

  inline hfp2d::Mesh currentMesh() const { return currentmesh_; };

  inline double time() const { return time_; };
  inline double timestep() const { return timestep_; }

  inline double errFront() const { return err_front_; }
  inline double errOpening() const { return err_openingDD_; }
  inline double errShear() const { return err_shearDD_; }
  inline double errPressure() const { return err_P_; }
  inline il::int_t frontIts() const { return frontIts_; }
  inline il::int_t ehlIts() const { return ehlIts_; }

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  // some set functions
  inline void setRibbonDistances(const il::Array<double> &srt) {
    ribbon_tips_s_ = srt;
  };

  inline void setTipsLocation(const il::Array2D<double> &tips_xy) {
    tipsLocation_ = tips_xy;
  };

  inline void setErrorFront(const double errF) { err_front_ = errF; };

  inline void setItsFront(const il::int_t its) { frontIts_ = its; };

  inline void setTipsVelocity(const il::Array<double> &tips_vel) {
    tips_velocity_ = tips_vel;
  };

  inline void setTimeStep(const double dt) { timestep_ = dt; };

  inline void setActiveElts(const il::Array<int> &act_set_elmnts) {
    active_set_elements_ = act_set_elmnts;
  };

  inline void setFrontPositIters(const il::int_t iters_fract_posit) {
    frontIts_ = iters_fract_posit;
  };

  il::Array<int> activeSetElements(Mesh &theMesh, Solution &SolutionAtTn,
                                   SolidEvolution &SolidEvolution,
                                   il::Array2D<double> &from_edge_to_coll_press,
                                   il::Array<double> &press_old) {
    // Move pore pressure from nodal points to coll points because elasticity
    // is evaluated at collocation points (-> MC criterion is evaluated at
    // collocation points!)
    il::Array<double> press_coll{2 * theMesh.numberOfElts(), 0};
    auto p_coll = il::dot(from_edge_to_coll_press, press_old);
    for (il::int_t i = 0, k = 1; i < press_coll.size(); ++i, k = k + 2) {
      press_coll[i] = p_coll[k];
    }

    // Get the set of 'failed' collocation points by checking the MC criterion
    il::Array<int> failed_set_collpoints{0};
    failed_set_collpoints.reserve(2 * theMesh.numberOfElts());
    for (int j = 0, k = 0; j < 2 * theMesh.numberOfElts(); ++j) {
      if (SolutionAtTn.tau(j) > SolidEvolution.getFricCoeff(j) *
                                    (SolutionAtTn.sigmaN(j) - press_coll[j])) {
        failed_set_collpoints.resize(k + 1);
        failed_set_collpoints[k] = j;
        k = k + 1;
      }
    }

    // Get active set of elements
    il::Array2D<int> dof_single_dd{theMesh.numberOfElts(),
                                   (theMesh.interpolationOrder() + 1), 0};
    for (int i = 0; i < theMesh.numberOfElts(); i++) {
      for (int j = 0; j < 1 * (theMesh.interpolationOrder() + 1); j++) {
        dof_single_dd(i, j) = i * 1 * (theMesh.interpolationOrder() + 1) + j;
      }
    }

    il::Array<int> set_elements{0};
    set_elements.reserve(2 * theMesh.numberOfElts());
    for (int l = 0, k = 0; l < failed_set_collpoints.size(); ++l) {
      set_elements.resize(k + 1);
      set_elements[l] =
          hfp2d::find_2d_integer(dof_single_dd, failed_set_collpoints[l])[0];
      k = k + 1;
    }

    auto active_set_elmnts = hfp2d::delete_duplicates_integer(set_elements);

    il::Array<int> active_set_elements{active_set_elmnts.size()};

    if (failed_set_collpoints.size() == 2 * active_set_elements.size()) {
      active_set_elements = active_set_elmnts;
    } else {
      active_set_elements.resize(active_set_elmnts.size() - 2);
      for (int i = 0, k = 1; i < active_set_elements.size(); ++i, ++k) {
        active_set_elements[i] = active_set_elmnts[k];
      }
    }

    return active_set_elements;
  };

  /// write solution to file  -> json format
  // for convenience
  using json = nlohmann::json;

  int writeToFile(std::string &filename) {

    // we output the mesh
    json json_coord = json::array();
    for (il::int_t m = 0; m < currentmesh_.coordinates().size(0); ++m) {
      json_coord[m] = {currentmesh_.coordinates(m, 0),
                       currentmesh_.coordinates(m, 1)};
    }

    //  connectivity, and dofs array
    json json_connectivity = json::array();
    json json_dof_handle_dd = json::array();
    json json_dof_handle_pres = json::array();

    for (il::int_t m = 0; m < currentmesh_.numberOfElts(); ++m) {
      json_connectivity[m] = {currentmesh_.connectivity(m, 0),
                              currentmesh_.connectivity(m, 1)};
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

    // note all the loop below should be collapsed into 1
    json json_shearDD = json::array();
    for (il::int_t m = 0; m < shearDD_.size(); ++m) {
      json_shearDD[m] = shearDD_[m];
    }

    json json_openingDD = json::array();
    for (il::int_t m = 0; m < openingDD_.size(); ++m) {
      json_openingDD[m] = openingDD_[m];
    }

    json json_pressure = json::array();
    for (il::int_t m = 0; m < pressure_.size(); ++m) {
      json_pressure[m] = pressure_[m];
    }

    json json_shear_stress = json::array();
    for (il::int_t m = 0; m < tau_.size(); ++m) {
      json_shear_stress[m] = tau_[m];
    }

    json json_normal_stress = json::array();
    for (il::int_t m = 0; m < sigma_n_.size(); ++m) {
      json_normal_stress[m] = sigma_n_[m];
    }

    // tips
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
                   {"Node Coordinates", json_coord},
                   {"Connectivity", json_connectivity},
                   {"Dof handle DD", json_dof_handle_dd},
                   {"Dof handle P", json_dof_handle_pres}};

    json j_tips = {{"Tip coordinates", json_tip_pos},
                   {"Tip velocity", json_tip_vel},
                   {"Ribbon-tip distance", json_ribbon_s},
                   {"Tip elts", json_tip_elt}};

    json j_obj = {{"Time", time_},
                  {"Time step", timestep_},
                  {"Its frac. front ", frontIts_},
                  {"Error Fracture front", err_front_},
                  {"Its EHL", ehlIts_},
                  {"Error EHL pressure", err_P_},
                  {"Error EHL opening", err_openingDD_},
                  {"Error EHL opening", err_shearDD_},
                  {"Mesh", j_mesh},
                  {"Tips", j_tips},
                  {"Shear DD", json_shearDD},
                  {"Opening DD", json_openingDD},
                  {"Fluid Pressure", json_pressure},
                  {"Shear traction", json_shear_stress},
                  {"Normal traction", json_normal_stress}};

    // write prettified JSON to file
    std::ofstream output(filename);
    output << std::setw(4) << j_obj << std::endl;

    return 0;
  };

  // TODO: read from file for restart
};
}

#endif  // HFPX2DUNITTEST_SOLUTIONATT_H
