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

// Inclusion from Inside Loop library
#include <il/Array.h>

// Inclusion from the project
#include <src/Core/Mesh.h>
#include <src/Core/Utilities.h>
#include <src/Core_dev/SolidEvolution.h>

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

  il::Array<double> opening_dd_;  // opg DD (at nodes if P1)

  il::Array<double> shear_dd_;  // shear DD (at nodes if P1)

  il::Array<double> pressure_;  // fluid pressure (at nodes)

  il::Array<double> sigma_n_;  // total normal stress component to element (at
                               // collocation points)

  il::Array<double> tau_;  // shear stress component (at collocation points)

  il::Array<int> active_set_elements_;  // active set of elements ->
                                        // Remember: an element is
                                        // 'active' if and only if both
                                        // two collocation points fail the
                                        // M-C criterion

  il::Array2D<double> tipsLocation_;  // 2D coordinate of the location of the
                                      // tips (i.e. may be inside one element in
                                      // the case of an ILSA scheme )

  il::Array<double> ribbon_tip_s_;

  // note in the case where the solution vectors are only on sub-parts of the
  // currentmesh, which may happen for cohesive zone model
  // (or lefm thru a pre-existing mesh). 2 options: either pads with zero
  // the different solution arrays, OR stored an array with active elements
  // (both for mechanics and flow). In that last case, the proper way is most
  // probably to create a derived class

  il::int_t front_its_;  // number of fracture front iterations
  il::int_t ehl_its_;    // number of ElastoHydrodynamics iterations

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
           const il::Array<double> &sigma0, const il::Array<double> &tau0,
           il::Array<int> &act_set_elmnts) {
    // todo should have checks here on dimensions with mesh etc.

    time_ = t;
    currentmesh_ = mesh;
    opening_dd_ = width;
    shear_dd_ = sheardd;
    sigma_n_ = sigma0;
    tau_ = tau0;
    pressure_ = pressure;
    active_set_elements_ = act_set_elmnts;
  };

  //#2.
  Solution(hfp2d::Mesh &mesh, double t, double dt,
           const il::Array<double> &width, const il::Array<double> &sheardd,
           const il::Array<double> &pressure, const il::Array<double> &sigma0,
           const il::Array<double> &tau0, il::int_t itsFront, il::int_t itsEHL,
           double err_front, double err_width, double err_shear, double err_p,
           il::Array<int> &act_set_elmnts) {
    // have checks here on dimensions with mesh....

    time_ = t;
    timestep_ = dt;

    currentmesh_ = mesh;

    opening_dd_ = width;
    shear_dd_ = sheardd;
    sigma_n_ = sigma0;
    tau_ = tau0;
    pressure_ = pressure;
    active_set_elements_ = act_set_elmnts;

    front_its_ = itsFront;
    ehl_its_ = itsEHL;
    err_front_ = err_front;
    err_openingDD_ = err_width;
    err_shearDD_ = err_shear;
    err_P_ = err_p;
  };

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  // get functions
  inline il::Array<double> openingDD() const { return opening_dd_; };
  inline il::Array<double> shearDD() const { return shear_dd_; };
  inline il::Array<double> pressure() const { return pressure_; };
  inline il::Array<double> sigmaN() const { return sigma_n_; };
  inline il::Array<double> tau() const { return tau_; };

  inline double openingDD(il::int_t i) { return opening_dd_[i]; };
  inline double shearDD(il::int_t i) { return shear_dd_[i]; };
  inline double pressure(il::int_t i) { return pressure_[i]; };
  inline double sigmaN(il::int_t i) { return sigma_n_[i]; };
  inline double tau(il::int_t i) { return tau_[i]; };

  inline il::Array<int> activeSetElements() const {
    return active_set_elements_;
  };

  inline il::Array2D<double> TipsLocation() const { return tipsLocation_; };
  inline il::Array<double> RibbonsDistance() const { return ribbon_tip_s_; };

  inline hfp2d::Mesh CurrentMesh() const { return currentmesh_; };

  inline double time() const { return time_; };
  inline double timestep() const { return timestep_; }

  inline double err_front() const { return err_front_; }
  inline double err_opening() const { return err_openingDD_; }
  inline double err_shear() const { return err_shearDD_; }
  inline double err_pressure() const { return err_P_; }
  inline il::int_t front_its() const { return front_its_; }
  inline il::int_t ehl_its() const { return ehl_its_; }

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  // some set functions
  void setRibbonDistances(const il::Array<double> &srt) {
    ribbon_tip_s_ = srt;
  };

  void setTipsLocation(const il::Array2D<double> &tips_xy) {
    tipsLocation_ = tips_xy;
  };

  il::Array<int> activeSetElements(
      Mesh &theMesh, Solution &SolutionAtTn, SolidEvolution &SolidEvolution,
      il::Array2D<double> &from_edge_to_coll_press) {
    // Move pore pressure from nodal points to coll points because elasticity
    // is evaluated at collocation points (-> MC criterion is evaluated at
    // collocation points!)
    il::Array<double> press_coll{2 * theMesh.numberOfElements(), 0};
    auto p_coll = il::dot(from_edge_to_coll_press, this->pressure_);
    for (il::int_t i = 0, k = 1; i < press_coll.size(); ++i, k = k + 2) {
      press_coll[i] = p_coll[k];
    }

    // Get the set of 'failed' collocation points by checking the MC criterion
    il::Array<int> failed_set_collpoints{0};
    failed_set_collpoints.reserve(2 * theMesh.numberOfElements());
    for (int j = 0, k = 0; j < 2 * theMesh.numberOfElements(); ++j) {
      if (this->tau(j) >=
          SolidEvolution.getFricCoeff(j) * (this->sigmaN(j) - press_coll[j])) {
        failed_set_collpoints.resize(k + 1);
        failed_set_collpoints[k] = j;
        k = k + 1;
      }
    }

    // Get active set of elements
    il::Array2D<int> dof_single_dd{theMesh.numberOfElements(),
                                   (theMesh.interpolationOrder() + 1), 0};
    for (int i = 0; i < theMesh.numberOfElements(); i++) {
      for (int j = 0; j < 1 * (theMesh.interpolationOrder() + 1); j++) {
        dof_single_dd(i, j) = i * 1 * (theMesh.interpolationOrder() + 1) + j;
      }
    }

    il::Array<int> set_elements{0};
    set_elements.reserve(2 * theMesh.numberOfElements());
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

  // TODO: write solution to file  -> json format

  // TODO: read from file for restart
};
}

#endif  // HFPX2DUNITTEST_SOLUTIONATT_H
