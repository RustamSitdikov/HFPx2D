//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_SOLUTIONCLASS_H
#define HFPX2D_SOLUTIONCLASS_H


// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

#include <src/core/Mesh.h>
#include <src/core_dev/Simulation.h>

namespace hfp2d {

///// BASE solution class
class SolutionK {

private:

    hfp2d::Mesh currentmesh_;

    il::Array<double> DDvalues_;
    double pressure_;

    il::Array<double> stress_;  // current stress field
    il::Array<double> stress0_; // initial far field stress

    il::Array<il::int_t> activeList_;

    double time_;
    double timestep_;

    double minDeltaT_;
    double maxDeltaT_;

    double err_DDvalues_;
    double err_pressure_;
    double err_resDDval_;
    double err_resPress_;

    il::int_t fracFrontIter_;
    il::int_t nonLinSysIter_;
    il::int_t fracFrontMaxIter_;
    il::int_t nonLinSysMaxIter_;


public:

    SolutionK() {};

    //#1.
    SolutionK(const hfp2d::Mesh &mesh,
              const il::Array<double> &DDsolution,
              const double pressure,
              const il::Array<double> &stressAtColl,
              const il::Array<double> &sigma0,
              const il::Array<il::int_t> &actList,
              const simulationParams &simParams
    ) {

      // todo should have checks here on dimensions with mesh etc.
      currentmesh_ = mesh;

      DDvalues_ = DDsolution;
      pressure_ = pressure;

      stress_ = stressAtColl;
      stress0_ = sigma0;

      activeList_ = actList;

      time_ = simParams.t;
      timestep_ = simParams.deltat;
      minDeltaT_ = simParams.minDeltat;
      maxDeltaT_ = simParams.maxDeltat;

      err_DDvalues_ = simParams.errorOnDDs;
      err_pressure_ = simParams.errorOnPress;
      err_resDDval_ = simParams.errorOnResDDs;
      err_resPress_ = simParams.errorOnResPress;

      fracFrontIter_ = simParams.ffIter;
      fracFrontMaxIter_ = simParams.ffMaxIter;
      nonLinSysIter_ = simParams.nlIter;
      nonLinSysMaxIter_ = simParams.nlMaxIter;

    };

    il::Array<double> DDvalues() const { return DDvalues_; }
    double DDvalues(il::int_t i) const { return DDvalues_[i]; }
    double pressure() const { return pressure_; }

    il::Array<double> sigma0() const { return stress0_; }
    double sigma0(il::int_t i) const { return stress0_[i]; }
    il::Array<double> stress() const { return stress_; }
    double stress(il::int_t i) const { return stress_[i]; }

    hfp2d::Mesh CurrentMesh() const { return currentmesh_;}

    double time() const { return time_;}
    double timestep() const { return timestep_;}

    double err_DDs() const { return err_DDvalues_;}
    double err_press() const { return err_pressure_;}

    il::Array<il::int_t> activeList() const { return activeList_;}
    il::int_t activeList(il::int_t i) const { return activeList_[i];}

    SolutionK save(const hfp2d::Mesh &mesh,
              const il::Array<double> &DDsolution,
              const double pressure,
              const il::Array<double> &stressAtColl,
              const il::Array<double> &sigma0,
              const il::Array<il::int_t> &actList,
              const simulationParams &simParams
    ) {

      // todo should have checks here on dimensions with mesh etc.
      currentmesh_ = mesh;

      DDvalues_ = DDsolution;
      pressure_ = pressure;

      stress_ = stressAtColl;
      stress0_ = sigma0;

      activeList_ = actList;

      time_ = simParams.t;
      timestep_ = simParams.deltat;
      minDeltaT_ = simParams.minDeltat;
      maxDeltaT_ = simParams.maxDeltat;

      err_DDvalues_ = simParams.errorOnDDs;
      err_pressure_ = simParams.errorOnPress;
      err_resDDval_ = simParams.errorOnResDDs;
      err_resPress_ = simParams.errorOnResPress;

      fracFrontIter_ = simParams.ffIter;
      fracFrontMaxIter_ = simParams.ffMaxIter;
      nonLinSysIter_ = simParams.nlIter;
      nonLinSysMaxIter_ = simParams.nlMaxIter;

    };


};


///// BASE solution class
class Solution {
    // base  class for solution of coupled fluid driven fracture problem
private:

    double time_;  // current time (tn+1 = tn + timestep) solution
    double timestep_;  // time step increment for this solution

    hfp2d::Mesh currentmesh_;  // the associated mesh


    il::Array<double> openingDD_;  // opg DD (at nodes if P1)
    il::Array<double> shearDD_;   // shear DD (at nodes if P1)
    il::Array<double> pressure_;  // fluid pressure (at nodes)

    il::Array<double> sigma0_;  // initial normal stress at collocation points
    il::Array<double> tau0_;    // initial shear stress at collocation points

    il::Array2D<double> tipsLocation_;  // 2D coordinate of the tips location
    // (i.e. may be inside one element in the case of an ILSA scheme )

    il::Array<double> ribbon_tip_s_;

    // note in the case where the solution vectors are only on sub-parts of the
    // currentmesh, which may happen for cohesive zone model
    // (or lefm thru a pre-existing mesh). 2 options: either pads with zero
    // the different solution arrays, OR stored an array with active elements
    // (both for mechanics and flow). In that last case, the proper way is most
    // probably to create a derived class

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
    // 1. quick construction with mesh obj., width, shear and fluid pressure,
    // sigma0 and tau0
    // 2. complete construction with also errors etc.
    // 3. complete construction from an input file (i.e. needed for a restart
    // simulation)

    // for now, we implement #1 ans #2
    Solution(){};

    //#1.
    Solution(hfp2d::Mesh &mesh, double t, const il::Array<double> &width,
             const il::Array<double> &sheardd,
             const il::Array<double> &pressure,
             const il::Array<double> &sigma0, const il::Array<double> &tau0) {

      // todo should have checks here on dimensions with mesh etc.
      time_ = t;
      currentmesh_ = mesh;
      openingDD_ = width;
      shearDD_ = sheardd;
      sigma0_ = sigma0;
      tau0_ = tau0;
      pressure_=pressure;
    };

    //#2.
    Solution(hfp2d::Mesh &mesh, double t, double dt,
             const il::Array<double> &width, const il::Array<double> &sheardd,
             const il::Array<double> &pressure,
             const il::Array<double> &sigma0, const il::Array<double> &tau0,
             il::int_t itsFront, il::int_t itsEHL, double err_front,
             double err_width, double err_shear, double err_p) {

      // have checks here on dimensions with mesh....

      time_ = t;
      timestep_ = dt;

      currentmesh_ = mesh;

      openingDD_ = width;
      shearDD_ = sheardd;
      sigma0_ = sigma0;
      tau0_ = tau0;
      pressure_=pressure;

      frontIts_ = itsFront;
      ehlIts_ = itsEHL;
      err_front_ = err_front;
      err_openingDD_ = err_width;
      err_shearDD_ = err_shear;
      err_P_ = err_p;

    };

    // some set functions
    void setRibbonDistances(const il::Array<double> &srt){
      ribbon_tip_s_=srt;
    };

    void setTipsLocation(const il::Array2D<double> &tips_xy){
      tipsLocation_=tips_xy;
    };



    // get functions
    il::Array<double> openingDD() const { return openingDD_; };
    il::Array<double> shearDD() const { return shearDD_; };
    il::Array<double> pressure() const { return pressure_; };
    il::Array<double> sigma0() const { return sigma0_; };
    il::Array<double> tau0() const { return tau0_; };

    il::Array2D<double> TipsLocation() const {return tipsLocation_;};
    il::Array<double> RibbonsDistance() const { return ribbon_tip_s_;};

    hfp2d::Mesh CurrentMesh() const { return currentmesh_;};


    double time() const { return time_;};
    double timestep() const { return timestep_;}

    double err_front() const { return err_front_;}
    double err_opening() const { return err_openingDD_;}
    double err_shear() const { return err_shearDD_;}
    double err_pressure() const { return err_P_;}
    il::int_t  front_its() const { return frontIts_;}
    il::int_t  ehl_its() const { return ehlIts_;}

//////////////////////////////////////////////////////////////////////////////

    // methods:

    // write solution to file  -> json format

    // read from file for restart

};

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///// OLD solution class
class SolutionOLD {

  // TODO: create a suite of tests to run on class SolutionOLD
  // TODO: construct the derived classes
  // TODO: add sapxy to the class (to sum and multiply multiple solutions)
  // TODO: use (for many of the operations) the pointers rather than the copy itself
  // TODO: swap of solutions to be done only with pointers to two solution classes

private:

  //todo : we may remove this.
  // Private, number of displacement variables
  il::int_t sizeDisplacement_;
  // Private vector of pressure
  il::int_t sizePressure_;

  // Private global vector of solution
  il::Array<double> globalVector_;

  // Private errors on convergence
  double errConvergence_;
  double errConvergenceDispl_;
  double errConvergencePress_;

  // Private number of iterations
  il::int_t numIterations_;

 public:

  /// CONSTRUCTORS
  // Constructor with size: it reserve the space for the problem solution
  // (with number of elements, interpolation order BUT NO value)
  explicit SolutionOLD(const il::int_t numberElements, const il::int_t interpolationOrder) :

      // SolutionOLD vectors
  // TODO: correct this because it is not true with multiple fractures
  // if there are 3 fractures the number of pressure dofs is NumElements+3
      sizeDisplacement_ {2*(1+interpolationOrder)*numberElements}, // 2 is the number of degrees of freedom at the node
      sizePressure_ {interpolationOrder*numberElements+1},
      globalVector_ {sizeDisplacement_+sizePressure_},

      // Error on convergence
      errConvergence_{},
      errConvergenceDispl_{},
      errConvergencePress_{},

      // Number of iterations
      numIterations_{}
  {};

  // Constructor with number of elements, interpolation order AND value
  explicit SolutionOLD(const il::int_t numberElements, const il::int_t interpolationOrder, const double value) :

      // SolutionOLD vectors
  // TODO: correct this because it is not true with multiple fractures
      sizeDisplacement_ {2*(1+interpolationOrder)*numberElements}, // 2 is the number of degrees of freedom at the node
      sizePressure_ {interpolationOrder*numberElements+1},
      globalVector_ {sizeDisplacement_+sizePressure_, value},

      // Error on convergence
      errConvergence_{value},
      errConvergenceDispl_{value},
      errConvergencePress_{value},

      // Number of iterations
      numIterations_{0}
  {};

  // Constructor copying from other SolutionOLD object
  SolutionOLD(const SolutionOLD &s) :

  // SolutionOLD vectors
      sizeDisplacement_ {s.sizeDisplacement_}, // 2 is the number of degrees of freedom at the node
      sizePressure_ {s.sizePressure_},
      globalVector_ {s.globalVector_},

      // Error on convergence
      errConvergence_{s.errConvergence_},
      errConvergenceDispl_{s.errConvergenceDispl_},
      errConvergencePress_{s.errConvergencePress_},

      // Number of iterations
      numIterations_{s.numIterations_}
  {};


  //////////////////////////////// Interfaces ////////////////////////////////
  //// getter - READ ONLY
  // This gets the global vector values
  il::Array<double> globalVector() const {return globalVector_;};
  // This gets a particular item of the solution vector
  double globalVector(il::int_t i) const {return globalVector_[i];};

  // Displacement
  // This recall all displacement solution vector
  il::Array<double> displacement() const {
    il::Array<double> displacement_temp{sizeDisplacement_};

    for(il::int_t i=0; i < sizeDisplacement_; i++)
    {
      displacement_temp[i]=globalVector(i);
    }

    return displacement_temp;
  };

  // This recall an element in the displacement vector
  double displacement(il::int_t i) const {

    IL_EXPECT_FAST(i<sizeDisplacement_);

    return (globalVector(i));
  };

  // Pressure
  // This recall all pressure solution vector
  il::Array<double> pressure() const {
    il::Array<double> pressure_temp{sizeDisplacement_};

    for(il::int_t i=sizeDisplacement_; i < sizeDisplacement_+sizePressure_; i++)
    {
      pressure_temp[i]=globalVector(i);
    }

    return pressure_temp;
  };

  // This recall an element in the pressure vector
  double pressure(il::int_t i) const {

    IL_EXPECT_FAST(i<sizePressure_);

    return globalVector(i+sizeDisplacement_);
  };

  // Reading of errors and iteration numbers
  double errConvergence() const { return errConvergence_; };
  double errConvergenceDispl() const { return errConvergenceDispl_; };
  double errConvergencePress() const { return errConvergencePress_; };

  il::int_t numIterations() const {return numIterations_; };


  //// setter - WRITE but only if the velocity and pressure are written together
  // write the new solution all together
  void updateSolution(il::Array<double> const &newVector){

    // Check correct size of update
    IL_EXPECT_FAST(globalVector_.size()==newVector.size());

    for(il::int_t i=0; i<globalVector_.size(); i++)
    {
      globalVector_[i]=newVector[i];
    }

  };

  // update complete solution vector with displacement and pressures
  void updateSolution(il::Array<double> const &newDisplacement, il::Array<double> const &newPressure){

    IL_EXPECT_FAST(newDisplacement.size()==sizeDisplacement_ && newPressure.size()==sizePressure_)
    for(il::int_t i=0; i<sizeDisplacement_; i++)
    {
      globalVector_[i]=newDisplacement[i];
    }

    for(il::int_t i=sizeDisplacement_; i<sizeDisplacement_+sizePressure_; i++)
    {
      globalVector_[i]=newPressure[i];
    }

  };

  /// Placeholder for blas saxpy and other operations on the solution class


  // update the error values
  void updateError(double const error){
    errConvergence_ = error;
  };

  void updateErrorDispl(double const error){
    errConvergenceDispl_ = error;
  };
  void updateErrorPress(double const error){
    errConvergencePress_ = error;
  };

  // increase iteration number by one
  void increaseIteration(){
    numIterations_++;
  };

  // Reference (Address) of SolutionOLD
//  SolutionOLD *getAddress(){
//    return this; // we return the address of the solution object
//  }
//
//  SolutionOLD &getValue(){
//    return *this; // we return the value of the solution object
//  }
//
// This functions will be helpful later when, in the solution procedure,
// the addresses of solutions at n and n+1 will be changed, not their values

};


//////////////////////////////// Derived classes ////////////////////////////////

//// DERIVED solution class with cohesive zone model
class SolutionCohesive: SolutionOLD{

public:
  il::Array<double> cohesion_values;
  il::Array<double> friction_values;

private:
  il::Array<double> cohesion_values_;
  il::Array<double> friction_values_;

};





}


#endif //HFPX2D_SOLUTIONCLASS_H
