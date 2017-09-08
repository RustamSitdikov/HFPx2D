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

namespace hfp2d {

///// BASE solution class
class Solution {

  // TODO: create a suite of tests to run on class Solution
  // TODO: construct the derived classes
  // TODO: add sapxy to the class (to sum and multiply multiple solutions)
  // TODO: use (for many of the operations) the pointers rather than the copy itself
  // TODO: swap of solutions to be done only with pointers to two solution classes

private:

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
  explicit Solution(const il::int_t numberElements, const il::int_t interpolationOrder) :

      // Solution vectors
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
  explicit Solution(const il::int_t numberElements, const il::int_t interpolationOrder, const double value) :

      // Solution vectors
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

  // Constructor copying from other Solution object
  Solution(const Solution &s) :

  // Solution vectors
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

  // Reference (Address) of Solution
//  Solution *getAddress(){
//    return this; // we return the address of the solution object
//  }
//
//  Solution &getValue(){
//    return *this; // we return the value of the solution object
//  }
//
// This functions will be helpful later when, in the solution procedure,
// the addresses of solutions at n and n+1 will be changed, not their values

};


//////////////////////////////// Derived classes ////////////////////////////////

//// DERIVED solution class with cohesive zone model
class SolutionCohesive: Solution{

public:
  il::Array<double> cohesion_values;
  il::Array<double> friction_values;

private:
  il::Array<double> cohesion_values_;
  il::Array<double> friction_values_;

};

}


#endif //HFPX2D_SOLUTIONCLASS_H
