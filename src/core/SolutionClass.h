//
// Created by lorenzo on 8/30/17.
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

  // Constructor with size: it reserve the space for the problem solution (with number of elements, interpolation order and value)
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

  // Constructor with number of elements, interpolation order and value
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

  // Constructor with other Solution class
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

  //// getter - READ ONLY
  // This gets the global vector values
  il::Array<double> globalVector() const {return globalVector_;};
  // This gets a particular item of the solution vector
  double globalVector(il::int_t i) const {return globalVector_[i];};

  /// Displacement
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
    return (i<sizeDisplacement_ ? globalVector(i) : throw std::out_of_range ("Index out of range"));
  };


  /// Pressure
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
    return (i>sizeDisplacement_ ? globalVector(i) : throw std::out_of_range ("Index out of range"));
  };

  /// Reading of errors and iteration numbers
  double errConvergence() { return errConvergence_; };
  double errConvergenceDispl() { return errConvergenceDispl_; };
  double errConvergencePress() { return errConvergencePress_; };

  il::int_t numIterations() {return numIterations_; };


  /// setter - WRITE but only if the velocity and pressure are written together
  // write the new solution all together
  void updateSolution(il::Array<double> &newVector){

    // TODO: Add a check for sizeDisplacement_+sizePressure_ == newVector.size() ?? Also, other checks in the methods

    for(il::int_t i=0; i<sizeDisplacement_+sizePressure_; i++)
    {
      globalVector_[i]=newVector[i];
    }

  };


  void updateSolution(il::Array<double> &newDisplacement, il::Array<double> &newPressure){

    for(il::int_t i=0; i<sizeDisplacement_; i++)
    {
      globalVector_[i]=newDisplacement[i];
    }

    for(il::int_t i=sizeDisplacement_; i<sizeDisplacement_+sizePressure_; i++)
    {
      globalVector_[i]=newPressure[i];
    }

  };

  // update the error values
  void updateError(double error){
    errConvergence_ = error;
  };

  void updateErrorDispl(double error){
    errConvergenceDispl_ = error;
  };
  void updateErrorPress(double error){
    errConvergencePress_ = error;
  };

  // increase iteration number by one
  void increaseIteration(){
    numIterations_++;
  };

////////////////////////////////////////////////////////////////////////////////////////
/// RESET--REFACTOR
// First of all, we need only a global vector
// Initialization is done by creating a vector with the sum of all required variables
// we save only how many displacement variables and how many pressure variables are there
// then copy constructor is as easy as the other
// setter is super easy
// getter is more difficult because has to divide the vector in 2, but we do not need it many times
// Moreover, the address of each single vector is the address of the first value, so we can
// provide both copy of displ/pressure or their address.


  // Useful functions for the class

//  il::Array<double>* displacement() { return &displacement_;}

//  void displacement(int i, double val) {
//    displacement_[i]=val;
//    return;
//  }
//
//  il::Array<double> pressure() const { return pressure_;}
//  double errConvergence() const { return errConvergence_;}
//  double errConvergenceDispl() const { return errConvergenceDispl_;}
//  double errConvergencePress() const { return errConvergencePress_;}




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


// TODO: create a suite of tests to run on class Solution
// TODO: construct the derived class
// TODO: add checks for sizes of members
// TODO: add sapxy to the class (to sum and multiply multiple solutions)

};

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

// swap of solution can be done only with pointers to a solution class


#endif //HFPX2D_SOLUTIONCLASS_H
