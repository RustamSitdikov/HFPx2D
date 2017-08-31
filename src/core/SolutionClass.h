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
  // Private vector of displacements
  il::Array<double> displacement_;
  // Private vector of pressure
  il::Array<double> pressure_;

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
      displacement_{2*(1+interpolationOrder)*numberElements}, // 2 is the number of degrees of freedom at the node
      pressure_{interpolationOrder*numberElements+1},

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
      displacement_{2*(1+interpolationOrder)*numberElements, value}, // 2 is the number of degrees of freedom at the node
      pressure_{interpolationOrder*numberElements+1, value},

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
      displacement_{s.displacement_},
      pressure_{s.pressure_},

      // Error on convergence
      errConvergence_{s.errConvergence_},
      errConvergenceDispl_{s.errConvergenceDispl_},
      errConvergencePress_{s.errConvergencePress_},

      // Number of iterations
      numIterations_{s.numIterations_}
  {};

  //// getter - READ ONLY
  /// Displacement
  // This recall all displacement solution vector
  il::Array<double> displacement() const {return displacement_;};
  // This recall a particular item in the displacement vector
  double displacement(il::int_t i) const {return displacement_[i];};
  /// Pressure
  // This recall all displacement solution vector
  il::Array<double> pressure() const {return pressure_;};
  // This recall a particular item in the displacement vector
  double pressure(il::int_t i) const {return pressure_[i];};
  /// Global solution vector
  il::Array<double> global() const{
    il::int_t dsize = displacement_.size();
    il::int_t psize = pressure_.size();

    il::Array<double> globalVector {dsize+psize,0.0};
    for(il::int_t i=0; i<dsize; i++){
      globalVector[i]=displacement_[i];
    }
    for(il::int_t i=dsize; i<dsize+psize; i++){
      globalVector[i]=pressure_[i];
    }
  }

  /// setter - WRITE but only if the velocity and pressure are written together



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



};

//// DERIVED solution class with cohesive zone model
class SolutionCohesive: Solution{

public:
  il::Array<double> cohesion_values;
  il::Array<double> friction_values;

  void saveSolution(il::Array<double> displacement_)
  {
    for(il::int_t i=0; i<10; i++) {
      displacement_[i]=0.0;
    }
  }

private:
  il::Array<double> cohesion_values_;
  il::Array<double> friction_values_;

};

}

// swap of solution can be done only with pointers to a solution class


#endif //HFPX2D_SOLUTIONCLASS_H
