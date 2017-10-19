//
// Created by lorenzo on 10/2/17.
//

#ifndef HFPX2DUNITTEST_InSituStress_H
#define HFPX2DUNITTEST_InSituStress_H

#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

namespace hfp2d {

// distribution of stress and pore pressure....
// kept as a class for now... or move it as a struct ???

class InSituStress {

private:

  // need to a have reference to the mesh !!!

  // stressField should be a nrow by 3 arrays with column wise (sxx, syy, sxy )
  il::Array2D<double> stress_field_;
  // pore pressure is a nrow vector
  il::Array<double> pore_pressure_;

public:

  //  empty constructor
  InSituStress(){};

  // constructor from the 2 entries
  InSituStress(il::Array2D<double> &stressField,il::Array<double> &porePress){

     stress_field_ = stressField;
     pore_pressure_ = porePress;

  }

  // copy ......... we can do that with std::move
  InSituStress(InSituStress &theConditions){

     stress_field_ = theConditions.stress_field_;
     pore_pressure_ = theConditions.pore_pressure_;
  }

// Method  to compute normal and shear traction on element of the mesh.....


};

}
#endif //HFPX2DUNITTEST_CONDITIONS_H
