//
// Created by lorenzo on 10/2/17.
//

#ifndef HFPX2D_INSITUSTRESS_H
#define HFPX2D_INSITUSTRESS_H

#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp2d {

// distribution of stress and pore pressure.

class InSituStress {

private:

  // stressField should be a nrow by 3 arrays with column wise (sxx, syy, sxy )
  il::Array2D<double> stress_field_;
  // pore pressure is a nrow vector
  il::Array<double> pore_pressure_;

public:

  //  empty constructor
  InSituStress(){};

  // constructor from the 2 entries
  InSituStress(il::Array2D<double> &stressField,il::Array<double> &porePress){

    this->stress_field_ = stressField;
    this->pore_pressure_ = porePress;

  }

  // copy
  InSituStress(InSituStress &theConditions){

    this->stress_field_ = theConditions.stress_field_;
    this->pore_pressure_ = theConditions.pore_pressure_;
  }


};

}
#endif //HFPX2D_INSITUSTRESS_H
