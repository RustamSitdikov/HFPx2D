//
// Created by lorenzo on 10/2/17.
//

#ifndef HFPX2DUNITTEST_CONDITIONS_H
#define HFPX2DUNITTEST_CONDITIONS_H

#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

namespace hfp2D {

class Conditions {

private:
  il::Array2D<double> stress_field_;
  il::Array<double> pore_pressure_;

public:

  Conditions(il::Array2D<double> &stressField,il::Array<double> &porePress){

    stress_field_ = stressField;
    pore_pressure_ = porePress;

  }

  Conditions(Conditions &theConditions){

    this->stress_field_ = theConditions.stress_field_;
    this->pore_pressure_ = theConditions.pore_pressure_;

  }


};

}
#endif //HFPX2DUNITTEST_CONDITIONS_H
