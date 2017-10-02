//
// Created by lorenzo on 9/18/17.
//


//
// This file is part of HFPx2DUnitTest.
//
// Created by lorenzo on 9/18/17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2DUNITTEST_SOLIDEVOLUTION_H
#define HFPX2DUNITTEST_SOLIDEVOLUTION_H

#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

namespace hfp2d {

class SolidEvolution {

//// Based on the linear cohesive zone model
private:

  il::String type_;

  il::Array<double> maximum_stress_;
  il::Array<double> maximum_opening_;
  il::Array<double> last_saved_opening_;
  il::Array<double> fracture_energy_;

public:

  SolidEvolution(SolidEvolution theSolidEvolution){

    type_ = theSolidEvolution.type_;
    maximum_stress_ = theSolidEvolution.maximum_stress_;
    maximum_opening_ = theSolidEvolution.maximum_opening_;
    fracture_energy_ = theSolidEvolution.fracture_energy_;
    last_saved_opening_ = theSolidEvolution.last_saved_opening_;
  }

  SolidEvolution(il::Array<double> failureStress, il::Array<double> decohesionOpening){

    type_="Linear CZM";

    IL_EXPECT_FAST(failureStress.size() == decohesionOpening.size());
    maximum_stress_=failureStress;
    maximum_opening_= decohesionOpening;

    for(il::int_t i=0; i < failureStress.size(); i++) {
      fracture_energy_[i] = 0.5 * failureStress[i] * decohesionOpening[i];
      last_saved_opening_[i] = 0.0;
    }

  }

  /////////// RECODE RECODE RECODE ///////////
  /*
  bool isActive(il::int_t i, double stress) {
    return ( stress[i] > maximum_stress_[i]);
  }

  bool isLoading(double current_opening){
    return ( current_opening > last_saved_opening_);
  }

  double stressLoading(double opening){
    return ((opening < maximum_opening_) ?
            maximum_stress_ * ( 1.0 - opening/maximum_opening_ ) : 0.0);
  }

  double stressUnloading(double opening){
    return ((opening < maximum_opening_) ?
            maximum_stress_ * ( opening/maximum_opening_ ) : 0.0);
  }

  double tractionSeparationLaw(double opening){

    double stress;

    if(isLoading(opening)){

      stress = stressLoading(opening);
      last_saved_opening_ = opening;

    } else {

      stress = stressUnloading(opening);

    }

    return stress;
  }
*/

  il::String getType() { return type_; }

};

/*
class LinearCZM { // Example of linear CZM

private:

  double maximum_stress_;
  double maximum_opening_;
  double last_saved_opening_;
  double fracture_energy_;


  il::String type_;

public:

  LinearCZM(){};

  LinearCZM(double failureStress, double decohesionOpening){

    type_="Linear CZM";
    maximum_stress_=failureStress;
    maximum_opening_=decohesionOpening;
    fracture_energy_=0.5 * failureStress * decohesionOpening;
    last_saved_opening_=0.0;

  }

  bool isActive(double stress) {
    return ( stress > maximum_stress_);
  }

  bool isLoading(double current_opening){
    return ( current_opening > last_saved_opening_);
  }

  double stressLoading(double opening){
    return ((opening < maximum_opening_) ?
            maximum_stress_ * ( 1.0 - opening/maximum_opening_ ) : 0.0);
  }

  double stressUnloading(double opening){
    return ((opening < maximum_opening_) ?
            maximum_stress_ * ( opening/maximum_opening_ ) : 0.0);
  }

  double tractionSeparationLaw(double opening){

    double stress;

    if(isLoading(opening)){

      stress = stressLoading(opening);
      last_saved_opening_ = opening;

    } else {

      stress = stressUnloading(opening);

    }

    return stress;
  }

  il::String getType() { return type_; }


};
*/


}

#endif //HFPX2DUNITTEST_SOLIDEVOLUTION_H
