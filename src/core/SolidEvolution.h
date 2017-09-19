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

namespace hfp2d {



/*class SolidEvolution {
private:
  il::String type_;

  double stress_threshold_;

public:

  il::String getType() { return type_; }

  virtual bool isActive(double stressMeasure){
    return (stressMeasure > stress_threshold_);
  }

  *//*virtual SolidEvolution(const il::String &type, const double stress_threshold){
    type_=type;
    stress_threshold_=stress_threshold;
  }*//*


  SolidEvolution() {};

};*/


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



}

#endif //HFPX2DUNITTEST_SOLIDEVOLUTION_H
