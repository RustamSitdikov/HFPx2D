//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 07.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_PROPERTIES_H
#define HFPX2DUNITTEST_PROPERTIES_H

// Properties/Material ID
// - Rock elastic properties
// - Existing fluid properties
// - Solid non-linearity in fault
// - Flow in fault non-linearity
//   - Initial Permeability
//   - Permeability change
//   - Coupled opening/permeability effect

#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>


namespace hfp2d {

// for each fracture ID we will provide a properties object which is


class Properties {

private:
/*
  Solid fault_solid;
  Fluid

  SolidEvolution
  FluidEvolution




public:

*/

};



class Solid {
  private:
    // Global domain properties
    double young_; // Young modulus
    double poiss_; // Poisson ratio
    double bulkm_; // Bulk Modulus
    double lame1_; // Lame' first parameter
    double lame2_; // Shear modulus or Lame' second parameter

  /////////// SETTER OF PARAMETERS
  void setSolidParameters1(double YoungModulus, double PoissonRatio){

    young_ = YoungModulus;
    poiss_ = PoissonRatio;

    bulkm_ = young_ / (3.0* (1.0 - 2.0 * poiss_));
    lame1_ = young_ / ((1.0 + poiss_)*(1.0 - 2.0 * poiss_));
    lame2_ = young_ / (2.0 * (1 + poiss_));

  }

  void setSolidParameters2(double BulkModulus, double ShearModulus){

    bulkm_ = BulkModulus;
    lame2_ = ShearModulus;

    young_ = 9.0 * bulkm_ * lame2_ / (3.0 * bulkm_ + lame2_);
    poiss_ = (3.0 * bulkm_ - 2.0 * lame2_) / (2.0 * (3.0 * bulkm_ + lame2_ ));
    lame1_ = bulkm_ - 2.0 * lame2_ / 3.0;

  }

  void setSolidParameters3(double Lame1, double Lame2){

    lame1_ = Lame1;
    lame2_ = Lame2;

    bulkm_ = lame1_ + 2.0 * lame2_ / 3.0;
    young_ = lame2_ * (3.0 * lame1_ + 2.0 * lame2_) / (lame1_ + lame2_);
    poiss_ = lame1_ / (2.0 * (lame1_ + lame2_));

  }


  /////////// GETTER OF PARAMETERS
  double youngModulus(){ return young_; }
  double poissonRatio(){ return poiss_; }
  double bulkModulus(){ return bulkm_; }
  double lame1Parameter(){ return lame1_; }
  double lame2Parameter(){ return lame2_; }
  double shearModulus(){ return lame2_; }


};

class Fluid {

private:
  double fluid_density_;
  double fluid_viscosity_;
  double fluid_compressibility_;

public:

  void setFluidParameters(double density, double viscosity, double compressibility) {

    fluid_density_ = density;
    fluid_viscosity_ = viscosity;
    fluid_compressibility_ = compressibility;
  }

  double fluidDensity() { return fluid_density_; };
  double fluidViscosity() { return fluid_viscosity_; };
  double fluidCompressibility() { return fluid_compressibility_; };

};



class SolidEvolution {
private:
  il::String type_;

  double stress_threshold_;

public:

  il::String getType() { return type_; }

  virtual bool isActive(double stressMeasure){
    return (stressMeasure > stress_threshold_);
  }

  /*virtual SolidEvolution(const il::String &type, const double stress_threshold){
    type_=type;
    stress_threshold_=stress_threshold;
  }*/


  SolidEvolution() {};

};



class czm : SolidEvolution { // Example of linear CZM

private:

  double maximum_stress_;
  double maximum_opening_;
  double last_saved_opening_;
  double fracture_energy_;

public:

  czm(double failureStress, double decohesionOpening){

    maximum_stress_=failureStress;
    maximum_opening_=decohesionOpening;
    fracture_energy_=0.5 * failureStress * decohesionOpening;

  }

  bool isActive(double stress) override {
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

    }else{

      stress = stressUnloading(opening);

    }

    return stress;
  }

};


class FlowEvolution{  // Example of constant permeability
private:
  il::String type_;
  double permeability_;
  double initial_permeability_;

public:
  il::String getType() { return type_; }

  void setInitialPermeability(double k){
    permeability_=k;
    initial_permeability_=k;
  }

  double getPermeability(){
    return permeability_;
  }

};


};


#endif //HFPX2DUNITTEST_PROPERTIES_H
