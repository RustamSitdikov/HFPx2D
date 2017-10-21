//
// Created by lorenzo on 9/18/17.
//

//
// This file is part of HFPx2DUnitTest.
//
// Created by lorenzo on 9/18/17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_SOLIDEVOLUTION_H
#define HFPX2D_SOLIDEVOLUTION_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>
#include <il/base.h>

namespace hfp2d
{

class SolidEvolution
{

    //// Based on the linear cohesive zone model
private:
    il::String type_;

    il::Array<double> failure_stress_;
    il::Array<double> maximum_opening_;
    il::Array<double> last_saved_opening_;
    il::Array<double> fracture_energy_;

public:
    explicit SolidEvolution(const SolidEvolution &theSolidEvolution)
    {

        type_ = theSolidEvolution.type_;
        failure_stress_ = theSolidEvolution.failure_stress_;
        maximum_opening_ = theSolidEvolution.maximum_opening_;
        fracture_energy_ = theSolidEvolution.fracture_energy_;
        last_saved_opening_ = theSolidEvolution.last_saved_opening_;
    }

    explicit SolidEvolution(const double failureStress,
                            const double decohesionOpening,
                            const il::int_t numCollPoints)
    {

        type_ = "Linear CZM";

        double fracEnergy = 0.5 * failureStress * decohesionOpening;

        maximum_opening_.resize(numCollPoints);
        failure_stress_.resize(numCollPoints);
        fracture_energy_.resize(numCollPoints);
        last_saved_opening_.resize(numCollPoints);

        for (il::int_t i = 0; i < numCollPoints; i++)
        {
            failure_stress_[i] = failureStress;
            maximum_opening_[i] = decohesionOpening;
            fracture_energy_[i] = fracEnergy;
            last_saved_opening_[i] = 0.0;
        }
    }

    explicit SolidEvolution(const il::Array<double> failureStress,
                            const il::Array<double> decohesionOpening)
    {

        type_ = "Linear CZM";

        IL_EXPECT_FAST(failureStress.size() == decohesionOpening.size());
        failure_stress_ = failureStress;
        maximum_opening_ = decohesionOpening;

        // std::cout << "Here 2" << failureStress.size()<< " " <<
        // failureStress.capacity() << std::endl;
        for (il::int_t i = 0; i < failureStress.size(); i++)
        {
            fracture_energy_[i] = 0.5 * failureStress[i] * decohesionOpening[i];
            last_saved_opening_[i] = 0.0;
        }

        // std::cout << "Here 3" << std::endl;
    }

    /////////// RECODE RECODE RECODE ///////////
    bool isActive(double stress, il::int_t loc){
        return (stress > failure_stress_[loc]);
    }

    il::Array<bool> isActive(il::Array<double> stressVector)
    {
        IL_EXPECT_FAST(stressVector.size() == failure_stress_.size());
        il::Array<bool> actVec(stressVector.size());

        for(il::int_t i=0; i < stressVector.size(); i++) {
            actVec[i]=isActive(stressVector[i],i);
        }

        return actVec;
    }


    bool isLoading(double current_opening, il::int_t collPT)
    {
        return (current_opening > last_saved_opening_[collPT]);
    }

    double stressLoading(double opening, il::int_t collPT)
    {
        double stress;

        if (opening < maximum_opening_[collPT]){
            stress =
                failure_stress_[collPT] * (1.0 - opening / maximum_opening_[collPT]);
        } else {
            stress = 0.0;
        }

        return stress;
    }

    double stressUnloading(double opening,  il::int_t collPT)
    {
        double stress;

        if(opening < maximum_opening_[collPT])
        {
            stress =
                failure_stress_[collPT] * (opening / maximum_opening_[collPT]);
        } else {

            stress = 0.0;
        }

        return stress;
    }

    double tractionSeparationLaw(double opening, il::int_t collPT)
    {

        double stress;

        if (isLoading(opening,collPT))
        {

            stress = stressLoading(opening, collPT);
            last_saved_opening_[collPT] = opening;

        }
        else
        {

            stress = stressUnloading(opening,collPT);

        }

        return stress;
    }


    il::String getType() { return type_; }
    double getMaxOpening(il::int_t i) { return maximum_opening_[i]; }
    double getMaxStress(il::int_t i) { return failure_stress_[i]; }

    // This method sets the **HISTORICAL** maximum opening in the element in
    // the list crackElements, so that they are considered a crack that is
    // already open and they have cohesive stress equal to zero.
    void setInitialCrack(il::Array<il::int_t> crackElements, il::int_t DDxElem){

        for(il::int_t i=0; i < crackElements.size(); i++){

            for(il::int_t j=0; j < DDxElem; j++){

                last_saved_opening_[i * DDxElem + j] =
                                maximum_opening_[i * DDxElem + j];

            }

        }

    }

};

/*
class LinearCZM { // Example of linear CZM

private:

  double failure_stress_;
  double maximum_opening_;
  double last_saved_opening_;
  double fracture_energy_;


  il::String type_;

public:

  LinearCZM(){};

  LinearCZM(double failureStress, double decohesionOpening){

    type_="Linear CZM";
    failure_stress_=failureStress;
    maximum_opening_=decohesionOpening;
    fracture_energy_=0.5 * failureStress * decohesionOpening;
    last_saved_opening_=0.0;

  }

  bool isActive(double stress) {
    return ( stress > failure_stress_);
  }

  bool isLoading(double current_opening){
    return ( current_opening > last_saved_opening_);
  }

  double stressLoading(double opening){
    return ((opening < maximum_opening_) ?
            failure_stress_ * ( 1.0 - opening/maximum_opening_ ) : 0.0);
  }

  double stressUnloading(double opening){
    return ((opening < maximum_opening_) ?
            failure_stress_ * ( opening/maximum_opening_ ) : 0.0);
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

#endif // HFPX2D_SOLIDEVOLUTION_H
