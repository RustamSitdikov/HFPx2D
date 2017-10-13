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
  // local stress for force vector
  il::Array<double> local_stress_;

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

  double getSigmaXX(il::int_t i){
    return stress_field_(i,0);
  }

  double getSigmaYY(il::int_t i){
    return stress_field_(i,1);
  }

  double getSigmaXY(il::int_t i){
    return stress_field_(i,2);
  }

  double getPorePress(il::int_t i){
    return pore_pressure_[i];
  }

  il::Array<double> localStress(Mesh &theMesh){

    SegmentData elemSegData;

    il::StaticArray<double, 2> tempLocStr(0.);
    il::StaticArray2D<double, 2, 2> tempGlbStr(0.);

    il::Array<double> locStress(theMesh.numDisplDofs(),0.);

    for(il::int_t i; i < theMesh.numElems(); i++){

      elemSegData=get_segment_DD_data(theMesh, i, theMesh.interpOrd());

      tempGlbStr(0,0)=stress_field_(i,0);
      tempGlbStr(0,1)=stress_field_(i,2);
      tempGlbStr(1,0)=stress_field_(i,2);
      tempGlbStr(1,1)=stress_field_(i,1);

      // tempLocStr contains
      // first component: the stress normal to the segment
      // second component: the stress parallel to the segment
      tempLocStr=il::dot(tempGlbStr, elemSegData.n);

      for(il::int_t k=0; k<(theMesh.interpOrd()+1); k++){

        // collocation point k in element, component x
        locStress[theMesh.dofDispl(i,2*k)] = tempLocStr[1];

        // collocation point k in element, component y
        locStress[theMesh.dofDispl(i,2*k+1)] = tempLocStr[0];

      }

    }

    return locStress;
  }

};

}
#endif //HFPX2D_INSITUSTRESS_H
