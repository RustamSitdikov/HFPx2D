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

namespace hfp2d {

class Properties {

};

class Solid {

private:
  double young_; // Young modulus
  double poiss_; // Poisson ratio
  double bulkm_; // Bulk Modulus
  double lame1_; // Lame' first parameter
  double lame2_; // Shear modulus or Lame' second parameter

public:

  /////////// SETTER OF PARAMETERS
  void setSolidParamters1(double YoungModulus, double PoissonRatio){

    young_ = YoungModulus;
    poiss_ = PoissonRatio;

    bulkm_ = young_ / (3.0* (1.0 - 2.0 * poiss_));
    lame1_ = young_ / ((1.0 + poiss_)*(1.0 - 2.0 * poiss_));
    lame2_ = young_ / (2.0 * (1 + poiss_));

  }

  void setSolidParamters2(double BulkModulus, double ShearModulus){

    bulkm_ = BulkModulus;
    lame2_ = ShearModulus;

    young_ = 9.0 * bulkm_ * lame2_ / (3.0 * bulkm_ + lame2_);
    poiss_ = (3.0 * bulkm_ - 2.0 * lame2_) / (2.0 * (3.0 * bulkm_ + lame2_ ));
    lame1_ = bulkm_ - 2.0 * lame2_ / 3.0;

  }

  void setSolidParamters3(double Lame1, double Lame2){

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

}

#endif //HFPX2DUNITTEST_PROPERTIES_H
