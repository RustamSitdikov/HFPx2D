//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_ELASTICPROPERTIES_H
#define HFPX2D_ELASTICPROPERTIES_H

namespace hfp2d {

class ElasticProperties {

 private:
  double young_; // Young's modulus
  double poiss_; // Poisson's ratio
  double lame1_; // Lame's first parameter
  double lame2_; // Lame's second parameter or shear modulus
  double bulkm_; // Bulk modulus
  double youngPS_; // Plane-strain Young's modulus

 public:
  ElasticProperties(){};

  // Creation of elastic properties class from Young's modulus and Poisson's Ratio
  ElasticProperties(double YoungModulus, double PoissonRatio) {

    young_ = YoungModulus;
    poiss_ = PoissonRatio;

    bulkm_ = young_ / (3.0* (1.0 - 2.0 * poiss_));
    lame1_ = young_ / ((1.0 + poiss_)*(1.0 - 2.0 * poiss_));
    lame2_ = young_ / (2.0 * (1 + poiss_));

    youngPS_ = young_ / (1.0 - poiss_ * poiss_);

  }

  // Explicit creation from bulk and shear modulus
  void setElasKG(double BulkModulus, double ShearModulus){

    bulkm_ = BulkModulus;
    lame2_ = ShearModulus;

    young_ = 9.0 * bulkm_ * lame2_ / (3.0 * bulkm_ + lame2_);
    poiss_ = (3.0 * bulkm_ - 2.0 * lame2_) / (2.0 * (3.0 * bulkm_ + lame2_ ));
    lame1_ = bulkm_ - 2.0 * lame2_ / 3.0;

    youngPS_ = young_ / (1.0 - poiss_ * poiss_);

  }

  // Explicit creation from the two Lame's parameters
  void setElas12(double Lame1, double Lame2){

    lame1_ = Lame1;
    lame2_ = Lame2;

    bulkm_ = lame1_ + 2.0 * lame2_ / 3.0;
    young_ = lame2_ * (3.0 * lame1_ + 2.0 * lame2_) / (lame1_ + lame2_);
    poiss_ = lame1_ / (2.0 * (lame1_ + lame2_));

    youngPS_ = young_ / (1.0 - poiss_ * poiss_);

  }

  // Copy constructor
  ElasticProperties(const ElasticProperties& theSolid){
    young_ = theSolid.young_;
    poiss_ = theSolid.poiss_;
    bulkm_ = theSolid.bulkm_;
    lame1_ = theSolid.lame1_;
    lame2_ = theSolid.lame2_;
    youngPS_= theSolid.youngPS_;
  }


  /////////// GETTER OF PARAMETERS
  // Long version
  double youngModulus(){ return young_; }
  double poissonRatio(){ return poiss_; }
  double bulkModulus(){ return bulkm_; }
  double lame1Parameter(){ return lame1_; }
  double lame2Parameter(){ return lame2_; }
  double shearModulus(){ return lame2_; }
  double planeStrainE(){ return youngPS_; }

  // Short version
  double E() const {return young_;}
  double nu() const {return poiss_;}
  double K() const {return bulkm_;}
  double L1() const {return lame1_;}
  double L2() const {return lame2_;}
  double G() const {return lame2_;}
  double Ep() const {return youngPS_;}

};
}


#endif  // HFPX2DUNITTEST_ELASTICPROPERTIES_H

