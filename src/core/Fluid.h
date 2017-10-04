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

#ifndef HFPX2DUNITTEST_FLUID_H
#define HFPX2DUNITTEST_FLUID_H

namespace hfp2d {

class Fluid {

private:
  double density_;
  double viscosity_;
  double compressibility_;

public:

  explicit Fluid(const Fluid& theFluid){

    density_ = theFluid.density_;
    viscosity_ = theFluid.viscosity_;
    compressibility_ = theFluid.compressibility_;

  }
  
  explicit Fluid(double density, double viscosity, double compressibility) {

    density_ = density;
    viscosity_ = viscosity;
    compressibility_ = compressibility;
  }

//  void setFluidParameters(double density, double viscosity, double compressibility) {
//
//    density_ = density;
//    viscosity_ = viscosity;
//    compressibility_ = compressibility;
//  }

  double fluidDensity() { return density_; };
  double fluidViscosity() { return viscosity_; };
  double fluidCompressibility() { return compressibility_; };

};

}


#endif //HFPX2DUNITTEST_FLUID_H
