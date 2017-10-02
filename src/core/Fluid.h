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
  double fluid_density_;
  double fluid_viscosity_;
  double fluid_compressibility_;

public:

  explicit Fluid(Fluid theFluid){

    fluid_density_ = theFluid.fluid_density_;
    fluid_viscosity_ = theFluid.fluid_viscosity_;
    fluid_compressibility_ = theFluid.fluid_compressibility_;

  }


  explicit Fluid(double density, double viscosity, double compressibility) {

    fluid_density_ = density;
    fluid_viscosity_ = viscosity;
    fluid_compressibility_ = compressibility;
  }

//  void setFluidParameters(double density, double viscosity, double compressibility) {
//
//    fluid_density_ = density;
//    fluid_viscosity_ = viscosity;
//    fluid_compressibility_ = compressibility;
//  }

  double fluidDensity() { return fluid_density_; };
  double fluidViscosity() { return fluid_viscosity_; };
  double fluidCompressibility() { return fluid_compressibility_; };

};

}


#endif //HFPX2DUNITTEST_FLUID_H
