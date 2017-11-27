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

#ifndef HFPX2DUNITTEST_FLUID_H
#define HFPX2DUNITTEST_FLUID_H

namespace hfp2d {

class FluidProperties {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////

 private:
  double density_;
  double viscosity_;
  double compressibility_;

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////

 public:
  FluidProperties() = default;

  FluidProperties(double density, double viscosity, double compressibility) {
    density_ = density;
    viscosity_ = viscosity;
    compressibility_ = compressibility;
  }

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

 public:
  double fluidDensity() { return density_; };
  double fluidViscosity() { return viscosity_; };
  double fluidCompressibility() { return compressibility_; };

  //////////////////////////////////////////////////////////////////////////
  //        METHODS
  //////////////////////////////////////////////////////////////////////////

  void setFluidProp(double dens, double visc, double compress) {
    density_ = dens;
    viscosity_ = visc;
    compressibility_ = compress;
  }
};
}

#endif  // HFPX2DUNITTEST_FLUID_H
