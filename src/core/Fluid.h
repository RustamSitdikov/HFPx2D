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

  // could be a struct.....
//

private:
  double density_;
  double viscosity_;
  double compressibility_;

public:

   Fluid(double density, double compressibility, double viscosity) {

    density_ = density;
    viscosity_ = viscosity;
    compressibility_ = compressibility;
  }

  double fluidDensity()  const { return density_; };
  double fluidViscosity() const { return viscosity_; };
  double fluidCompressibility() const { return compressibility_; };

};

}


#endif //HFPX2DUNITTEST_FLUID_H
