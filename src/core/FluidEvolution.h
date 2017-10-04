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

#ifndef HFPX2DUNITTEST_FLUIDEVOLUTION_H
#define HFPX2DUNITTEST_FLUIDEVOLUTION_H

namespace hfp2d {

class FluidEvolution {

  // Example of constant permeability
private:

  il::String type_;

  il::Array<double> permeability_;
  il::Array<double> initial_permeability_;

public:

  explicit FluidEvolution(const FluidEvolution &theFluidEvolution){

    type_ = theFluidEvolution.type_;
    initial_permeability_ = theFluidEvolution.initial_permeability_;
    permeability_ = theFluidEvolution.permeability_;

  }

  explicit FluidEvolution(const il::Array<double> &initialPermeability) {

    type_ = il::toString("Constant permeability");
    initial_permeability_ = initialPermeability;
    permeability_ = initialPermeability;

  };

  il::String getType() { return type_; }

  double getPermeability(il::int_t k) { return initial_permeability_[k]; }

  /*void setInitialPermeability(double k) {
    permeability_ = k;
    initial_permeability_ = k;
  }*/

//  double getPermeability() {
//    return permeability_;
//  }

};

}

#endif //HFPX2DUNITTEST_FLUIDEVOLUTION_H
