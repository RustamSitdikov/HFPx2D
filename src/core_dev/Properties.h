//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 07.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_PROPERTIES_H
#define HFPX2D_PROPERTIES_H

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
#include "src/core/ElasticProperties.h"
#include "src/core/Fluid.h"
#include "src/FluidFlow/PermeabilityEvolution.h"
#include "src/core_dev/SolidEvolution.h"



namespace hfp2d {

// for each fracture ID we will provide a properties object which is

class Properties{ //}; : Solid, Fluid, SolidEvolution, PermeabilityEvolution {

private:

  ElasticProperties* solid_;
  Fluid* fluid_;

  SolidEvolution* solid_evolution_;
  PermeabilityEvolution* fluid_evolution_;

public:

  Properties(){};

  Properties(const Properties &theProperties){

    // pointer to pointer copy
    this->solid_ = theProperties.solid_;
    this->fluid_ = theProperties.fluid_;

    this->solid_evolution_ = theProperties.solid_evolution_;
    this->fluid_evolution_ = theProperties.fluid_evolution_;

  };

  explicit Properties(ElasticProperties &theSolid,
                      Fluid &theFluid,
                      SolidEvolution &theSolidEvolution,
                      PermeabilityEvolution &theFluidEvolution){

    *solid_=theSolid;
    *fluid_=theFluid;
    *solid_evolution_=theSolidEvolution;
    *fluid_evolution_=theFluidEvolution;

  };


};


};

#endif //HFPX2D_PROPERTIES_H
