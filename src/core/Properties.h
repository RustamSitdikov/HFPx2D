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
#include "Solid.h"
#include "Fluid.h"
#include "SolidEvolution.h"
#include "FluidEvolution.h"


namespace hfp2d {

// for each fracture ID we will provide a properties object which is

class Properties {

private:


public:

  Solid fault_solid;
  Fluid fault_fluid;

  //SolidEvolution
  LinearCZM evolutionCZM;
  //FluidEvolution
  ConstChannel evolutionFlow;


};


};


#endif //HFPX2D_PROPERTIES_H
