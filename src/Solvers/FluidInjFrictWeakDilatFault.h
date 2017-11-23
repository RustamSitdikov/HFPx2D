//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 08.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_FLUIDINJFRICTWEAKDILATFAULT_H
#define HFPX2D_FLUIDINJFRICTWEAKDILATFAULT_H

// Inclusion from the project
#include <src/Core/ElasticProperties.h>
#include <src/Core/FluidProperties.h>
#include <src/Core/Mesh.h>
#include <src/Core/SimulationParameters.h>
#include <src/Core/Solution.h>
#include <src/Core/Sources.h>
#include <src/Core_dev/FractureEvolution.h>
#include <src/Core_dev/SolidEvolution.h>

namespace hfp2d {

void fluidInjFrictWeakDilatFault(int argc, char const *argv[]);

Solution fractFrontPosition(il::Array2D<double> &elast_matrix,
                        il::Array2D<double> &fetc_dds,
                        il::Array2D<double> &fetc_dd,
                        il::Array2D<double> &fetc_press, Mesh &theMesh,
                        FluidProperties &FluidProperties,
                        SimulationParameters &SimulationParameters,
                        SolidEvolution &SolidEvolution,
                        FractureEvolution &FractureEvolution, Sources &Source,
                        Solution &SolutionAtTn);
}

#endif  // HFPX2D_FLUIDINJFRICTWEAKDILATFAULT_H
