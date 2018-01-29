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
#include <src/core/ElasticProperties.h>
#include <src/core/FluidProperties.h>
#include <src/core/Mesh.h>
#include <src/core/SimulationParameters.h>
#include <src/core/Solution.h>
#include <src/core/Sources.h>
#include <src/core_dev/FractureEvolution.h>
#include <src/core_dev/SolidEvolution.h>

namespace hfp2d {

void fluidInjFrictWeakDilatFault(int argc, char const *argv[]);

Solution fractFrontPosition(
    il::Array2D<double> &elast_matrix, il::Array2D<double> &fetc_dds,
    il::Array2D<double> &fetc_dd, il::Array2D<double> &fetc_press,
    il::Array2D<il::int_t> &dof_single_dd, Mesh &theMesh,
    FluidProperties &FluidProperties,
    SimulationParameters &SimulationParameters, SolidEvolution &SolidEvolution,
    FractureEvolution &FractureEvolution, Sources &Source,
    Solution &SolutionAtTn, bool expl_impl, bool damping_term,
    double damping_coeff, double dilat_plast);
}

#endif  // HFPX2D_FLUIDINJFRICTWEAKDILATFAULT_H
