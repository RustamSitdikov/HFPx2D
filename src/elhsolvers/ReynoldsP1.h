//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 17.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/linear_algebra/dense/norm.h>

// Inclusion from the project
#include <src/solvers/FluidInjFrictWeakDilatFault.h>

#ifndef HFPX2D_REYNOLDSP1_H
#define HFPX2D_REYNOLDSP1_H

namespace hfp2d {

Solution reynoldsP1(Mesh &theMesh, il::Array2D<double> &elast_matrix,
                    il::Array2D<double> &fetc_dds, il::Array2D<double> &fetc_dd,
                    il::Array2D<double> &fetc_press, Solution &SolutionAtTn,
                    Solution &SolutionAtTn_k, il::Array<double> &incrm_shearDD,
                    il::Array<double> &incrm_openingDD,
                    SimulationParameters &SimulationParameters,
                    FluidProperties &FluidProperties,
                    SolidEvolution &SolidEvolution,
                    FractureEvolution &FractureEvolution, Sources &Source,
                    il::Array<int> &dof_active_elmnts, il::Status &status,
                    il::Norm &norm, bool damping_term, double damping_coeff,
                    double dilat_plast,
                    hfp2d::InSituStress &BackgroundLoadingConditions);
}

#endif  // HFPX2D_REYNOLDSP1_H
