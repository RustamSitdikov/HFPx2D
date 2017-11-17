//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 17.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from the project
#include <src/Solvers/FluidInjFrictWeakDilatFault.h>
#include <il/linear_algebra/dense/norm.h>

#ifndef HFPX2D_REYNOLDSP1_H
#define HFPX2D_REYNOLDSP1_H

namespace hfp2d {

Solution reynoldsP1(Mesh &theMesh, il::Array2D<double> &elast_submatrix,
                    il::Array2D<double> &fetc_dd,
                    il::Array2D<double> &fetc_press, Solution &SolutionAtTn,
                    SimulationParameters &SimulationParameters,
                    SolidEvolution &SolidEvolution,
                    FractureEvolution &FractureEvolution, Sources &Source,
                    il::Array<int> &dof_active_elmnts, il::Status &status, il::Norm &norm);
}

#endif  // HFPX2D_REYNOLDSP1_H
