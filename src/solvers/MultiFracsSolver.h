//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 15.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_MULTIFRACSSOLVER_H
#define HFPX2D_MULTIFRACSSOLVER_H

#include <src/solvers/MultiFracsSolution.h>
#include <src/util/RootFinder.h>

namespace hfp2d{




int MultipleFracsPropagation();


hfp2d::MultiFracsSolution wellHFsSolver_fixedpts(
    hfp2d::MultiFracsSolution &Sol_n,
    double dt,
    hfp2d::WellMesh &w_mesh,
    hfp2d::WellInjection &w_inj,
    hfp2d::Fluid &fracfluid,
    hfp2d::SolidProperties &rock,
    hfp2d::SimulationParameters &frac_solver_p,
    hfp2d::SimulationParameters &well_solver_p,
    hfp2d::SimulationParameters &coupling_p,
    double frac_height,
    bool mute,
    il::io_t, il::Array2D<double> &K);

// parameters structure for the interface to the root finder
struct IFParamEntryFriction : imf::IFParameters {
  double dp=0.;  // well pressure - frac entry pressure
  double fp = 0.0;  // perf. friction coef.
  double ft = 0.0; // wb tort. friction coef.
  double beta_t =0; //wb tort. friction exponent
};
/////////////////////



}


#endif //HFPX2D_MULTIFRACSSOLVER_H
