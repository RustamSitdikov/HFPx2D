//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2DUNITTEST_REYNOLDSP0_H

#define HFPX2DUNITTEST_REYNOLDSP0_H

#include <il/Array.h>

#include <src/core/Fluid.h>
#include <src/core/Mesh.h>
#include <src/core/Solution.h>
#include <src/core/SimulationParameters.h>
#include <src/core/Sources.h>
#include <src/wellbore/WellMesh.h>
#include <src/wellbore/WellSolution.h>

// todo: move this namespace to Utilities.h & .cpp ?
namespace imf {

//    class ImplicitFunction {
//    protected:
//        // implicit function parameters
//        virtual struct params_ {};
//    public:
//        ImplicitFunction() {};
//        virtual double fun(double s) = 0;
//    };

//////////////////////////////////////////////////////////////////////////
//        IMPLICIT FUNCTION TOOLS
//////////////////////////////////////////////////////////////////////////

// implicit function parameters (virtual)
struct IFParameters {};

// (virtual) implicit residual function to minimize (set to zero)
typedef double (*ImplicitFunction)(double s, IFParameters &params);

//////////////////////////////////////////////////////////////////////////
//        ROUTINES (bracketing (virtual) & root finder)
//////////////////////////////////////////////////////////////////////////

// bracket structure
struct BracketS {
  double l_b;
  double u_b;
};

//    // (virtial) bracketing the root
//    typedef BracketS (*Bracketing)
//            (ImplicitFunction fun,
//             IFParameters &params,
//             double lw_bound, //
//             double up_bound, // wellMesh size
//             const int max_iter,
//             bool mute);

// Brent root finder
// template <class IFParameters>
double brent(ImplicitFunction fun, IFParameters &params, double a0, double b0,
             const double epsilon, const int max_iter, bool mute);
}

namespace hfp2d {

//////////////////////////////////////////////////////////////////////////
//        CONSTANTS  -
//////////////////////////////////////////////////////////////////////////

const double epsilon = 1e-9;
const int max_iter = 100;

//////////////////////////////////////////////////////////////////////////
//        (virtual) implicit residual function - friction factor Ff(Re_num)
//////////////////////////////////////////////////////////////////////////

// todo: expand this to tip namespace (see tip.h & .cpp)
// implicit residual function parameters for HD solvers
struct IFParametersHD : imf::IFParameters {
  double Re_num;       // Reynolds number
  double rough = 0.0;  // surface roughness (zero default)
};

//    class ImplicitFunctionHD : protected imf::ImplicitFunction {
//        friend class imf::ImplicitFunction;
//    protected:
//        struct params_ : imf::ImplicitFunction.params_ {
//            double Re_; // Reynolds number
//        };
//    public:
//        ImplicitFunctionHD() {};
//        ImplicitFunctionHD(double Re_num) {
//            struct par {
//                double Re_ = Re_num;
//            };
//            params_ = par;
//        };
//        setParams(double Re_num) {
//            params_.Re_ = Re_num;
//        };
//        Re_num() {return params_.Re_;};
//    };

// bracketing routine for HD solver
imf::BracketS bracketingHD(imf::ImplicitFunction fun, IFParametersHD &params,
                           double lw_bound,  //
                           double up_bound,  // wellMesh size
                           const int max_iter, bool mute);

// overload: bracketing the root (friction factor)
imf::BracketS bracketingHD(imf::ImplicitFunction fun, IFParametersHD &params,
                           const int max_iter, bool mute);

//////////////////////////////////////////////////////////////////////////
//        FRICTION FACTOR MODELS
//////////////////////////////////////////////////////////////////////////

// Yang & Dou model (implicit)
double imYangDou(double s, IFParametersHD &params);

// Mixed model: Laminar / Yang-Dou / Reichert
double ffYDRmixed(IFParametersHD &params);

// Mixed model: Laminar / Yang-Dou / Gauckler-Manning-Strickler
double ffYDSmixed(IFParametersHD &params);

// Churchill model (explicit)
double ffChurchill(IFParametersHD &params);

// Haaland model (explicit)
double ffHaaland(IFParametersHD &params);

// Maximum drag reduction (MDR) asymptote (implicit)
double imMDR(double s, IFParametersHD &params);

// Maximum drag reduction (MDR) approximation, explicit
double ffMDRappr(IFParametersHD &params);

// Mixed model: Laminar / MDR
double ffMDRmixed(IFParametersHD &params);
double ffMDRmixed2(IFParametersHD &params);
double ffMDRmixed3(IFParametersHD &params);

//////////////////////////////////////////////////////////////////////////
//        FLOW SOLVER ROUTINES
//////////////////////////////////////////////////////////////////////////

// nodal (edge) edge conductivities based on Darcy friction factor
// Ff(Re_num, rough)
// is this needed outsite ?
il::Array<double> edgeConductivitiesP0(
    hfp2d::WellMesh &w_mesh, il::Array<double> &velocity, hfp2d::Fluid &fluid,
    double (*ffFunction)(IFParametersHD &params));

// Finite Volume (Finite Difference) Matrix L from edge conductivities
// is this needed outsite ?
il::Array2D<double> buildWellFiniteDiffP0(hfp2d::WellMesh &w_mesh,
                                          il::Array<double> &edg_cond,
                                          double coef);

// function for current cell (element) volume increments
// is this needed outsite ?
il::Array<double> wellVolumeCompressibilityP0(
    hfp2d::Mesh &mesh, hfp2d::Fluid &fluid, il::Array<double> &hydraulic_diameter);



// Solver of the well Hydrodynamics over one time-step.
hfp2d::WellSolution wellFlowSolverP0(hfp2d::WellSolution &well_soln,
                                     hfp2d::WellMesh &w_mesh,
                                     hfp2d::WellInjection &w_inj,
                                     hfp2d::Sources &out_flow,
                                     double (*ffFunction)(IFParametersHD &),
                                     double timestep,
                                     hfp2d::SimulationParameters &simul_params,
                                     bool mute,
                                     hfp2d::Fluid &fluid);
}

#endif  // HFPX2DUNITTEST_REYNOLDSP0_H
