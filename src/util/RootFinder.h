//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 18.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_ROOTFINDER_H
#define HFPX2D_ROOTFINDER_H


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



#endif //HFPX2D_ROOTFINDER_H
