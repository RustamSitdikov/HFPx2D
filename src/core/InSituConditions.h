//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 13.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
// refactor to use background mesh Brice Lecampion 13.12.17.

#ifndef HFPX2D_INSITUCONDITIONS_H
#define HFPX2D_INSITUCONDITIONS_H

#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/math.h>

#include <src/core/DomainMesh.h>

namespace hfp2d {

// WE SHOULD NOW have the information of the in-situ stress field mapped on a
// background mesh

// we store the basic model parameters of a
class InSituConditions {
 private:
  //  stress at center of elements of the background mesh
  //   sxx = sxx_center + D_x_sxx (x-xcenter) +  D_y_sxx (y-ycenter)
  //   syy = syy_center + D_x_syy (x-xcenter) +  D_y_syy (y-ycenter)
  //   sxy =  sxy_center + NO gradient for now (because we are lazy)
  // Note that we don t check if  we  have Div Sig_0 + b = 0 !
  //

  il::Array<double> pp_o_;
  il::Array<double> pp_grad_x_;
  il::Array<double> pp_grad_y_;

  il::Array<double> sxx_o_;
  il::Array<double> syy_o_;
  il::Array<double> sxy_o_;

  il::Array<double> syy_grad_o_x_;  // for syy
  il::Array<double> syy_grad_o_y_;  // for syy

  il::Array<double> sxx_grad_o_y_;  // for sxx
  il::Array<double> sxx_grad_o_x_;  // for sxx

 public:
  //////////////////////////////////////////////////////////////////////////////
  // constructor
  // we should pass the background domain mesh here
  // to have some checks !!

  // shall we have a matrix instead of individual stress components vectors ?

  InSituConditions(il::Array<double> &sxx, il::Array<double> &syy,
                   il::Array<double> &sxy,
                   il::Array<double> &sxx_grad_x,
                   il::Array<double> &sxx_grad_y,
                   il::Array<double> &syy_grad_x,
                   il::Array<double> &syy_grad_y,
                   hfp2d::DomainMesh &background) {
    // by definition the stresses (& grads if present) are defined at centroid,
    //
    IL_EXPECT_FAST(background.numberOfElts() == sxx.size());
    IL_EXPECT_FAST(background.numberOfElts() == syy.size());
    IL_EXPECT_FAST(background.numberOfElts() == sxy.size());
    IL_EXPECT_FAST(background.numberOfElts() == stress_grad_x.size());
    IL_EXPECT_FAST(background.numberOfElts() == stress_grad_y.size());

    sxx_o_ = sxx;
    syy_o_ = syy;
    sxy_o_ = sxy;
    sxx_grad_o_x_ = sxx_grad_x;
    sxx_grad_o_y_ = sxx_grad_y;

    syy_grad_o_x_ = syy_grad_x;
    syy_grad_o_y_ = syy_grad_y;
  };
  // todo constructor with pp_ !

  //////////////////////////////////////////////////////////////////
  //                  METHODS
  //////////////////////////////////////////////////////////////////

  // Obtain the in-situ stress at a given xy point
  il::StaticArray<double, 3> insituStressAtxy(hfp2d::DomainMesh &domain,
                                              il::Array<double> &xy) {

    il::StaticArray<double, 3> stress;

    il::int_t kelt = domain.locate(xy);
    il::Array<double> center = domain.elementCentroid(kelt);

    stress[0] = sxx_o_[kelt] + sxx_grad_o_x_[kelt]* (xy[0] - center[0])
        + sxx_grad_o_y_[kelt] * (xy[1] - center[1]);
    stress[1] = syy_o_[kelt] + syy_grad_o_x_[kelt] * (xy[0] - center[0])
        + syy_grad_o_y_[kelt] * (xy[1] - center[1]);
    stress[2] = sxy_o_[kelt];

    return stress;
  };


  // Get Normal and shear Traction at all collocation point of a given segment ?



  // todo :

};
}

#endif  // HFPX2D_INSITUCONDITIONS_H
