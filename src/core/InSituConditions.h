//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 13.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

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
  // switch to Array of double
  //
  il::Array<double> pp_o_;
  il::Array<double> pp_grad_x_;
  il::Array<double> pp_grad_y_;

  il::Array<double> sxx_o_;
  il::Array<double> syy_o_;
  il::Array<double> sxy_o_;

  il::Array<double> stress_grad_x_;  // for syy  (such that Div Sig )
  il::Array<double> stress_grad_y_;  // for sxx

 public:
  //////////////////////////////////////////////////////////////////////////////
  // constructor
  // we should pass the background domain mesh here
  // to have some checks

  InSituConditions(il::Array<double> &sxx, il::Array<double> &syy,
                   il::Array<double> &sxy, il::Array<double> &stress_grad_x,
                   il::Array<double> &stress_grad_y) {
    sxx_o_ = sxx;
    syy_o_ = syy;
    sxy_o_ = sxy;
    stress_grad_x = stress_grad_x;
    stress_grad_y_ = stress_grad_y;
  };
  // todo constructor with pp_

  // Obtain the in-situ stress at a given xy point

  il::StaticArray<double, 3> insituStressAtxy(hfp2d::DomainMesh &domain,
                                              il::Array<double> &xy) {
    il::StaticArray<double, 3> stress;

    il::int_t kelt = domain.locate(xy);


    return stress;
  };


};

}

#endif  // HFPX2D_INSITUCONDITIONS_H
