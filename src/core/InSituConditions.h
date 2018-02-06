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

#include <src/core/DomainMesh.h>
#include <src/core/Mesh.h>
#include <src/core/SegmentData.h>

namespace hfp2d {

// WE  NOW have the information of the in-situ stress field mapped on a
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

  // uniform stress field constructor
  InSituConditions(double &sxx, double &syy, double &sxy, double &p0) {
    sxx_o_ = il::Array<double>(1, sxx);
    syy_o_ = il::Array<double>(1, syy);
    sxy_o_ = il::Array<double>(1, sxy);

    sxx_grad_o_x_ = il::Array<double>(1, 0.);
    sxx_grad_o_y_ = il::Array<double>(1, 0.);

    syy_grad_o_x_ = il::Array<double>(1, 0.);
    syy_grad_o_y_ = il::Array<double>(1, 0.);

    pp_o_ = il::Array<double>(1, p0);

    pp_grad_x_ = il::Array<double>(1, 0.);
    pp_grad_y_ = il::Array<double>(1, 0.);
  };

  // we should pass the background domain mesh here
  // to have some checks !!

  // shall we have a matrix instead of individual stress components vectors ?

  InSituConditions(il::Array<double> &sxx, il::Array<double> &syy,
                   il::Array<double> &sxy, il::Array<double> &sxx_grad_x,
                   il::Array<double> &sxx_grad_y, il::Array<double> &syy_grad_x,
                   il::Array<double> &syy_grad_y,
                   hfp2d::DomainMesh &background) {
    // by definition the stresses (& grads if present) are defined at centroid,
    //
    IL_EXPECT_FAST(background.numberOfElts() == sxx.size());
    IL_EXPECT_FAST(background.numberOfElts() == syy.size());
    IL_EXPECT_FAST(background.numberOfElts() == sxy.size());
    IL_EXPECT_FAST(background.numberOfElts() == sxx_grad_x.size());
    IL_EXPECT_FAST(background.numberOfElts() == sxx_grad_y.size());
    IL_EXPECT_FAST(background.numberOfElts() == syy_grad_x.size());
    IL_EXPECT_FAST(background.numberOfElts() == syy_grad_y.size());

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
  //                  Getter functions
  //////////////////////////////////////////////////////////////////

  il::Array<double> getLocalInSituPorePressure() { return pp_o_; }
  double getLocalInSituPorePressure(il::int_t k) { return pp_o_[k]; }

  //////////////////////////////////////////////////////////////////
  //                  METHODS
  //////////////////////////////////////////////////////////////////

  // Obtain the in-situ stress at a given xy point
  il::StaticArray<double, 3> insituStressAtxy(hfp2d::DomainMesh &domain,
                                              il::Array<double> &xy) {
    il::StaticArray<double, 3> stress;

    il::int_t kelt = domain.locate(xy);
    il::Array<double> center = domain.elementCentroid(kelt);

    stress[0] = sxx_o_[kelt] + sxx_grad_o_x_[kelt] * (xy[0] - center[0]) +
                sxx_grad_o_y_[kelt] * (xy[1] - center[1]);
    stress[1] = syy_o_[kelt] + syy_grad_o_x_[kelt] * (xy[0] - center[0]) +
                syy_grad_o_y_[kelt] * (xy[1] - center[1]);
    stress[2] = sxy_o_[kelt];

    return stress;
  };

  // Get Normal and shear Traction at all collocation point of a given segment
  // case of an uniform stress field
  il::StaticArray<double, 2> uniformInsituTractions(hfp2d::SegmentData &seg) {
    IL_EXPECT_FAST(sxx_o_.size() == 1);

    // normal and tangential to the segment

    il::StaticArray<double, 2> n = seg.n();
    il::StaticArray<double, 2> s = seg.s();

    il::StaticArray<double, 2> traction;

    // s.Sig.n
    traction[0] = n[0] * s[0] * sxx_o_[0] +
                  (n[0] * s[1] + n[1] * s[0]) * sxy_o_[0] +
                  n[1] * s[1] * syy_o_[0];  // ts
    // n.Sig.n
    traction[1] = n[0] * n[0] * sxx_o_[0] + 2. * n[0] * n[1] * sxy_o_[0] +
                  n[1] * n[1] * syy_o_[0];  // tn

    return traction;
  };

  // populate all traction from an uniform stress field......
  il::Array<double> uniformAllInSituTractions(hfp2d::Mesh &mesh) {
    il::StaticArray<double, 2> aux;

    il::Array<double> all_tractions{mesh.numberDDDofs(), 0.};

    IL_EXPECT_FAST(mesh.numberDDDofs() ==
                   mesh.numberDDDofsPerElt() * mesh.numberOfElts());

    //    il::int_t p=mesh.interpolationOrder();
    il::int_t dd_e = mesh.numberDDDofsPerElt();
    for (il::int_t e = 0; e < mesh.numberOfElts(); e++) {
      hfp2d::SegmentData seg = mesh.getElementData(e);
      aux = uniformInsituTractions(seg);
      all_tractions[e * dd_e] = aux[0];
      all_tractions[e * dd_e + 1] = aux[1];
      if ((mesh.interpolationOrder() == 1)) {
        all_tractions[e * dd_e + 2] = aux[0];
        all_tractions[e * dd_e + 3] = aux[1];
      }
    }

    return all_tractions;
  }

  // populate all traction from an uniform stress field......
  il::Array<double> uniformShearInSituTractions(hfp2d::Mesh &mesh) {
    il::StaticArray<double, 2> aux;

    il::Array<double> all_shear_tractions{mesh.numberDDDofs() / 2, 0.};

    IL_EXPECT_FAST((mesh.numberDDDofs() / 2) ==
                   (mesh.numberDDDofsPerElt() / 2) * mesh.numberOfElts());

    //    il::int_t p=mesh.interpolationOrder();
    il::int_t dd_e = mesh.numberDDDofsPerElt() / 2;
    for (il::int_t e = 0; e < mesh.numberOfElts(); e++) {
      hfp2d::SegmentData seg = mesh.getElementData(e);
      aux = uniformInsituTractions(seg);
      all_shear_tractions[e * dd_e] = aux[0];
      if ((mesh.interpolationOrder() == 1)) {
        all_shear_tractions[e * dd_e +1] = aux[0];
      }
    }

    return all_shear_tractions;
  }

  // populate all traction from an uniform stress field......
  il::Array<double> uniformNormalInSituTractions(hfp2d::Mesh &mesh) {
    il::StaticArray<double, 2> aux;

    il::Array<double> all_normal_tractions{mesh.numberDDDofs() / 2, 0.};

    IL_EXPECT_FAST((mesh.numberDDDofs() / 2) ==
                   (mesh.numberDDDofsPerElt() / 2) * mesh.numberOfElts());

    //    il::int_t p=mesh.interpolationOrder();
    il::int_t dd_e = mesh.numberDDDofsPerElt() / 2;
    for (il::int_t e = 0; e < mesh.numberOfElts(); e++) {
      hfp2d::SegmentData seg = mesh.getElementData(e);
      aux = uniformInsituTractions(seg);
      all_normal_tractions[e * dd_e] = aux[1];
      if ((mesh.interpolationOrder() == 1)) {
        all_normal_tractions[e * dd_e +1] = aux[1];
      }
    }

    return all_normal_tractions;
  }

  // todo : more general cases of piece-wise gradient stress over the domain
  // mesh.
};
}

#endif  // HFPX2D_INSITUCONDITIONS_H
