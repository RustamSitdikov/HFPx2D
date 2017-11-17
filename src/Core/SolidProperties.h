//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 11.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2DUNITTEST_ROCKPROPERTIES_H
#define HFPX2DUNITTEST_ROCKPROPERTIES_H

#include <il/Array.h>
#include <src/Core/ElasticProperties.h>

namespace hfp2d {

class SolidProperties {
  //////////////////////////////////////////////////////////////////////////
  //        CLASS MEMBERS
  //////////////////////////////////////////////////////////////////////////
 private:
  hfp2d::ElasticProperties elastic_properties_;  // elastic properties object

  il::Array<double> fracture_toughness_;  // list of fracture toughness  (should
                                          // match the number of differents
                                          // material in the matid vector of the
                                          // mesh)
  // todo :: add fracture energy array....

  il::Array<double> wh_o_;  // residual hydraulic width when the mechanical
                            // width is zero (should match the number of
                            // differents material in the matid vector of the
                            // mesh

  il::Array<double> carter_leakoff_coef_;  //  list of carter leak-off
                                           //  coefficient  (should match the
                                           //  number of differents material in
                                           //  the matid vector of the mesh)

  // add here a solid evolution object.

  //////////////////////////////////////////////////////////////////////////
  //        CONSTRUCTORS/DESTRUCTORS
  //////////////////////////////////////////////////////////////////////////
 public:
  // constructor

  SolidProperties(hfp2d::ElasticProperties &elas,
                  const il::Array<double> &toughness,
                  const il::Array<double> wh_0, const il::Array<double> &Cl) {
    elastic_properties_ = elas;
    fracture_toughness_ = toughness;
    wh_o_ = wh_0;
    carter_leakoff_coef_ = Cl;
  };

  //////////////////////////////////////////////////////////////////////////
  //        GETTER FUNCTIONS
  //////////////////////////////////////////////////////////////////////////

  // get functions
  hfp2d::ElasticProperties elasticProperties() const {
    return elastic_properties_;
  };

  il::Array<double> whO() const { return wh_o_; };

  double whO(il::int_t k) const { return wh_o_[k]; }

  il::Array<double> kIc() const { return fracture_toughness_; };

  double kIc(il::int_t k) const { return fracture_toughness_[k]; }

  il::Array<double> cl() const { return carter_leakoff_coef_; };

  double cl(il::int_t k) const { return carter_leakoff_coef_[k]; }
};
}

#endif  // HFPX2DUNITTEST_ROCKPROPERTIES_H
