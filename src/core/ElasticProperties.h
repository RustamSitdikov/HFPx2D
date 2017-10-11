//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_ELASTICPROPERTIES_H
#define HFPX2D_ELASTICPROPERTIES_H

namespace hfp2d {

class ElasticProperties {
  // class to store Elastic properties of the medium - only isotropy for now
  // should generalized...
 private:
  double Young_; // Young;s modulus
  double nu_; // Poisson's ratio
  double G_; // shear modulus
  double K_; // Bulk modulus
  double Ep_; // Plane-strain modulus

 public:

  ElasticProperties(){};

// we construct from value of Young and PR
  ElasticProperties(double Young, double nu) {

    Young_ = Young;
    nu_ = nu;
    G_ = Young / (2. * (1 + nu));
    Ep_ = Young / (1. - nu * nu);   // plane strain Ep
    K_ = Young / (3. * (1. - 2. * nu));
  }

  double Young() const {return Young_;}
  double nu() const {return nu_;}
  double G() const {return G_;}
  double Ep() const {return Ep_;}
  double K() const {return K_;}

// should set other
// should set other

};
}

#endif  // HFPX2DUNITTEST_ELASTICPROPERTIES_H

