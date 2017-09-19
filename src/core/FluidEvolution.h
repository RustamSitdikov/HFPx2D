//
// Created by lorenzo on 9/18/17.
//


//
// This file is part of HFPx2DUnitTest.
//
// Created by lorenzo on 9/18/17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2DUNITTEST_FLUIDEVOLUTION_H
#define HFPX2DUNITTEST_FLUIDEVOLUTION_H

namespace hfp2d {

class FlowEvolution {};


class ConstChannel{

  // Example of constant permeability
private:
  il::String type_;
  double permeability_;
  double initial_permeability_;

public:
  il::String getType() { return type_; }

  void setInitialPermeability(double k) {
    permeability_ = k;
    initial_permeability_ = k;
  }

  double getPermeability() {
    return permeability_;
  }

};

}

#endif //HFPX2DUNITTEST_FLUIDEVOLUTION_H
