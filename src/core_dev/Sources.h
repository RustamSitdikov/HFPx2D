//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2DUNITTEST_SOURCES_H
#define HFPX2DUNITTEST_SOURCES_H

#include <il/Array.h>


// todo : will evolve .... we may not need a class for that.
namespace hfp2d {
class Sources {

private:
  il::Array<double> injection_rate_;

public:
  Sources() {};

  Sources(il::Array<double> injectionRate) {

    // we need the corresponding location of the injection
    // element number (for P0) or nodes number for P1 ....
    //

    // can be a list if there is multiple rate....

    this->injection_rate_ = injectionRate;


  };

};
}
#endif //HFPX2DUNITTEST_SOURCES_H
