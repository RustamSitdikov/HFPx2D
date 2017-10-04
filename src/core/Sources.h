//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2DUNITTEST_SOURCES_H
#define HFPX2DUNITTEST_SOURCES_H

#include <il/Array.h>
namespace hfp2d {
class Sources {

private:
  il::Array<double> injection_rate_;

public:
  Sources() {};

  Sources(il::Array<double> injectionRate) {

    this->injection_rate_ = injectionRate;

  };

};
}
#endif //HFPX2DUNITTEST_SOURCES_H
