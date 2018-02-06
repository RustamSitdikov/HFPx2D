//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2DUNITTEST_SOURCES_H
#define HFPX2DUNITTEST_SOURCES_H

// Inclusion from Inside Loop library
#include <il/Array.h>

namespace hfp2d {
class Sources {
 protected:
  // todo : move it to a struct ?

  // we start simple with one rate and one location (given by the element
  // number)
  // to modify later to 2 arrays for different sources... ? probably not needed
  // as this will be handle via the wellbore flow solver and flux partitionning
  // between fracs.
  // we are far away of injecting simultaneously in 2 boreholes!

  // then when rate is variable  will have input table of rates and time etc.

  il::Array<il::int_t> source_location_;
  il::Array<double> injection_rate_;

  //  il::Array<double> injection_rate_;

  //  il::Array

 public:
  ////////////////////////////////////////////////////////////////////////////
  // constructor
  ////////////////////////////////////////////////////////////////////////////
  Sources() {};

  Sources(il::Array<il::int_t> &sourceLoc, il::Array<double> &injectionRate) {
    // we need the corresponding location of the injection
    // element number (for P0) or nodes number for P1 ....
    //
    // can be a list if there is multiple rate.... ???

    IL_EXPECT_FAST(sourceLoc.size() == injectionRate.size());

    injection_rate_ = injectionRate;

    // in P0 we want element, and in P1 it's a node
    source_location_ = sourceLoc;
  };

  il::Array<double> InjectionRate() const { return injection_rate_; };
  double InjectionRate(il::int_t k) const { return injection_rate_[k]; };

  il::Array<il::int_t> SourceElt() const { return source_location_; };
  il::int_t SourceElt(il::int_t k) const { return source_location_[k]; };

  // just overloading here.... in case the source is at a node
  il::Array<il::int_t> SourceNodes() const { return source_location_; };
  il::int_t SourceNodes(il::int_t k) const { return source_location_[k]; };

  ////////////////////////////////////////////////////////////////////////////
  // METHODS
  ////////////////////////////////////////////////////////////////////////////

  void setInjectionRates(il::Array<double> &influx) {
    IL_EXPECT_FAST(influx.size() == source_location_.size());
    IL_EXPECT_FAST(influx.size() == injection_rate_.size());
    injection_rate_ = influx;
  };
};
}
#endif  // HFPX2DUNITTEST_SOURCES_H
