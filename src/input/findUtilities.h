//
// Created by lorenzo on 9/18/17.
//

#ifndef HFPX2DUNITTEST_FINDUTILITIES_H
#define HFPX2DUNITTEST_FINDUTILITIES_H

#include <iostream>
#include <typeinfo>
#include <il/String.h>
#include "il/base.h"
#include "il/toml.h"

namespace hfp2d {
il::String findString(const il::String &keyword,
                      const il::MapArray<il::String, il::Dynamic> &inputMap,
                      const il::String &inputFileName);

il::int_t findInteger(const il::String &keyword,
                      const il::MapArray<il::String, il::Dynamic> &inputMap,
                      const il::String &inputFileName);

double findDouble(const il::String &keyword,
                  const il::MapArray<il::String, il::Dynamic> &inputMap,
                  const il::String &inputFileName);

bool findBool(const il::String &keyword,
              const il::MapArray<il::String, il::Dynamic> &inputMap,
              const il::String &inputFileName);

}

#endif //HFPX2DUNITTEST_FINDUTILITIES_H
