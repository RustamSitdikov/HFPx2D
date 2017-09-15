//
// Created by lorenzo on 9/13/17.
//

#ifndef HFPX2DUNITTEST_AUTOMESHUTILITIES_H
#define HFPX2DUNITTEST_AUTOMESHUTILITIES_H

#include <iostream>
#include "il/toml.h"

double findXC(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findYC(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findX1(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findY1(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findX2(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findY2(const il::String &inputFileName,
                    il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findAngle(const il::String &inputFileName,
                       il::int_t fractureID,
                 const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findLength(const il::String &inputFileName,
                        il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findNumElem(const il::String &inputFileName,
                            il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findInterpOrder(const il::String &inputFileName,
                                il::int_t fractureID,
                          const il::MapArray<il::String, il::Dynamic> &autoCreationMap)

il::String findSource(const il::String &inputFileName,
                      il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findMaterialID(const il::String &inputFileName,
                               il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findFarFieldID(const il::String &inputFileName,
                               il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findPorePresID(const il::String &inputFileName,
                               il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

#endif //HFPX2DUNITTEST_AUTOMESHUTILITIES_H
