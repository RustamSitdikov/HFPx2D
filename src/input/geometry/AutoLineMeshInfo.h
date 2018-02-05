//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 13.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_AUTOMESHUTILITIES_H
#define HFPX2DUNITTEST_AUTOMESHUTILITIES_H

// Inclusion from standard library
#include <iostream>

// Inclusion from the project
#include <src/input/findUtilities.h>

double findX1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findY1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findX2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

double findY2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findNumElem(const il::String &inputFileName,
                      const il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

il::int_t findMaterialID(const il::String &inputFileName,
                         const il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap);

#endif //HFPX2DUNITTEST_AUTOMESHUTILITIES_H
