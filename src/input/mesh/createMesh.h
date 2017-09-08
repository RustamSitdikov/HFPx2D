//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 06.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2DUNITTEST_CREATEMESH_H
#define HFPX2DUNITTEST_CREATEMESH_H


#include <il/base.h>
#include <il/Array.h>
#include <il/Array2D.h>
#include <cmath>
#include "src/core/Mesh.h"

namespace hfp2d {

void createVerticalMesh(double& x_c,
                        double& y_c,
                        double& length,
                        il::int_t& numElements,
                        il::int_t& materialID,
                        il::io_t,
                        Mesh& theMesh);

void createHorizontalMesh(double& x_c,
                          double& y_c,
                          double& length,
                          il::int_t& numElements,
                          il::int_t& materialID,
                          il::io_t,
                          Mesh& theMesh);

void createDiagonalMesh(double& x_c,
                        double& y_c,
                        double& angle,
                        double& length,
                        il::int_t& numElements,
                        il::int_t& materialID,
                        il::io_t,
                        Mesh& theMesh);

void createCustomMesh(double& x_1,
                      double& y_1,
                      double& x_2,
                      double& y_2,
                      il::int_t& numElements,
                      il::int_t& materialID,
                      il::io_t,
                      Mesh& theMesh);

}
#endif //HFPX2DUNITTEST_CREATEMESH_H
