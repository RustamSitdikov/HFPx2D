//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 07.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_SIMULATION_H
#define HFPX2D_SIMULATION_H

namespace hfp2d {


//class Simulation {
////placeholder for simulation parameters
//};

//// Structure of simulation parameters
struct simulationParams{

    double minDeltat;
    double maxDeltat;

    il::int_t ffMaxIter;
    il::int_t nlMaxIter;

    double tolX1;
    double tolX2;

    double relaxParam;

};

}
#endif //HFPX2_SIMULATION_H
