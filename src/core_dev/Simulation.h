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


class Simulation {
//placeholder for simulation parameters

};

//// Structure of simulation parameters
struct simulationParams{

    double t;
    double deltat;

    double minDeltat;
    double maxDeltat;

    double errorOnDDs;
    double errorOnPress;

    double errorOnResDDs;
    double errorOnResPress;

    il::int_t ffIter;
    il::int_t ffMaxIter;

    il::int_t nlIter;
    il::int_t nlMaxIter;

};

simulationParams initSimParams(double time,
                               double deltaTime,
                               double minDeltaTime,
                               double maxDeltaTime,
                               double errDDs,
                               double errPress,
                               double errResDDs,
                               double errResPress,
                               il::int_t fracfrontIter,
                               il::int_t fracfrontMaxIter,
                               il::int_t nonlinIter,
                               il::int_t nonlinMaxIter){

    simulationParams simParams = {};

    simParams.t = time;
    simParams.deltat = deltaTime;
    simParams.minDeltat = deltaTime;
    simParams.maxDeltat = deltaTime;
    simParams.errorOnDDs = errDDs;
    simParams.errorOnPress = errPress;
    simParams.errorOnResDDs = errResDDs;
    simParams.errorOnResPress = errResPress;
    simParams.ffIter = fracfrontIter;
    simParams.ffMaxIter = fracfrontMaxIter;
    simParams.nlIter = nonlinIter;
    simParams.nlMaxIter = nonlinMaxIter;

    return simParams;
};

void saveSimParams(simulationParams &simParams,
                   il::io_t,
                   double time,
                   double deltaTime,
                   double errDDs,
                   double errPress,
                   double errResDDs,
                   double errResPress,
                   il::int_t fracfrontIter,
                   il::int_t nonlinIter){

    simParams.t = time;
    simParams.deltat = deltaTime;

    simParams.errorOnDDs = errDDs;
    simParams.errorOnPress = errPress;
    simParams.errorOnResDDs = errResDDs;
    simParams.errorOnResPress = errResPress;

    simParams.ffIter = fracfrontIter;
    simParams.nlIter = nonlinIter;

};

}
#endif //HFPX2_SIMULATION_H
