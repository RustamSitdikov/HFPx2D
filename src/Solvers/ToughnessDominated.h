//
// Created by lorenzo on 10/13/17.
//

#ifndef HFPX2D_TOUGHNESSDOMINATED_H
#define HFPX2D_TOUGHNESSDOMINATED_H

//#include <cmath>
//#include <complex>
//#include <iostream>
//
//#include <il/Array.h>
//#include <il/Array2D.h>
//#include <il/linear_algebra/dense/factorization/LU.h>
//#include <il/StaticArray.h>
//#include <il/Timer.h>
//#include <il/linear_algebra.h>
//#include <il/norm.h>
//
//
//#include "SimpleElastic.h"
//#include "src/Elasticity/Simplified3D.h"
//#include "src/Elasticity/AssemblyDDM.h"
//#include "src/Elasticity/PlaneStrainInfinite.h"
//#include "src/core/Mesh.h"
//#include "src/core/DOF_Handles.h"
//#include "src/core/ElasticProperties.h"

//#include <cmath>
//#include <complex>
#include <iostream>
//#include <string>
#include <fstream>
#include <iomanip>

#include <il/Array.h>
//#include "il/math.h"
//#include <il/Array2C.h>
//#include <il/StaticArray.h>
#include <il/Timer.h>
//#include <il/linear_algebra.h>
#include <il/norm.h>

//#include "SimpleElastic.h"
//#include "src/Elasticity/Simplified3D.h"
#include "src/Elasticity/AssemblyDDM.h"
#include "src/Elasticity/PlaneStrainInfinite.h"
#include "src/core/DOF_Handles.h"
#include "src/core/ElasticProperties.h"
#include "src/core/Mesh.h"
#include "src/core_dev/SolidEvolution.h"
#include "src/devt/FromEdgeToCol.h"
#include "src/core_dev/Simulation.h"
#include "src/core_dev/SolutionClass.h"

namespace hfp2d {

double ToughnessDominated(int nelts);

bool isElemAlreadyActive(il::Array<il::int_t> activeList, il::int_t element);

bool checkActOpening(double strShear, double strOpening);

il::StaticArray<double, 2> tractionSeparation(il::int_t i, double u_x,
                                              double u_y);

double normDD(il::Array<double> R, il::int_t dd_dofs);

double normSplit(il::Array<double> R, il::int_t begin, il::int_t end);

il::Array2D<double> takeMatFromList(il::Array2D<double> const &fullMat,
                                    il::Array<il::int_t> const &elemList,
                                    Mesh const &mesh);

il::Array<double> extractVecFromList(const il::Array<double> &fullVec,
                                     const il::Array<il::int_t> &elemList,
                                     const Mesh &mesh);

il::Array<double> extractFfromList(il::Array<double> &inSituStr,
                                   double pressure,
                                   il::Array<il::int_t> &elemList,
                                   Mesh &mesh);

il::Array2D<double> extractMatFromList(const il::Array2D<double> &fullMat,
                                       const il::Array<il::int_t> &elemList,
                                       const Mesh &mesh);

void prepareK_tough(const il::Array<il::int_t> &activeList,
                    const il::Array2D<double> &globalK_DD,
                    Mesh &mesh,
                    il::io_t,
                    il::Array2D<double> &Kact);

void prepareF_tough(const double source,
                    const SolutionK &solutionAtN,
                    const il::Array2D<double> &fetc,
                    const il::Array2D<double> &globalK_DD,
                    const il::Array<il::int_t> &activeList,
                    const Mesh &mesh,
                    const il::Array<double> &inSituStress,
                    il::io_t,
                    il::Array<double> &Fact_k,
                    SolidEvolution &CZM,
                    il::Array<double> &globalDDs);

simulationParams initSimParams(double minDeltaTime,
                               double maxDeltaTime,
                               il::int_t fracfrontMaxIter,
                               il::int_t nonlinMaxIter,
                               double tolX1,
                               double tolX2,
                               double relaxParam);

il::Array<double> insertVecInList(const il::Array<double> &smallVec,
                                  const il::Array<il::int_t> &elemList,
                                  const Mesh &mesh);

/*void initialSolution(Mesh &mesh,
                     double &initPress,
                     il::Array<double> &inSituStr,
                     il::Array2D<double> &globalK_DD,
                     il::io_t,
                     SolidEvolution &linearCZM,
                     il::Array<il::int_t> &activeList,
                     il::Array2D<double> &Kact,
                     il::Array<double> &Fact,
                     il::Array<double> &initialDDSol);*/

}

#endif // HFPX2D_TOUGHNESSDOMINATED_H
