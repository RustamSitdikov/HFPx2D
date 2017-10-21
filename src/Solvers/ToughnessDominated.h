//
// Created by lorenzo on 10/13/17.
//

#ifndef HFPX2DUNITTEST_TOUGHNESSDOMINATED_H
#define HFPX2DUNITTEST_TOUGHNESSDOMINATED_H

//#include <cmath>
//#include <complex>
//#include <iostream>
//
//#include <il/Array.h>
//#include <il/Array2D.h>
//#include <il/linear_algebra/dense/factorization/LU.h>
//#include <il/StaticArray.h>
//#include <il/Timer.h>
////#include <il/linear_algebra.h>
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

//#include <il/Array.h>
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
#include "src/core/Mesh.h"
#include "src/core/DOF_Handles.h"
#include "src/core/ElasticProperties.h"
#include "src/devt/FromEdgeToCol.h"

namespace hfp2d{

double ToughnessDominated(int nelts);

bool isElemAlreadyActive(il::Array<il::int_t> activeList, il::int_t element){

  return (std::find(activeList.begin(), activeList.end(), element) != activeList.end());

}

bool checkActOpening(double strShear, double strOpening){

  double strThreshold = 1.0e6;


  return (strOpening > strThreshold);

}

il::StaticArray<double ,2> tractionSeparation(il::int_t i,double u_x, double u_y){

  il::StaticArray<double, 2> cohForces;

  return cohForces;
};

double normDD(il::Array<double> R, il::int_t dd_dofs){

  il::Array<double> DDres(dd_dofs);
  double theNorm;

  for(il::int_t i=0; i<dd_dofs; i++){
    DDres[i] = R[i];
  }

  theNorm = il::norm(DDres,il::Norm::L2);

  return theNorm;
}

double normSplit(il::Array<double> R, il::int_t begin, il::int_t end){

  IL_EXPECT_FAST(begin > 0);
  IL_EXPECT_FAST(end > 0);
  IL_EXPECT_FAST(begin-end < R.size());

  il::int_t newSize = end-begin;

  il::Array<double> newR(newSize);
  double theNorm;

  for(il::int_t i=0; i<newSize; i++){
    newR[i] = R[begin+i];
  }

  theNorm = il::norm(newR,il::Norm::L2);

  return theNorm;
}

}

#endif //HFPX2DUNITTEST_TOUGHNESSDOMINATED_H
