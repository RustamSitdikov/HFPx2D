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
#include "src/core/Mesh.h"
#include "src/core/DOF_Handles.h"
#include "src/core/ElasticProperties.h"
#include "src/devt/FromEdgeToCol.h"
#include "src/core_dev/SolidEvolution.h"

namespace hfp2d{

double ToughnessDominated(int nelts);

bool isElemAlreadyActive(il::Array<il::int_t> activeList, il::int_t element);

bool checkActOpening(double strShear, double strOpening);

il::StaticArray<double ,2> tractionSeparation(il::int_t i,double u_x, double
u_y);

double normDD(il::Array<double> R, il::int_t dd_dofs);

double normSplit(il::Array<double> R, il::int_t begin, il::int_t end);

}

#endif //HFPX2D_TOUGHNESSDOMINATED_H
