//
// Created by Brice Lecampion on 10.12.16.
//

#ifndef HFPX2D_ELASTICITY2D_H
#define HFPX2D_ELASTICITY2D_H


//#include <il/Array2D.h>
#include <il/StaticArray2D.h>
//#include <il/Array.h>
#include <il/StaticArray.h>


il::StaticArray2D<double,2,4>  StressesKernelLinearDD(const double h,const double Ep,const double x,const double y);


void NormalShearStressKernel_LinearDD(il::StaticArray2D<double,2,4>& St,const il::StaticArray<double,2> xe, const double & h,const il::StaticArray<double,2> s, const il::StaticArray<double,2> n,
                                      const double  Ep );


#endif //HFPX2D_ELASTICITY2D_H