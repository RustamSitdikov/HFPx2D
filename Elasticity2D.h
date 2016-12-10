//
// Created by Brice Lecampion on 10.12.16.
//

#ifndef HFPX2D_ELASTICITY2D_H
#define HFPX2D_ELASTICITY2D_H


#include <il/Array2D.h>
#include <il/Array.h>


il::Array2D<double> StressesKernelLinearDD(const double h,const double Ep,const double x,const double y);


void NormalShearStressKernel_LinearDD(il::Array2D<double>& St,const il::Array<double> xe, const double & h,const il::Array<double>& s, const il::Array<double>& n,
                                      const double&  Ep );


#endif //HFPX2D_ELASTICITY2D_H