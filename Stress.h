//
// Created by DONG LIU on 1/30/17.
//

#ifndef HFPX2D_STRESS_H
#define HFPX2D_STRESS_H

#include "AssemblyDDM.h"
#include "Mesh.h"
#include "DOF_Handles.h"
#include "Elasticity2D.h"

il::Array2D<double> stresscalculation(Mesh mesh,il::Array2D<int> id,int p,const double Ep,double xg, double yg);
il::Array<double> stressoutput(Mesh mesh,il::Array2D<int> id,int p,const double Ep,double xg, double yg, il::Array<double> width);


#endif //HFPX2D_STRESS_H
