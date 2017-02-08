//
// Created by DONG LIU on 2/7/17.
//
//considering the couplage of fluid

#ifndef HFPX2D_VISCOSITY_H
#define HFPX2D_VISCOSITY_H


#include "Stress.h"
#include "Coh_Propagation.h"
#include <il/linear_algebra.h>

void matrix_lw(il::Array<double> widthB,  double visco,
               Mesh mesh_total,int p,int i, int j, int i0, int j0);

void matrix_edge_col(Mesh mesh_total,int i, int j, int i0, int j0,int p);

void matrix_vp(il::Array<double> widthB,Mesh mesh_total,double beta,int i,int j, int i0, int j0,int p);

void matrix_vw(Mesh mesh_total,int i,int j, int i0, int j0,int p);







#endif //HFPX2D_VISCOSITY_H
