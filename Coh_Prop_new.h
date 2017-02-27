//
// Created by DONG LIU on 2/23/17.
//

#ifndef HFPX2D_COH_PROP_NEW_H
#define HFPX2D_COH_PROP_NEW_H


#include "Stress.h"
#include <il/linear_algebra.h>
#include "Coh_Propagation.h"

void stress_criteria_new(il::Array<int> &index,il::Array2D<double> kmat,il::Array2D<int> id,
                         Material material,Initial_condition initial_condition,
                         il::Array<double> widthB,
                         il::Array<double> delta_w, int i, int j,int p);
il::Array<double> width_edge_col(il::Array<double> width);

il::Array<double> write_history(il::Array<double> width,il::Array<double> width_history,int i,Material material);

il::Array<double> cohesive_force_new(Material material,
                                     il::Array<double> width, int i, int j, int i0, int j0,il::Array<double> width_history);

void cc_length(double &length_coh,double &crack_length,Mesh mesh_total,
               Material material,il::Array<double> width,int i,int j,
               il::Array<double> width_history,int p);

void initialwidth_new(il::Array<double> &delta_w_ini, il::Array2D<double> kmat,
                      il::Array2D<int> id,int p,
                      Material material,Initial_condition initial_condition,
                      int i,int j,il::Status &status);
void plasticity_loop_new(il::Array<double> &delta_width,double &pressure_change,
                         Material material, Initial_condition initial_condition,
                         int i,int j, int i0, int j0,
                         il::Array2D<double> kmat,Mesh mesh_total, il::Array2D<int> id,int p,
                         il::Array<double> widthB, double pressure_f, il::Array<double> width_history);
void propagation_loop_new(il::Array2D<double> &widthlist,il::Array<double> &plist,
                          il::Array<double> &l_coh,il::Array<double> &l,
                          Mesh mesh_total,il::Array2D<int> id,int p,
                          Material material, Initial_condition initial_condition,
                          int i0, int j0, int nstep,il::Status status);







#endif //HFPX2D_COH_PROP_NEW_H
