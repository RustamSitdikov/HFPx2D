//
// Created by DONG LIU on 1/30/17.
//
//Propagation of the crack with the assumption of the volume control

#ifndef HFPX2D_COH_PROPAGATION_H
#define HFPX2D_COH_PROPAGATION_H


#include "Stress.h"
#include <il/linear_algebra.h>

struct Material{
    double wc;
    double sigma_t;
    double Ep;
    double U;
};

struct Initial_condition{
    double pini;
    double Q0;
    double timestep;
    double sigma0;
};

void material_condition(Material &material, Initial_condition &initial_condition,
                        double wc1, double sigma_t1,double Ep1,double U1,
                        double pini1,double Q01,double timestep1,double sigma01);

void take_submatrixint(il::Array2D<int>& sub, int i0, int i1, int j0, int j1,
                       const il::Array2D<int>& A);
Mesh take_submesh(int n1, int n2, Mesh mesh_total);

void stress_criteria(il::Array<int> &index,Mesh mesh_total,il::Array2D<int> id,
                     Material material,Initial_condition initial_condition,
                     il::Array<double> widthB,
                     il::Array<double> delta_w, int i, int j,int p);
il::Array<double> former_pressure(double pressure_f,int i,int j,int i0,int j0);

il::Array<double> former_k(il::Array2D<double> kmatC,int i,int j,int i0,int j0,il::Array<double> widthB);
il::Array<double> cohesive_force(Material material,
                                 il::Array<double> width, int i, int j, int i0, int j0);

void construct_matrix(il::Array<double> &width, double &pressure,
                      il::Array2D<double> kmatC,il::Array<double> vwc,
                      il::Array<double> cohf,Material material,
                      il::Array<double> unit,Initial_condition initial_condition,il::Status &status,int i,int j,int i0,int j0,double pressure_f,il::Array<double> widthB);
void get_vwc_vc(il::Array<double> &vwc,Mesh mesh_total,int p);
void initialwidth(il::Array<double> &delta_w_ini, Mesh mesh_total,
                  il::Array2D<int> id,int p,
                  Material material,Initial_condition initial_conditon,
                  int i,int j,il::Status &status);
void plasticity_loop(il::Array<double> &delta_width,double &pressure_change,
                     Material material, Initial_condition initial_condition, int i,
                     int j, int i0, int j0,Mesh mesh_total, il::Array2D<int> id,int p,
                     il::Array<double> widthB, double pressure_f);
void propagation_loop(il::Array2D<double> &widthlist,il::Array<double> &plist,
                      Mesh mesh_total,il::Array2D<int> id,int p,
                      Material material, Initial_condition initial_conditon,
                      int i0, int j0, int nstep,il::Status status);
void get_xlist(il::Array<double> &xlist,Mesh mesh_total);



#endif //HFPX2D_COH_PROPAGATION_H
