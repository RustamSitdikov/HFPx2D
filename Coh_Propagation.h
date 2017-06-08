//
// Created by DONG LIU on 1/30/17.
//
//Propagation of the crack with the assumption of the volume control

#ifndef HFPX2D_COH_PROPAGATION_H
#define HFPX2D_COH_PROPAGATION_H


#include "Stress.h"
#include <il/linear_algebra.h>

namespace hfp2d{

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



void
material_condition(const double &wc1, const double &sigma_t1, const double &Ep1, const double &U1,
                   const double &pini1, const double &Q01, const double &timestep1,
                   const double &sigma01,il::io_t,Material &material, Initial_condition &initial_condition);

//    void
//    take_submatrixint(il::Array2D<int> &sub, il::int_t i0, il::int_t i1, il::int_t j0, il::int_t j1,
//                      const il::Array2D<int> &A);

//    Mesh take_submesh(il::int_t n1, il::int_t n2, Mesh mesh_total);

//    void
//    stress_criteria(il::Array<int> &index, Mesh mesh_total, il::Array2D<int> id,
//                    Material material, Initial_condition initial_condition,
//                    il::Array<double> widthB,
//                    il::Array<double> delta_w, int i, int j, int p);

//    il::Array<double>
//    former_pressure(double pressure_f, il::int_t i, il::int_t j, il::int_t i0, il::int_t j0, int p);
//
//il::Array<double> former_k(il::Array2D<double> kmatC, il::int_t i, il::int_t j,
//                           il::int_t i0, il::int_t j0, il::Array<double> widthB,
//                           int p) ;

//    il::Array<double> cohesive_force(Material material,
//                                     il::Array<double> width, il::int_t i, il::int_t j,
//                                     il::int_t i0, il::int_t j0);

//void construct_matrix(il::Array2D<double> &kmatC, il::Array<double> vwc,
//                      il::Array<double> cohf, const Material &material,
//                      il::Array<double> unit, const Initial_condition &initial_condition,
//                      il::Status &status, il::int_t i, il::int_t j, il::int_t i0, il::int_t j0,
//                     double pressure_f, il::Array<double> widthB, int p,il::io_t,il::Array<double> &width, double &pressure);

    void get_vwc_vc(il::Array<double> &vwc, Mesh mesh_total, int p);

//    void initialwidth(il::Array<double> &delta_w_ini, Mesh mesh_total,
//                      il::Array2D<int> id, int p,
//                      Material material, Initial_condition initial_conditon,
//                      il::int_t i, il::int_t j, il::Status &status);

//    void
//    plasticity_loop(il::Array<double> &delta_width, double &pressure_change,
//                    Material material, Initial_condition initial_condition,
//                    il::int_t i,
//                    il::int_t j, il::int_t i0, il::int_t j0, Mesh mesh_total, il::Array2D<int> id,
//                    int p,
//                    il::Array<double> widthB, double pressure_f);

//    void
//    propagation_loop(il::Array2D<double> &widthlist, il::Array<double> &plist,
//                     Mesh mesh_total, il::Array2D<int> id, int p,
//                     Material material, Initial_condition initial_conditon,
//                     il::int_t i0, il::int_t j0, il::int_t nstep, il::Status status);

    void get_xlist(il::Array<double> &xlist, Mesh mesh_total);

}

#endif //HFPX2D_COH_PROPAGATION_H
