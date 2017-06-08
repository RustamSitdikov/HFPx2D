//
// Created by DONG LIU on 2/23/17.
//

#ifndef HFPX2D_COH_PROP_NEW_H
#define HFPX2D_COH_PROP_NEW_H


#include "src/AssemblyDDM.h"
#include "src/DOF_Handles.h"
#include "src/Mesh.h"
#include "src/Elasticity2D.h"
#include <il/linear_algebra.h>

//#include "Coh_Propagation.h"
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

    void get_vwc_vc(il::Array<double> &vwc, Mesh mesh_total, int p);

il::Array<il::int_t> stress_criteria_new(const il::Array2D<double> &kmat,const il::Array2D<int> &id,
                                         const Material &material,const Initial_condition &initial_condition,
                                         il::Array<double> widthB,
                                         il::Array<double> delta_w, il::int_t i, il::int_t j,int p);

il::Array<double> width_edge_col(il::Array<double> width,int p);

il::Array<double> write_history(il::Array<double> width,il::Array<double> width_history,il::int_t i,const Material &material,int p);

il::Array<double> cohesive_force_new(const Material &material,
                                     il::Array<double> width, il::int_t i, il::int_t j,il::Array<double> width_history,int p);

void cc_length(Mesh mesh_total,
               const Material &material,il::Array<double> width,il::int_t i,il::int_t j,
               il::Array<double> width_history,int p,il::io_t,double &length_coh,double &crack_length);

il::Array<double> initialwidth_new(const il::Array2D<double> &kmat,
                                   const il::Array2D<int> &id,int p,
                                   const Material &material,const Initial_condition &initial_condition,
                                   il::int_t i,il::int_t j,il::Status &status);

void construct_matrix_new(il::Array2D<double> &kmatC, il::Array<double> vwc,//here not sure whether there should be a reference for kmatC
                          il::Array<double> cohf, const Material &material,
                          il::Array<double> unit, const Initial_condition &initial_condition,
                          il::Status &status, il::int_t i, il::int_t j, il::int_t i0, il::int_t j0,
                          double pressure_f, il::Array<double> widthB, int p,il::io_t,il::Array<double> &width, double &pressure);

void plasticity_loop_new(const Material &material, const Initial_condition &initial_condition,
                         il::int_t i,il::int_t j, il::int_t i0, il::int_t j0,
                         const il::Array2D<double> &kmat,Mesh mesh_total, const il::Array2D<int> &id,int p,
                         il::Array<double> widthB, double pressure_f, il::Array<double> width_history,il::io_t,il::Array<double> &delta_width,double &pressure_change,il::Array<double> &coht);

void propagation_loop_new(Mesh mesh_total,il::Array2D<int> &id,int p,
                          const Material &material,
                          const Initial_condition &initial_condition,
                          il::int_t i0, il::int_t j0, int nstep,
                          il::Status &status,il::io_t,
                          il::Array2C<double> &widthlist,
                          il::Array<double> &plist,
                          il::Array<double> &l_coh,il::Array<double> &l,il::Array2C<double> &coh_list);

    void get_xlist(il::Array<double> &xlist, Mesh mesh_total);

}




#endif //HFPX2D_COH_PROP_NEW_H
