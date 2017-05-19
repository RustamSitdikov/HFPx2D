//
// Created by DONG LIU on 2/23/17.
//

#ifndef HFPX2D_COH_PARTIAL_H
#define HFPX2D_COH_PARTIAL_H


#include "src/AssemblyDDM.h"
#include "src/DOF_Handles.h"
#include "src/Mesh.h"
#include "src/Elasticity2D.h"
#include <il/linear_algebra.h>
#include "src/FVM.h"
#include <Coh_Prop_Col.h>

namespace hfp2d {


    il::Array2D<double>partial_matrix(il::Array<double> coh_f);

    void construct_matrix_col_partial(il::Array2D<double> &kmatC,
                                      il::Array2D<double> vwc,
                                      il::Array<double> cohf, const Material &material,
                                      il::Array<double> unit,
                                      const Initial_condition &initial_condition,
                                      il::Status &status,
                                      double pressure_f, il::Array<double> widthB,
                                      const int &p, const il::int_t &dof_dim,
                                      il::Array2D<double> partial_m, il::io_t,
                                      il::Array<double> &width, double &pressure);

    void plasticity_loop_col_partial(const Material &material,
                                     const Initial_condition &initial_condition,
                                     il::int_t c_i, il::int_t c_j,
                                     const il::Array2D<double> &kmat,
                                     const il::Array<double> &vwc0,
                                     const il::Array2D<int> &id,const int &p,
                                     il::Array<double> widthB, double pressure_f,
                                     il::Array<double> width_history,
                                     const il::int_t &dof_dim,
                                     const il::Array2D<int>&col_matrix, il::io_t,
                                     il::Array<double> &delta_width,
                                     double &pressure_change,
                                     il::Array<double> &coht,int &mm);

    void
    propagation_loop_col_partial(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                                 const Material &material,
                                 const Initial_condition &initial_condition,
                                 il::int_t i0, il::int_t j0, int nstep,
                                 il::Status &status, il::io_t,
                                 il::Array2C<double> &widthlist,
                                 il::Array<double> &plist, il::Array<double> &l_coh,
                                 il::Array<double> &l, il::Array2C<double> &coh_list,
                                 il::Array<int> &mvalue,int &break_time,
                                 il::Array2C<double> &stress_list,
                                 il::Array<double> &energy_g);



}



#endif //HFPX2D_COH_PROP_NEW_H
