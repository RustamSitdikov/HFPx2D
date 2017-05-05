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
#include "src/FVM.h"

namespace hfp2d {

    struct Material {
        double wc;
        double sigma_t;
        double Ep;
        double U;
    };

    struct Initial_condition {
        double pini;
        double Q0;
        double timestep;
        double sigma0;
    };

    void
    material_condition_col(const double &wc1, const double &sigma_t1,
                           const double &Ep1, const double &U1,
                           const double &pini1, const double &Q01,
                           const double &timestep1,
                           const double &sigma01, il::io_t, Material &material,
                           Initial_condition &initial_condition);

    il::Array<double>
    get_vwc_vc_col(Mesh mesh_total, const int &p, const il::int_t &dof_dim);


    il::Array<il::int_t> stress_criteria_col(const il::Array2D<double> &kmat,
                                             const il::Array2D<int> &id,
                                             const Material &material,
                                             const Initial_condition &initial_condition,
                                             il::Array<double> widthB,
                                             il::Array<double> delta_w,
                                             il::int_t i, il::int_t j,
                                             il::int_t c_i, il::int_t c_j,
                                             const il::int_t &dof_dim,
                                             const il::Array2D<int> &col_matrix,
                                             const int &p);


    il::Array<double>
    width_edge_col_col(il::Array<double> width, il::Array2D<int> col_row_i,
                       il::Array2D<int> col_row_j, const il::int_t &dof_dim,
                       const il::Array2D<int> &id);

    il::Array<double>
    write_history_col(il::Array<double> width, il::Array<double> width_history,
                      il::int_t c_i, il::int_t c_j, il::Array2D<int> col_matrix,
                      const Material &material, const int &p,
                      const il::int_t &dof_dim, const il::Array2D<int> &id);


    il::Array<double> cohesive_force_col(const Material &material,
                                         il::Array<double> width,
                                         il::Array2D<int> col_row_i,
                                         il::Array2D<int> col_row_j,
                                         il::Array<double> width_history,
                                         const int &p, const il::int_t &dof_dim,
                                         const il::Array2D<int> &id);

    void cc_length_col(Mesh mesh_total,
                       const Material &material, il::Array<double> width_large,
                       il::int_t i, il::int_t j, il::int_t c_i, il::int_t c_j,
                       il::Array<double> width_history, const int &p,
                       const il::Array2D<int> &id,
                       const il::Array2D<int> &col_matrix, const il::int_t &dof_dim,
                       il::io_t,
                       double &length_coh, double &crack_length);


    il::Array<double> initialwidth_col(const il::Array2D<double> &kmat,
                                       const il::Array2D<int> &id, int p,
                                       const Initial_condition &initial_condition,
                                       il::int_t i, il::int_t j,
                                       il::Status &status,
                                       const il::int_t &dof_dim);

    void construct_matrix_col(il::Array2D<double> &kmatC,
                              il::Array<double> vwc,//here not sure whether there should be a reference for kmatC
                              il::Array<double> cohf, const Material &material,
                              il::Array<double> unit,
                              const Initial_condition &initial_condition,
                              il::Status &status,
                              double pressure_f, il::Array<double> widthB,
                              const int &p, const il::int_t &dof_dim, il::io_t,
                              il::Array<double> &width,
                              double &pressure);


    void plasticity_loop_col(const Material &material,
                             const Initial_condition &initial_condition,
                             il::int_t c_i, il::int_t c_j,
                             const il::Array2D<double> &kmat, const il::Array<double> &vwc0,
                             const il::Array2D<int> &id, const int &p,
                             il::Array<double> widthB, double pressure_f,
                             il::Array<double> width_history,
                             const il::int_t &dof_dim,
                             const il::Array2D<int> &col_matrix,
                             il::io_t,
                             il::Array<double> &delta_width,
                             double &pressure_change,
                             il::Array<double> &coht, int &mm);

    il::Array2D<int> collocation_matrix(int &nelts, const int &p);

    void
    propagation_loop_col(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                         const Material &material,
                         const Initial_condition &initial_condition,
                         il::int_t i0, il::int_t j0, int nstep,
                         il::Status &status, il::io_t,
                         il::Array2C<double> &widthlist,
                         il::Array<double> &plist,
                         il::Array<double> &l_coh, il::Array<double> &l,
                         il::Array2C<double> &coh_list, il::Array<int> &mvalue,
                         int &break_time,il::Array2C<double> &stress_list);


    void get_xlist_col(il::Array<double> &xlist, Mesh mesh_total);

}



#endif //HFPX2D_COH_PROP_NEW_H
