//
// Created by DONG LIU on 6/2/17.
//

#ifndef HFPX2D_VISCOSITY_ALL_NODESELEMT_H
#define HFPX2D_VISCOSITY_ALL_NODESELEMT_H



#include "Stress.h"
#include "Coh_Prop_Col.h"
#include <il/linear_algebra.h>
#include <src/FVM.h>
#include <src/ConductivitiesNewtonian.h>
#include <src/FromEdgeToCol.h>
#include "Coh_Linear_softening.h"


namespace hfp2d {

    il::Array<il::int_t> stress_criteria_visco(const il::Array2D<double> &kmat,
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

    il::Array<double> conductivities_newtonian_open(const il::Array<double> &rho,
                                                    const il::Array<double> &vector,
                                                    il::Array<double> EltSizes,
                                                    Parameters_fluid &fluid_parameters,
                                                    double kf, il::io_t);

    il::Array <double> element_size(Mesh mesh_total, const int &p );

    il::Array<double> average_open(il::Array<double> width,
                                   il::Array2D<int> col_row_i,
                                   il::Array2D<int> col_row_j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim);

    il::Array<double> quarter_open(il::Array<double> width,
                                   il::Array2D<int> col_row_i,
                                   il::Array2D<int> col_row_j,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   il::Array2D<int> col_matrix);

    il::Array2D<double> matrix_lw(il::Array<double> widthB,
                                  il::Array2D<int> col_row_i,il::Array2D<int> col_row_j,
                                  const il::int_t &dof_dim, const int &p,
                                  il::Array<double> element_size_all,
                                  Parameters_fluid &fluid_parameters,
                                  const il::Array2D<int> &id);

    il::Array2D<double> matrix_vpw(Parameters_fluid &fluid_parameters,
                                   il::Array<double> element_size_all,
                                   il::Array2D<int>col_row_i,il::Array2D<int> col_row_j,
                                   il::Array<double> widthB,
                                   const il::Array2D<int> &id,
                                   const il::int_t &dof_dim,
                                   const int &p,const il::Array2D<int> &col_matrix,il::io_t);

    il::Array2D<double> matrix_vw(Parameters_fluid &fluid_parameters,
                                  il::Array<double> element_size_all,
                                  il::Array2D<int>col_row_i,il::Array2D<int> col_row_j,
                                  il::Array<double> widthB,
                                  const il::Array2D<int> &id,
                                  const il::int_t &dof_dim, const int &p,
                                  const il::Array2D<int> &col_matrix,il::io_t);

    il::int_t find_source(il::int_t i0, il::int_t j0);

    void construct_matrix_visco(il::Array2D<double> &kmatC,
                                il::Array2D<double> vwc, il::Array2D<double> vp,il::Array2D<double> ll,
                                il::Array<double> cohf, const Material &material,
                                il::Array2D<double> matrix_edge_to_col,il::Array<double> m_source,
                                const Initial_condition &initial_condition,
                                il::Status &status,
                                il::Array<double> pressure_f, il::Array<double> widthB,
                                const int &p, const il::int_t &dof_dim, il::Array2D<int> col_row_i,il::Array2D<int> col_row_j,il::io_t,
                                il::Array<double> &width, il::Array<double> &pressure,il::Array<double> &volume_vary,il::Array<double>&elastic_vary);

    void plasticity_loop_visco(const Material &material,
                               const Initial_condition &initial_condition,
                               il::int_t c_i, il::int_t c_j,
                               const il::Array2D<double> &kmat,
                               const il::Array2D<int> &id,const int &p,
                               il::Array<double> widthB, il::Array<double> pressure_f,
                               il::Array<double> width_history,
                               const il::int_t &dof_dim,
                               const il::Array2D<int>&col_matrix,
                               Parameters_fluid fluid_parameters,
                               il::Array<double> element_size_all,
                               const il::Array2D<double> &matrix_edge_to_col_all,
                               il::int_t input,il::io_t,
                               il::Array<double> &delta_width,
                               il::Array <double> &pressure_change,
                               il::Array<double> &coht,int &mm,il::Array<double> &volume_vary,
                               il::Array<double> &elastic_vary,il::Array<double>&error_w_list);
    void
    propagation_loop_visco(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                           const Material &material,
                           const Initial_condition &initial_condition,
                           il::int_t i0, il::int_t j0, int nstep,
                           il::Status &status, Parameters_fluid fluid_parameters, il::io_t,
                           il::Array2C<double> &widthlist,
                           il::Array2D<double> &plist, il::Array<double> &l_coh,
                           il::Array<double> &l, il::Array2C<double> &coh_list,
                           il::Array<int> &mvalue,int &break_time,
                           il::Array2C<double> &stress_list, il::Array<double> &energy_g,
                           il::Array2D<double> &volume_vary_list,il::Array2D<double> &elastic_vary_list,
                           il::Array2D<double> &error_matrix);

}




#endif //HFPX2D_VISCOSITY_ALL_NODESELEMT_H
