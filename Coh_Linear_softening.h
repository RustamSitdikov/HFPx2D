//
// Created by DONG LIU on 5/24/17.
//

#ifndef HFPX2D_COH_LINEAR_SOFTENING_H
#define HFPX2D_COH_LINEAR_SOFTENING_H

#include "Coh_Prop_Col.h"

namespace hfp2d {

    il::Array<double> cohesive_force_linear(const Material &material,
                                         il::Array<double> width,
                                         il::Array2D<int> col_row_i,
                                         il::Array2D<int> col_row_j,
                                         il::Array<double> width_history,
                                         const int &p, const il::int_t &dof_dim,
                                         const il::Array2D<int> &id);

    void plasticity_loop_linear(const Material &material,
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

    void
    propagation_loop_linear(Mesh mesh_total, il::Array2D<int> &id, const int &p,
                         const Material &material,
                         const Initial_condition &initial_condition,
                         il::int_t i0, il::int_t j0, int nstep,
                         il::Status &status, il::io_t,
                         il::Array2C<double> &widthlist,
                         il::Array<double> &plist,
                         il::Array<double> &l_coh, il::Array<double> &l,
                         il::Array2C<double> &coh_list, il::Array<int> &mvalue,
                         int &break_time,il::Array2C<double> &stress_list,il::Array<double> &energy_g);


}





#endif //HFPX2D_COH_LINEAR_SOFTENING_H
