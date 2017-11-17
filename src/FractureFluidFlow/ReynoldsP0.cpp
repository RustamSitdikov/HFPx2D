//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <algorithm>

#include <il/linear_algebra/dense/norm.h>
#include "il/Array.h"

#include "ReynoldsP0.h"

#include <src/core/SolidProperties.h>
#include <src/core/Solution.h>
#include <src/core/Sources.h>
#include "src/core/Fluid.h"
#include "src/core/Mesh.h"
#include <src/core/SimulationParameters.h>


namespace hfp2d {

il::Array<double> EdgeConductivitiesP0Newtonian(
    il::Array2D<il::int_t> &edgeAdj, il::Array<double> &hydraulic_width,
    hfp2d::Fluid &fluid) {
  // input should be an array of int : row edges, column the 2 adjacent elements
  // . (and we assume that we element# = width(dof#) for simplicity
  //  hydraulic_width array with hydraulic width in each elements.

  double poiseuille_ct = 1. / (12. * fluid.fluidViscosity());
  double cr, cl;

  il::Array<double> Cond{edgeAdj.size(0), 0};
  // loop on the inner edges

  // harmonic mean   2*w_l^3*w_r^3/(w_l^3 + w_r^3)
  for (il::int_t i = 0; i < edgeAdj.size(0); i++) {
    cl = pow(hydraulic_width[edgeAdj(i, 0)], 3.);
    cr = pow(hydraulic_width[edgeAdj(i, 1)], 3.);

    Cond[i] = poiseuille_ct * 2. * (cr * cl) / (cr + cl);
  }
  return Cond;
};
//

////////////////////////////////////////////////////////////////////////////////
// FV Matrix L
// should return a sparse matrix....
// for now , for quick debug assemble a dense ;(
il::Array2D<double> BuildFD_P0(hfp2d::Mesh &mesh, hfp2d::Fluid &fluid,
                               il::Array<double> &hydraulicwidth, double coef) {
  //
  // mesh :: Mesh object
  // fluid :: Fluid properties object
  // hydraulicwidth:: array containing the hydraulic width of the different
  // cells
  // coef:: coefficient factor of all entries (typically the time-step)
  //

  // checks here....

  // loop on edges connecting 2 elements

  il::Array2D<il::int_t> edgecommon =
      mesh.getNodesSharing2Elts();  // this would be better outside the function

  il::Array2D<double> L{mesh.numberOfElts(), mesh.numberOfElts(),
                        0.};  // should be a sparse matrix..... ;(

  // compute conductivities.....
  il::Array<double> CurrentCond =
      EdgeConductivitiesP0Newtonian(edgecommon, hydraulicwidth, fluid);

  il::int_t er, el, ml;

  double hi;
  ml = edgecommon.size(0);

  for (il::int_t i = 0; i < edgecommon.size(0); i++) {
    er = edgecommon(i, 0);
    el = edgecommon(i, 1);

    hi = (mesh.eltSize(er) + mesh.eltSize(el)) / 2.;

    L(er, er) += coef * CurrentCond[i] / hi;
    L(el, el) += coef * CurrentCond[i] / hi;
    L(er, el) = -coef * CurrentCond[i] / hi;
    L(el, er) = -coef * CurrentCond[i] / hi;
  };

  return L;  // such tha L.P is Div.q in FV sense
}
////////////////////////////////////////////////////////////////////////////////
// function for current cell volume -> returning a vector
il::Array<double> P0VolumeCompressibility(hfp2d::Mesh &mesh,
                                          hfp2d::Fluid &fluid,
                                          il::Array<double> &hydraulic_width) {
  il::Array<double> volume{mesh.numberOfElts(), 0.};

  for (il::int_t e = 0; e < mesh.numberOfElts(); e++) {
    volume[e] =
        mesh.eltSize(e) * (hydraulic_width[e]) * fluid.fluidCompressibility();
  }

  return volume;
}
////////////////////////////////////////////////////////////////////////////////
Solution ReynoldsSolverP0(
    hfp2d::Solution &soln, il::Array2D<double> &ElasMat, hfp2d::Fluid &fluid,
    const hfp2d::SolidProperties &rock, const hfp2d::Sources &source, double timestep,
    bool imp_tip_width, il::Array<il::int_t> &tip_region_elt,
    il::Array<double> &tip_width, // will need to add leak-off volume...
    hfp2d::SimulationParameters &simulParams, bool mute) {
  // Solution of the Elasto-Hydrodynamics Lubrication over a time step / known
  // mesh
  // PICARD / Fixed Pt Iterations SCHEME
  // Solve for both increment of DD and fluid pressure over the time step
  //
  // P0 elements
  //
  // soln:: solution object at time tn (converged solution) [contains the mesh]
  // ElasMat:: elasticity matrix  organized in [ all normal dofs   - all shear
  // dofs ]
  // fluid    :: Fluid properties object
  // rock     :: Solid properties object
  // source   :: Source / injection object
  // timestep :: double, current size of the times step
  // imp_tip_width :: flag if width (at tn+1) in the tip region are imposed
  // tip_region_elt :: vector of element # located in the tip region
  // tip_width :: corresponding imposed tip width at tn+1
  // return a Solution at tn+1 object

  // todo: add the case of imposed tip width in tip elements as option ! ...
  // currently no inequality constraints on negative width are enforced....

  // check consistency of tip region ct.
  IL_EXPECT_FAST(tip_region_elt.size() == tip_width.size());

  hfp2d::Mesh meshn = soln.CurrentMesh();
  il::int_t n_elts = soln.CurrentMesh().numberOfElts();
  il::int_t tot_dofs = 3 * n_elts;
  //
  il::Array<double> DX_k{tot_dofs, 0.}, DX_k_1{tot_dofs, 0.},
      errDX{tot_dofs, 0.};
  il::Array<double> DW_k{n_elts, 0.};

  il::Array<double> Wn = soln.openingDD();
  il::Array<double> Vn = soln.shearDD();
  il::Array<double> Pn = soln.pressure();

  il::Array<double> sig0 = soln.sigma0();
  il::Array<double> tau0 = soln.tau0();

  il::Array<il::int_t> alldof{tot_dofs,0};
  for (il::int_t i=0;i<tot_dofs;i++){
    alldof[i]=i;
  };

  // imposed tip width increment if any (from imposed tip width at tn+1 and solution at tn)
  il::Array<double> Dw_tip{tip_width.size(),0.};
  il::Array<il::int_t> tip_w_dof{tip_width.size()};
  if (imp_tip_width) {
    for (il::int_t i=0;i<tip_width.size();i++){
      tip_w_dof[i]=meshn.dofDD(tip_region_elt[i],1); // global dof number in the DD vectors
      Dw_tip[i]=tip_width[i]-Wn[tip_region_elt[i]]; // imposed tip width increment over the time step.
    }
    };

  // reglue DDs vector at tn
  il::Array<double> DDn{2 * n_elts, 0.};
   for (il::int_t i = 0; i <  n_elts; i++) {
    DDn[2*i] = Vn[i];     // shear dd
    DDn[2*i + 1] = Wn[i];  // normal dd
  }

  // INITIALIZE THE part of the tangent matrix Xi that won't change during
  // iterations
  // initialize elasticity part of the tangent matrix.
  il::Array2D<double> Xi{tot_dofs, tot_dofs, 0.};

  for (il::int_t j = 0; j < 2 * n_elts; j++) {
    for (il::int_t i = 0; i < 2 * n_elts; i++) {
      Xi(i, j) = ElasMat(i, j);
    }
  }
  // effect of net pressure increment on elasticity eq.
  for (il::int_t j = 0; j < n_elts; j++) {
    Xi(2*j+1, j + 2*n_elts) = 1.;
  }

  //
  il::Array<double> AllCellSizes = meshn.allEltSize();

  for (il::int_t j = 0; j < n_elts; j++) {
    Xi( j+ 2*n_elts, 2*j+1) = AllCellSizes[j];
  }

  // RHS part that does not change....
  il::Array<double> Fn_elas = il::dot(ElasMat, DDn);
  il::Array<double> Gamma{tot_dofs,0.};  // tangent rhs add

   for (il::int_t i = 0; i < n_elts; i++) {
    Gamma[2*i] = tau0[i] + Fn_elas[2*i];
    Gamma[2*i + 1] = sig0[i] - Pn[i] +Fn_elas[2*i+1]  ;
  }

  il::Array2D<il::int_t> sharedEdges = meshn.getNodesSharing2Elts();

  // hydraulic width, vectors of rel-errors
  il::Array<double> Wh{n_elts, 0.}, err_Dw{n_elts, 1.}, err_Dv{n_elts, 1.},
      err_Dp{n_elts, 1.};

  il::Array<double> Residuals{tot_dofs, 1.};

  double res_norm = 0.;

  double betarela = simulParams.EHL_relaxation;
  il::int_t k = 0;

  il::int_t neq = tot_dofs - tip_width.size();
  il::Array2D<double> Xi_aux{neq,neq};
  il::Array<double> G_aux{neq};
  il::Array<il::int_t> dof_eq{neq};

  if (imp_tip_width) {
    il::int_t j2=0,j3=0,i3=0; // find dof_eq complement...
    for (il::int_t j1=0;j1<tot_dofs;j1++){
      for (j2=0;j2<tip_width.size();j2++){
        if (j1==tip_w_dof[j2]){
          break;
        }
      }
      if (j2==tip_width.size()){
        dof_eq[j3]=j1;
        j3++;
      }
    }
  };


  //-------------------------------------
  // Fixed Point Iteration Solver.
  while ((k < simulParams.EHL_max_its) &&
         (il::norm(err_Dp, il::Norm::L2) > simulParams.EHL_tolerance) &&
         (il::norm(err_Dw, il::Norm::L2) > simulParams.EHL_tolerance) &&
         (il::norm(err_Dv, il::Norm::L2) > simulParams.EHL_tolerance)) {
    k++;

    // update hydraulic width.... at tn+dt
    for (il::int_t i = 0; i < n_elts; i++) {
      Wh[i] = rock.Wh_O(0) + Wn[i] + DW_k[i];
    }

    // compute the updated Finite Diff Matrix
    il::Array2D<double> L = BuildFD_P0(meshn, fluid, Wh, timestep);
    il::Array<double> cfV = P0VolumeCompressibility(meshn, fluid, Wh);
    il::Array<double> LdotPn = il::dot(L, Pn);

    // construct the tangent system (matrix + RHS)
    // add the effect of L mat (bottom right part of Xi), the rest being
    // constant
    for (il::int_t j = 0; j < n_elts; j++) {
      for (il::int_t i = 0; i < n_elts; i++) {
        if (i == j) {
          Xi(i + 2 * n_elts, j + 2 * n_elts) = L(i, j) + cfV[i];
        } else {
          Xi(i + 2 * n_elts, j + 2 * n_elts) = L(i, j);
        }
      }
      Gamma[j + 2 * n_elts] = -LdotPn[j];
    }
    // don t forget to add sources to the RHS
    for (il::int_t i = 0; i < source.SourceElt().size(); i++) {
      Gamma[2 * n_elts + source.SourceElt(i)] += source.InjectionRate(i) * timestep;
    }

    // compute current residuals with  previous iteration solution
    // (before solution of the new iteration)
    Residuals=Gamma;
    il::blas(1.0,Xi,DX_k,-1.0,il::io,Residuals);
    res_norm = il::norm(Residuals, il::Norm::L2);

    // Solve the tangent system

    if (imp_tip_width) {
      //      std::cout <<"imposing tip width " << k <<  "\n";
      // imposing tip width increment

      // take submatrix... brute force tbc
      for (il::int_t j1=0;j1<neq;j1++){
        for (il::int_t i1=0;i1<neq;i1++){
          Xi_aux(i1,j1)=Xi(dof_eq[i1],dof_eq[j1]);
        }
      };// LHS with the contribution of the constraints.
      for (il::int_t i1=0;i1<neq;i1++) {
        G_aux[i1]=Gamma[dof_eq[i1]];
        for (il::int_t j1=0;j1<tip_width.size();j1++){
          G_aux[i1]=G_aux[i1]-Xi(dof_eq[i1],tip_w_dof[j1])*Dw_tip[j1];
        }
      }
      il::Array<double> DX_aux{neq,0.};
      il::Status status;
      // use a direct solver
      DX_aux = il::linearSolve(Xi_aux, G_aux, il::io, status);
      status.ok();
      //reglue all system
      for (il::int_t i1=0;i1<neq;i1++){
        DX_k[dof_eq[i1]]=DX_aux[i1];
      };
      for (il::int_t i1=0;i1<tip_width.size();i1++){
        DX_k[tip_w_dof[i1]]=Dw_tip[i1];
      };

    }
    else {  // no tip width imposed
      il::Status status;
      // use a direct solver
      DX_k = il::linearSolve(Xi, Gamma, il::io, status);
      status.ok();

      };

    // under-relaxation of the solution
    il::blas((1.-betarela),DX_k_1,betarela,il::io,DX_k);

    // compute relative difference between sub-its....
    for (il::int_t i = 0; i < 3 * n_elts; i++) {
      errDX[i] = abs((DX_k[i] - DX_k_1[i]) / DX_k[i]);  // note that in case
                                                        // where shear dd is
                                                        // zero, this relative
                                                        // error on the shear
                                                        // dof will be large and
                                                        // insignificant.
    }

    DX_k_1 = DX_k;  // old is new

    //  rel-error estimates splitted back in width, slip, and pressure...
    // & get back DW_k from DX
     for (il::int_t i = 0; i <  n_elts; i ++) {
      err_Dv[i] = errDX[2*i];
      err_Dw[i] = errDX[2*i + 1];
      err_Dp[i] = errDX[2*n_elts + i];
      DW_k[i] = DX_k[2*i + 1];
    }

    if (mute == false) {
      std::cout << " its " << k
                << " rel. err. dw: " << il::norm(err_Dw, il::Norm::L2) <<
                " rel. err dp: " << il::norm(err_Dp, il::Norm::L2)
                << " rel.err. Dv:" <<
                il::norm(err_Dv, il::Norm::L2) << " rel. err. Dx:" <<
                il::norm(errDX, il::Norm::L2) << "\n";
    };

  }

  if (mute == false) {
    std::cout << " end of Picard Scheme for Reynolds, after " << k << " its "
              << " rel. err. dw: " << il::norm(err_Dw, il::Norm::L2)
              << " rel. err dp: " << il::norm(err_Dp, il::Norm::L2)
              << " norm of residuals: " << res_norm << "\n";
  };

  // update total width, pressure and shear to the end values
   for (il::int_t i = 0; i < n_elts; i++) {
    Vn[i]+= DX_k[2*i];
    Wn[i] += DX_k[2*i +1 ];
    Pn[i] += DX_k[2*n_elts+i ];
  }

  return hfp2d::Solution(
      meshn, soln.time() + timestep, timestep, Wn, Vn, Pn, sig0, tau0,
      soln.front_its(), k, soln.err_front(), il::norm(err_Dv, il::Norm::L2),
      il::norm(err_Dw, il::Norm::L2), il::norm(err_Dp, il::Norm::L2));

};
}
////////////////////////////////////////////////////////////////////////////////
