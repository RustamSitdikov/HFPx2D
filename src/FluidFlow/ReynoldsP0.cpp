//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "ReynoldsP0.h"
#include <src/core/RockProperties.h>
#include <src/core/SolutionAtT.h>
#include <src/core_dev/Sources.h>

#include <il/linear_algebra/dense/norm.h>
#include "il/Array.h"

#include "src/core/Fluid.h"
#include "src/core/Mesh.h"

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

il::Array2D<il::int_t> GetEdgesSharing2(hfp2d::Mesh &mesh) {
  // get element sharing common edges

  // loop over all nodes of the mesh
  // find corresponding elements
  // if n of elt sharing that node ==2 -> put in array...
  //

  // maximum number of vertex with 2 elements is
  // this will not work for the case of more than 2 elements sharing the node.
  // we don t care of that case for now.

  // format is col1: element1, col2: element2   (note we don t store the
  // corresponding node here....)

  il::Array2D<il::int_t> edge{mesh.numberOfNodes(),2, 0.};

  il::StaticArray<il::int_t, 2> temp;

  int j = 0;
  int k = 0;

  for (il::int_t i = 0; i < mesh.numberOfNodes(); i++) {
    j = 0;
    //    temp[j]=i;j++;
    for (il::int_t e = 0; e < mesh.numberOfElements(); e++) {
      if (mesh.connectivity(e, 0) == i) {
        temp[j] = e;
        j++;
      } else if (mesh.connectivity(e, 1) == i) {
        temp[j] = e;
        j++;
      };

      if (j == 2) {
        edge(k, 0) = temp[0];
        edge(k, 1) = temp[1];
        //         edge(k, 2) = temp[2];
        //         temp[0]=0;temp(0,1)=0;
        j = 0;
        k++;
        break;
      };
    }

    j = 0;
  }

  il::Array2D<il::int_t> outputedge{k, 2, 0.};

  for (il::int_t i = 0; i < outputedge.size(0); i++) {
    outputedge(i, 0) = edge(i, 0);
    outputedge(i, 1) = edge(i, 1);
    //    outputedge(i,2)=edge(i,2);
  }

  return outputedge;
}
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

  // checks here....

  // loop on edges.....
  il::Array2D<il::int_t> edgecommon = hfp2d::GetEdgesSharing2(mesh);

  il::Array2D<double> L{mesh.numberOfElements(), mesh.numberOfElements(),
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

    hi = (mesh.eltsize(er) + mesh.eltsize(el)) / 2.;

    L(er, er) += coef * CurrentCond[i] / hi;
    L(el, el) += coef * CurrentCond[i] / hi;
    L(er, el) = -coef * CurrentCond[i] / hi;
    L(el, er) = -coef * CurrentCond[i] / hi;
  };

  return L;  // such tha L.P is Div.q in FV sense
}
////////////////////////////////////////////////////////////////////////////////
// function for current cell volume -> returning a vector
//
il::Array<double> P0VolumeCompressibility(hfp2d::Mesh &mesh,
                                          hfp2d::Fluid &fluid,
                                          il::Array<double> &hydraulic_width) {
  il::Array<double> volume{mesh.numberOfElements(), 0.};

  for (il::int_t e = 0; e < mesh.numberOfElements(); e++) {
    volume[e] =
        mesh.eltsize(e) * (hydraulic_width[e]) * fluid.fluidCompressibility();
  }

  return volume;
}

////////////////////////////////////////////////////////////////////////////////
// function for elt size  -> returning a vector -> could be moved as method of
// the mesh class ?
il::Array<double> EltSize(hfp2d::Mesh &mesh) {
  il::Array<double> all_eltsize{mesh.numberOfElements(), 0.};

  for (il::int_t e = 0; e < mesh.numberOfElements(); e++) {
    all_eltsize[e] = mesh.eltsize(e);
  }
  return all_eltsize;
}
////////////////////////////////////////////////////////////////////////////////

hfp2d::SolutionAtT ReynoldsSolverP0(hfp2d::SolutionAtT &soln,
                                    il::Array2D<double> &ElasMat,
                                    hfp2d::Fluid &fluid,
                                    hfp2d::RockProperties &rock,
                                    hfp2d::Sources &source, double timestep) {
  // Solution of the Elasto-Hydrodynamics Lubrication over a time step / known
  // mesh
  // PICARD/ Fixed Pt Iterations SCHEME
  // Solve for both increment of DD and fluid pressure over the time step
  //
  // P0 elements
  //
  // soln:: solution object at time tn (converged solution) [contains the mesh]
  // ElasMat:: elasticity matrix  organized in [ all normal dofs   - all shear
  // dofs ]
  // fluid :: Fluid properties object
  // source :: Sourcce / injection object
  // timestep:: double, current size of the times step
  // return a Solution at tn+1 object

  // todo: add simulation parameters data as input
  // todo: add the case of imposed tip width in tip elements as option...
  // currently no inequality constraints on negative width are enforced....

  hfp2d::Mesh meshn = soln.CurrentMesh();

  il::int_t n_elts = soln.CurrentMesh().numberOfElements();
  il::int_t tot_dofs = 3 * n_elts;
  //
  il::Array<double> DX_k{tot_dofs, 0.}, DX_k_1{tot_dofs, 0.},
      errDX{tot_dofs, 0.};
  il::Array<double> DW_k{n_elts, 0.};

  il::Array<double> Wn = soln.openingDD();
  il::Array<double> shn = soln.shearDD();
  il::Array<double> Pn = soln.pressure();

  il::Array<double> sig0 = soln.sigma0();
  il::Array<double> tau0 = soln.tau0();


  // reglue DDs vector.
  il::Array<double> DDn{2 * n_elts, 0.};
  for (il::int_t i = 0; i < n_elts; i++) {
    DDn[i] = shn[i];  // shear dd
    DDn[i + n_elts] = Wn[i]; // normal dd
  }

  // INITIALIZE THE part of the tangent matrix that won't change during
  // iterations
  // initialize elasticity part of the tangent matrix.
  il::Array2D<double> Xi{tot_dofs, tot_dofs, 0.};

  for (il::int_t j = 0; j < 2 * n_elts; j++) {
    for (il::int_t i = 0; i < 2 * n_elts; i++) {
      Xi(i, j) = ElasMat(i, j);
    }
  }
  // effect of net pressure increment on elasticity eq.
  for (il::int_t j = n_elts; j < 2 * n_elts; j++) {
    Xi(j, j + n_elts) = 1.;
  }

  //
  il::Array<double> CellSize = EltSize(meshn);
  for (il::int_t j = 0; j < n_elts; j++) {
    Xi(2 * n_elts + j, j) = CellSize[j];
  }

  // RHS
  il::Array<double> Fn_elas = il::dot(ElasMat, DDn);
  il::Array<double> Gamma{tot_dofs, 0.};  // tangent rhs
  for (il::int_t i = 0; i < n_elts; i++) {
    Gamma[i] = tau0[i] - Fn_elas[i];
    Gamma[i + n_elts] = sig0[i] - Pn[i] - Fn_elas[i + n_elts];
  }

  il::Array2D<il::int_t> sharedEdges = GetEdgesSharing2(meshn);//

  // hydraulic width;
  il::Array<double> Wh{n_elts, 0.}, err_Dw{n_elts, 1.}, err_Dv{n_elts, 1.},
      err_Dp{n_elts, 1.};

  il::Array<double> Residuals{tot_dofs,1.};
  double res_norm=0.;

  double tolerance = 1.e-6;
  double betarela = 0.7;
  il::int_t k = 0, itermax = 100;

  // Fixed Point Iteration Solver.
  while ((k < itermax) && (il::norm(err_Dp, il::Norm::L2) > tolerance) &&
         (il::norm(err_Dw, il::Norm::L2) > tolerance) &&
         (il::norm(err_Dv, il::Norm::L2) > tolerance)) {
    k++;

    // update hydraulic width....
    for (il::int_t i = 0; i < n_elts; i++) {
      Wh[i] = rock.Wh_O(0) + Wn[i] + DW_k[i];
    }

    // compute the updated Finite Diff Matrix
    il::Array2D<double> L = BuildFD_P0(meshn, fluid, Wh, timestep);
    il::Array<double> cfV = P0VolumeCompressibility(meshn, fluid, Wh);
    il::Array<double> LdotPn = il::dot(L, Pn);

    // construct the tangent system (matrix + RHS)
    // add the effect of L mat (bottom right part of Xi), the rest being constant
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
    for (il::int_t i=0;i<source.SourceElt().size();i++) {
      Gamma[2 * n_elts + source.SourceElt(i)] += source.InjectionRate(i) * timestep;
    }

    // compute current residuals.. (before solution of the system)
    Residuals=il::dot(Xi,DX_k);
    for (il::int_t i=0;i<tot_dofs;i++){
      Residuals[i]=Residuals[i]-Gamma[i];
    }
    res_norm=il::norm(Residuals,il::Norm::L2);

    // Solve the tangent system
    il::Status status;
    // use a direct solver
    DX_k = il::linearSolve(Xi, Gamma, il::io, status);
    status.ok();

    // under-relaxation of the solution
    for (il::int_t i = 0; i < 3 * n_elts; i++) {
      DX_k[i] = betarela * DX_k[i] + (1. - betarela) * DX_k_1[i];
      errDX[i] = abs((DX_k[i] - DX_k_1[i]) / DX_k[i]); // note that in case where shear dd is zero, this relative error on the shear dof will be large and insignificant.
      DX_k_1[i] = DX_k[i];  // old is new
    }

    // compute rel-error estimates ...
    // & get back DW_k from DX
    for (il::int_t i = 0; i < n_elts; i++) {
      err_Dv[i] = errDX[i];
      err_Dw[i] = errDX[i + n_elts];
      err_Dp[i] = errDX[i + 2*n_elts];
      DW_k[i] = DX_k[i + n_elts];
    }

//    std::cout << " its " << k
//              <<    " rel. err. dw: "<< il::norm(err_Dw,il::Norm::L2) <<  " rel. err dp: " << il::norm(err_Dp,il::Norm::L2)<< " rel. err. Dv:" <<
//             il::norm(err_Dv,il::Norm::L2) << " rel. err. Dx:" << il::norm(errDX,il::Norm::L2)  << "\n" ;

  }

  std::cout << " end of Picard Scheme for Reynolds, after " << k
            << " its " <<  " rel. err. dw: "<< il::norm(err_Dw,il::Norm::L2) <<  " rel. err dp: " << il::norm(err_Dp,il::Norm::L2)
            <<  " norm of residuals: " << res_norm
            <<                                                                  "\n" ;

  // update total width, pressure and shear to the end values
  for (il::int_t i=0;i<n_elts;i++){
    Wn[i]+=DX_k[i + n_elts];
    Pn[i]+=DX_k[i + 2*n_elts];
    shn[i]+=DX_k[i];
  }

  return hfp2d::SolutionAtT(meshn,soln.time()+timestep,timestep,Wn,shn,Pn,sig0,tau0,soln.front_its(),k,soln.err_front(),
                            il::norm(err_Dv,il::Norm::L2) ,il::norm(err_Dw,il::Norm::L2) ,il::norm(err_Dp,il::Norm::L2)  );

  //   will need to create a deep copy of soln at some point
};
}
////////////////////////////////////////////////////////////////////////////////
