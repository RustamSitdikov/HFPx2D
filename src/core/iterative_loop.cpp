//
// Created by lorenzo on 9/20/17.
//

#include "iterative_loop.h"

namespace hfp2d {

Solution iterative_solution(Mesh &theMesh,
                            Properties &theProperties,
                            Source &theSource,
                            Solution &theSolutionAtN,
                            Simulation &SimParam) {

  Solution theSolutionAtN1;

  for (il::int_t timeStep = 0; timeStep < SimParam.maxTimeStep; timeStep++) {

    //// SETUP TIMESTEP (but not declaration and constant variables)

    il::int_t totalNumberOfDofs = theMesh.numberOfPressDofs() + theMesh.numberOfDisplDofs();

    // Create source vector
    il::Array<double> f{totalNumberOfDofs, 0.};
    // Fill the source vector with injection... and (SIGMA0 - P)
    for (il::int_t i = 0; i < theSource.numberOfSources(); i++) {
      f[theSource.sourceNode[i]] = theSource.QatNode[i];
    }

    il::Array<double> iterativeSolution{totalNumberOfDofs, 0.};

    // Remember to add matrix of edge to collocation points
    // Remember to add segment characteristic (an array?)

    while (!theSolutionAtN1.isConverged) {

      // Update matrices
      il::Array2D<double> K = matrix_K(theMesh, theProperties, theSolutionAtN, theSource);
      // create it from basic_assembly( mesh, id, p, Ep). Also, using theSource or the historical variables to determine
      // which nodes/collocation points are active

      il::Array2D<double> Vw = matrix_Vw(theMesh, theProperties, theSolutionAtN, theSource);

      il::Array2D<double> N = matrix_N(theMesh, theProperties, theSolutionAtN, theSource);

      il::Array2D<double> Vp = matrix_Vp(theMesh, theProperties, theSolutionAtN, theSource);

      il::Array2D<double> L = matrix_L(theMesh, theProperties, theSolutionAtN, theSource);


      // Update force vectors

      // Solve system of equations (as the submatrices)

      // Update historical/depending variables

      // Compute errors and increase iterations

      // Swap variables

    }

  }

};


il::Array2D<double> matrix_lw(il::Array<double> widthB,
                              il::Array2D<int> col_row_i, il::Array2D<int> col_row_j,
                              const il::int_t &dof_dim, const int &p,
                              il::Array<double> element_size_all,
                              Parameters_fluid &fluid_parameters,
                              const il::Array2D<int> &id) {

  //widthB should be the ajusted-width profile at the previous moment
  // the size is (N+1)*(N+1)
  //int dof = 2 * (p + 1);
  il::int_t n = col_row_j(0, 0) - col_row_i(0, 0) + 1;//number of related elements
  if (col_row_i(0, 1) == 1) {
    n = n - 1; //left node is in the right half of the element
  };
  if (col_row_j(0, 1) == 0) {
    n = n - 1;
  };
  il::Array2D<double> LL{n + 1, n + 1, 0.};//contained node number
  il::Array<double> element_size{col_row_j(0, 0) - col_row_i(0, 0) + 1, 0.};//size of related element number
  il::Array<double> w_mid = average_open(widthB, col_row_i, col_row_j, id, dof_dim);//size of related element number
  for (il::int_t m = 0; m < col_row_j(0, 0) - col_row_i(0, 0) + 1; m++) {
    element_size[m] = element_size_all[m + col_row_i(0, 0)];
  };
  il::Array2D<double> density{col_row_j(0, 0) - col_row_i(0, 0) + 1, 2, 0.};
  take_submatrix(density,
                 col_row_i(0, 0),
                 col_row_j(0, 0),
                 0,
                 1,
                 fluid_parameters.density);//this should be size of related element number

  il::Array<double> rho_mid = average(density, il::io);
  il::Array<double> Kk = conductivities_newtonian_open(rho_mid, w_mid,
                                                       element_size,
                                                       fluid_parameters,
                                                       1., il::io);//size of related element number

  il::int_t dofj = 0;
  il::Array2D<int> Dofp = dofhandle_cg2d(2, col_row_j(0, 0) - col_row_i(0, 0) + 1, il::io);


  // Loop over the pressure nodes
  for (il::int_t s = 0; s < n + 1; ++s) {

    il::Array2D<int> ed = search(Dofp, int(s + col_row_i(0, 1)),
                                 il::io);//possible rows and coloumns

    for (il::int_t sc = 0; sc < ed.size(0); ++sc) { //number of related elements

      il::int_t ej = ed(sc, 0);//row number-element
      il::int_t ec = ed(sc, 1);//coloumn number-left or right node
      if (ec == 0) {//dofj corresponds to the other point of the element
        dofj = Dofp(ej, 1) - col_row_i(0, 1);
      } else {
        dofj = Dofp(ej, 0) - col_row_i(0, 1);
      }

      LL(s, s) = LL(s, s) - Kk[ej];
      if (dofj >= 0 && dofj < n + 1) {
        LL(s, dofj) = LL(s, dofj) + Kk[ej];
      }
    }
  }

  return LL;
}//should times the timestep to build the whole matrix

}