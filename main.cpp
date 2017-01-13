//
// HFPx2D project.
//
// Created by Federico Ciardo on 09.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>

#include <il/Array.h>
#include <il/math.h>
//#include <il/Array2C.h>
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "FromEdgeToCol.h"
#include "FVM.h"
#include "Friction.h"
#include "Mesh.h"
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>

/*////////////////////////////////////////////////////////////////////////////////
// analytical solution of the griffith-crack (ct pressure)
il::Array<double> griffithcrack(const il::Array<double> &x, double a, double Ep,
                                double sig) {
  double coef = 4. * sig / (Ep);
  il::Array<double> wsol{x.size(), 0.};

  for (int i = 0; i < x.size(); ++i) {
    if (std::abs(x[i]) < a) {
      wsol[i] = coef * sqrt(pow(a, 2) - pow(x[i], 2));
    }
  }
  return wsol;
}*/

/*////////////////////////////////////////////////////////////////////////////////
int main() {
  int n = 200, p = 1;
  double h = 2. / (n - 1); //  element size

  il::Array<double> x{n};

  il::Array2D<double> xy{n, 2, 0.0};
  il::Array2D<int> myconn{n - 1, 2, 0.0};
  il::Array2D<int> id{n - 1, 4, 0};

  int ndof = (n - 1) * (p + 1) * 2; // number of dofs
  double Ep = 1.;                   // Plane strain Young's modulus

  //  std::complex(double re = 0.0, double im = 0.0) myC2;
  //  myC.real(2.);
  //  myC.imag(1.);
  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  }

  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  }

  // create mesh object
  Mesh mesh;
  mesh.set_values(xy, myconn);

  dofhandle_DG2D(id, 2, n - 1, p); // dof handle for DDs

  // some definitions needed for matrix assembly
  il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  //  SegmentCharacteristic mysege,mysegc;

  il::Array2D<double> K{ndof, ndof, 0.};

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
            << (n - 1) * (p + 1) * 2 << "---" << ndof << "\n";
  std::cout << myconn.size(0) << "\n";

  std::cout << "------\n";
  std::time_t result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  BasicAssembly(K, mesh, id, p, Ep); // passing p could be avoided here.

  std::cout << "------\n";
  result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, -1.};
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }

  il::Status status;
  // example if LU decomposition
  //  il::LU<il::Array2D<double>> lu_decomposition(K, il::io, status);
  //  if (!status.ok()) {
  //    // The matrix is singular to the machine precision. You should deal with
  //    the error.
  //  }
  // il::Array<double> dd = lu_decomposition.solve(f);

  // use a direct solver
  il::Array<double> dd =
      linear_solve(K, f, il::io, status); // lu_decomposition.solve(f);

  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  for (int e = 0; e < n - 1;
       ++e) { // this piece of codes gets 1D mesh of x doubling the nodes of
              // adjacent elements (for comparison with analytical solution)
    thex[i] = mesh.Coor(mesh.conn(e, 0), 0);
    thex[i + 1] = mesh.Coor(mesh.conn(e, 1), 0);
    i = i + 2;
  }

  wsol = griffithcrack(thex, 1., 1., 1.); // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  double rel_err;
  for (int j = 0; j < ndof / 2; ++j) {
    rel_err = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    //    std::cout << "x : " << thex[j] <<"..w anal:" << wsol[j] << " w num: "
    //    << dd[j*2+1]<<  " rel error: " << rel_err << "\n";
  }

  return 0;
}*/

il::Array2D<double>
AnalyticalSolution(double Dp, Mesh mesh,
                   il::Array2D<int> myconn);          // Function prototype
il::int_t find(il::Array<double> arr, double_t seek); // Function prototype
il::Array<int> find2D(il::Array2D<double> &arr2D,
                      double_t seek);       // Function prototype
double_t max2D(il::Array2D<double> &arr2D); // Function prototype
double_t min2D(il::Array2D<double> &arr2D); // Function prototype



////////////////////////////////////////////////////////////////////////////////
int main() {
  // Declare the input variables
  double Ep, Density, Cohes, CompressFluid, Visc, Init_dil, Incr_dil, d_wf,
      d_wd, Peak_fric, Resid_fric, sigma_n0, P0, Source;
  int Nnodes;

  std::ifstream read_file(
      "/Users/federicociardo/ClionProjects/HFPx2D/InputData.txt"); // Read the
                                                                   // file that
                                                                   // contains
                                                                   // the input
                                                                   // data
  IL_ASSERT(
      read_file
          .is_open()); // Check whether the input data file is opened or not
  read_file >> Nnodes >> Ep >> Density >> Cohes >> CompressFluid >> Visc >>
      Init_dil >> Incr_dil >> d_wf >> d_wd >> Peak_fric >> Resid_fric >>
      sigma_n0 >> P0 >> Source; // Read the data
  read_file.close();            // Close the data file

  // Create the vector of density, initial stress state (uniform for now),
  // initial/ambient pore pressure

  il::Array2D<double> rho{Nnodes - 1, 2, 0.}; // fluid density vector
  // {{rho1_left,rho1_right},{rho2_left,rho2_right}..}   Remember:
  // continuos linear varioation -> rho1_right = rho2_left and so on ..

  for (il::int_t i{0}; i < rho.size(0); ++i) {
    for (il::int_t j{0}; j < rho.size(1); ++j) {

      rho(i, j) = Density;
    }
  }

  double sigma_s0 = 0.55 * sigma_n0,
         Dp = 0.5 * sigma_n0; // Declare some variables -> sigma_s0 is the
                              // ambient shear stress (which is related to the
                              // ambient normal stress)
  int Nelts = Nnodes - 1;

  int NcollPoints = 2 * Nelts; // Number of collocation points

  il::Array2D<double> Sigma0{NcollPoints, 2, 0.}; // matrix -> {Ambient shear,
                                                  // Ambient normal stress} for
                                                  // each collocation points

  for (il::int_t i{0}; i < Sigma0.size(0); ++i) {
    Sigma0(i, 0) = sigma_s0;
    Sigma0(i, 1) = sigma_n0;
  }

  il::Array2D<double> Amb_press{2 * Nnodes - 2, 2,
                                0.}; // Ambient pore pressure at nodal points
  // {{P01_left,P01_right},{P02_left,P02_right}..}
  // Remember -> Pressure varies linearly and
  // continuously over the elements -> P01_right =
  // P02_left

  for (il::int_t k{0}; k < Amb_press.size(0); ++k) {
    for (il::int_t i{0}; i < Amb_press.size(1); ++i) {
      Amb_press(k, i) = P0;
    }
  }

  int p = 1;                  // interpolation order
  double l = 1.;              // half crack/fault length
  double h = (2 * l) / Nelts; // element size

  // Declare variables for mesh
  //  il::Array<double> x{Nnodes}; // vector that contains x-coord of 1D mesh
  il::Array2D<double> xy{
      Nnodes, 2, 0.0}; // matrix that contains x- and y- coordinates of 1D mesh
  il::Array2D<int> myconn{Nelts, 2, 0.0}; // connectivity matrix
  il::Array2D<int> id{Nelts, 4, 0};       // slip/opening dof

      int Ndof = (Nnodes - 1) * (p + 1) * 2; // number of dofs

  // create a basic 1D mesh ....
  for (il::int_t i{0}; i < xy.size(0); ++i) {
    xy(i, 0) = -l + i * h;
    xy(i, 1) = 0.;
  }

  for (il::int_t i{0}; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  }

  // create mesh object
  Mesh mesh;
  mesh.set_values(xy, myconn);

  il::Array2D<double> Pinit{(mesh.Coor).size(0), 2, 0.};
  Pinit = AnalyticalSolution(Dp, mesh, myconn); // Initial pore pressure profile

  double_t MaxPinit;
  MaxPinit = max2D(Pinit); // Find the max value in Pinit

  il::Array<int> idx;
  idx = find2D(Pinit, MaxPinit); // Find the position of the max value in Pinit

  il::Array<double> S{Nnodes, 0.};
  S[idx[0]] = Source; // Source vector

    dofhandle_DG2D(id, 2, Nelts, p); // dof handle for DDs

 /*   il::Array<int> t{2,0.};
    t[0] = 3;
    t[1] = 4;
    il::int_t de;

    for (il::int_t m = 0; m < t.size(); ++m) {

        if(t[m] != 3) de = t[m];

    }*/



//    test = ConductivitiesNewtonian(Visc,mesh,rho,Pinit,Incr_dil,d_wd,Init_dil);
//    test = Quarter(xy);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Analytical solution for pore pressure distribution in a fracture/fault with
// constant permeability subjected to a constant overpressure Dp
// It return a matrix {{Pinit1_left,Pinit1_right},{Pinit2_left,Pinit2_right}..}
// Remember -> pressure varies linearly and continuously over the elements ->
// Pinit1_right = Pinit2_left

il::Array2D<double> AnalyticalSolution(double Dp, Mesh mesh,
                                       il::Array2D<int> myconn) {

  il::Array2D<double> Pinit{((mesh.Coor).size(0)) - 1, 2, 0.};

  for (il::int_t i{0}; i < Pinit.size(0); ++i) {

    for (il::int_t j{0}; j < Pinit.size(1); ++j) {

      Pinit(i, j) =
          Dp * (erfc(fabs(mesh.Coor(myconn(i, j), 0) / sqrt(10 * 0.5)) / 2));
    }
  }

  return Pinit;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the index of a given value ("seek") in an array.
// arr -> array in which we want to find the index of  given value
// seek ->  value for which we want to find out the index

il::int_t find(il::Array<double> arr, double_t seek) {

  for (il::int_t i = 0; i < arr.size(); ++i) {

    if (arr[i] == seek)
      return i;
  }

  return -1;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the index ( IndRow,IndColumn ) of a given value ("seek")
// in an array2D.
// arr2D -> array2D in which we want to find the index of  given value
// seek -> value for which we want to find out the index
// It return an array that contain the {N.row, N.col} of the seek value

il::Array<int> find2D(il::Array2D<double> &arr2D, double_t seek) {

  il::Array<int> outp{2};

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (arr2D(i, j) == seek)
        outp[0] = i, outp[1] = j;
    }
  }

  return outp;
}


////////////////////////////////////////////////////////////////////////////////
// Function that return the max value in an array2D.
// arr2D -> array2D in which we want to find the max
// seek -> seek that we look for

double_t max2D(il::Array2D<double> &arr2D) {

  double_t max;
  max = arr2D(0., 0.);

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (max < arr2D(i, j)) {

        max = arr2D(i, j);
      }
    }
  }

  return max;
}

////////////////////////////////////////////////////////////////////////////////
// Function that return the min value in an array2D.
// arr2D -> array2D in which we want to find the min
// seek -> seek that we look for

double_t min2D(il::Array2D<double> &arr2D) {

  double_t min;
  min = arr2D(0., 0.);

  for (il::int_t i = 0; i < arr2D.size(0); ++i) {

    for (il::int_t j = 0; j < arr2D.size(1); ++j) {

      if (min > arr2D(i, j)) {

        min = arr2D(i, j);
      }
    }
  }

  return min;
}


