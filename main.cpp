//
// HFPx2D project.
//
// Created by Federico Ciardo on 09.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cmath>
#include <complex>
#include <iostream>
#include <sys/stat.h>

#include "il/Array.h"
#include "il/math.h"
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/Timer.h>
#include <il/Toml.h>
#include <il/linear_algebra.h>
#include <il/Timer.h>

//#include "il/linear_algebra/dense/factorization/LU.h"

#include "src/AssemblyDDM.h"
#include "src/DOF_Handles.h"
#include "src/Mesh.h"

////////////////////////////////////////////////////////////////////////////////
// analytical solution of the griffith-crack (ct pressure)
il::Array<double> griffithcrack(const il::Array<double>& x, double a, double Ep,
                                double sig) {
  double coef = 4. * sig / (Ep);
  il::Array<double> wsol{x.size(), 0.};

  for (int i = 0; i < x.size(); ++i) {
    if (std::abs(x[i]) < a) {
      wsol[i] = coef * sqrt(pow(a, 2) - pow(x[i], 2));
    }
  }
  return wsol;
}

////////////////////////////////////////////////////////////////////////////////
int main() {
  int nelts = 1000, p = 1 ;
  double h = 2. / (nelts);  //  element size

  il::Array<double> x{nelts+1};

  il::Array2D<double> xy{nelts+1, 2, 0.0};
  il::Array2D<int> myconn{nelts, 2, 0.0};
  il::Array2D<int> id{nelts, 4, 0};

  int ndof = (nelts) * (p + 1) * 2;  // number of dofs
  double Ep = 1.;                    // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
  for (int i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (int i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  // create mesh object
  hfp2d::Mesh mesh;
  mesh.set_values(xy, myconn);

  id= hfp2d::dofhandle_dp(2, nelts, p,il::io);  // dof handle for DDs

  // Read fluid parameters

  double CompressFluid;
  i = config.search("Fluid_compressibility");
  if (config.found(i) && config.value(i).is_floating_point()) {
    CompressFluid = config.value(i).to_floating_point();
  } else {
    CompressFluid = 0;
  }

  il::Array2D<double> K{ndof, ndof };

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
            << (nelts) * (p + 1) * 2 << "---" << ndof << "\n";
  std::cout << myconn.size(0) << "\n";

  // Read dilatancy parameters

  il::Timer timer{};
  timer.start();
  K=hfp2d::basic_assembly( mesh, id, p, Ep);  // passing p could be avoided here.
  timer.stop();
  std::cout << "------ " << timer.elapsed() << "  \n";
  std::cout << "---#---\n";
  result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));

  double d_wdilatancy;
  i = config.search("Slip_dw_for_dilatancy");
  if (config.found(i) && config.value(i).is_floating_point()) {
    d_wdilatancy = config.value(i).to_floating_point();
  } else {
    d_wdilatancy = 0;
  }

  // Read permeability parameters

  // use a direct solver
  il::Array<double> dd = il::linear_solve(K, f, il::io, status);  // lu_decomposition.solve(f);
//
//  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  // this piece of codes gets 1D mesh of x doubling the nodes of
  // adjacent elements (for comparison with analytical solution)
  for (int e = 0; e < nelts; ++e) {

    thex[i] = mesh.node(mesh.connectivity(e, 0), 0);
    thex[i + 1] = mesh.node(mesh.connectivity(e, 1), 0);
    i = i + 2;
  }

  wsol = griffithcrack(thex, 1., 1., 1.);  // call to analytical solution

  // Layer 1

//        std::cout << "x : " << thex[j] <<"..w anal:" << wsol[j] << " w num: "
//        << dd[j*2+1]<<  " rel error: " << rel_err << "\n";
  }

  std::cout << " end of code";
  status.ok();

  return 0;
}
