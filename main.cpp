//
// HFPx2D project.
//
// Created by Brice Lecampion on 06.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <cmath>
#include <complex>
#include <iostream>
#include <string>

#include <il/Array.h>
#include <il/math.h>
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/linear_algebra/dense/factorization/LU.h>
#include "AssemblyDDM.h"
#include "DOF_Handles.h"
#include "Mesh.h"
#include "Stress.h"
#include "Coh_Propagation.h"
#include "Coh_Prop_new.h"
#include <fstream>



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
  int n = 1001, p = 1;
  double h = 2. / (n - 1);  //  element size

  il::Array<double> x{n};

  il::Array2D<double> xy{n, 2, 0.0};
  il::Array2D<int> myconn{n - 1, 2, 0.0};
  il::Array2D<int> id{n - 1, 4, 0};

  int ndof = (n - 1) * (p + 1) * 2;  // number of dofs
  double Ep = 100.;                    // Plane strain Young's modulus

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

  dofhandle_DG2D(id, 2, n - 1, p);  // dof handle for DDs

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

  BasicAssembly(K, mesh, id, p, Ep);  // passing p could be avoided here.

  std::cout << "------\n";
  result = std::time(nullptr);
  std::cout << std::asctime(std::localtime(&result));
  // solve a constant pressurized crack problem...
  il::Array<double> f{ndof, -1.};
  // just opening dds - set shear loads to zero
  for (int i = 0; i < ndof / 2; ++i) {
    f[2 * i] = 0;
  }
//  il::Status status1;
//  il::Array2D<double> k1{4,4,1};
//  il::Array<double> f1{4,1};
//  il::Array<double> ddd=linear_solve(k1,f1,il::io,status1);
//  std::cout<<ddd[1];



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
      linear_solve(K, f, il::io, status);  // lu_decomposition.solve(f);//

//    //we add here to test the stress calculation code part;
//    il::Array<double> stress;
//    stress=stressoutput(mesh,id,p,Ep,2.0,0.0,dd);
//    std::cout << stress[1];

  //we add here to test the propagation code.
  il::Array2D<double> widthlist;
  il::Array<double> plist;
    il::Array<double> l_coh;
    il::Array<double> l_c;
  Material material;
  Initial_condition initial_condition;

  material_condition(material, initial_condition,
                    0.0001, 2.,100.0,0.,
          0.00001,0.001,0.01,0.);
//  void material_condition(Material &material, Initial_condition &initial_condition,
//                          double wc1, double sigma_t1,double Ep1,double U1,
//                          double pini1,double Q01,double timestep1,double sigma01);



    il::Array<double> widthB;

    il::Status status2;
  propagation_loop_new(widthlist,plist,l_coh,l_c,
                   mesh,id,p,material,initial_condition,490,510,25,status2);
  il::Array<double> xlist{2*mesh.nelts(),0.};
  get_xlist(xlist,mesh);

  std::ofstream fout;
  fout.open("output2.txt");
//    fout<<initial_condition.timestep<<"\n";
//    for(int pp=0;pp<plist.size();++pp){
//        fout<<plist[pp]<<"\t";
//    }
//    fout<<"\n";
//    for (int xx = 0; xx <xlist.size() ; ++xx) {
//        fout<<xlist[xx]<<"\t";
//    }
//    fout<<"\n";
  for(int qq=0;qq<widthlist.size(0);++qq){
    for(int qqq=0;qqq<widthlist.size(1);++qqq){
      fout<<widthlist(qq,qqq)<<"\t";
    }
    fout<<"\n";
  }
  fout.close();

  std::ofstream fout1;
  fout1.open("outputpressure.txt");
//    fout<<initial_condition.timestep<<"\n";
//    for(int pp=0;pp<plist.size();++pp){
//        fout<<plist[pp]<<"\t";
//    }
//    fout<<"\n";
//    for (int xx = 0; xx <xlist.size() ; ++xx) {
//        fout<<xlist[xx]<<"\t";
//    }
//    fout<<"\n";
  for(int mm=0;mm<plist.size();++mm){
      fout1<<plist[mm]<<"\n";
    }

  fout1.close();


  // Analytical solution at nodes
  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

  int i = 0;
  for (int e = 0; e < n - 1;
       ++e) {  // this piece of codes gets 1D mesh of x doubling the nodes of
               // adjacent elements (for comparison with analytical solution)
    thex[i] = mesh.Coor(mesh.conn(e, 0), 0);
    thex[i + 1] = mesh.Coor(mesh.conn(e, 1), 0);
    i = i + 2;
  }

  wsol = griffithcrack(thex, 1., 1., 1.);  // call to analytical solution

  // printing out the comparisons for each nodes (left and right values at each
  // nodes due to the piece-wise constant nature of the solution)...
  double rel_err;
  for (int j = 0; j < ndof / 2; ++j) {
    rel_err = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];

    //    std::cout << "x : " << thex[j] <<"..w anal:" << wsol[j] << " w num: "
    //    << dd[j*2+1]<<  " rel error: " << rel_err << "\n";
  }

  return 0;
}


