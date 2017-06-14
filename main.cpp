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
#include <fstream>

#include "il/Array.h"
#include "il/math.h"
//#include <il/Array2C.h>
#include <il/StaticArray.h>
#include <il/linear_algebra.h>
#include <il/Timer.h>

//#include "il/linear_algebra/dense/factorization/LU.h"

#include "src/AssemblyDDM.h"
#include "src/DOF_Handles.h"
#include "src/Mesh.h"
#include "Coh_Prop_Col.h"
#include "src/FVM.h"
#include "Coh_Col_Partial.h"
#include "Coh_Linear_softening.h"
//#include "Viscosity.h"
#include "Viscosity_all_nodeselemt.h"
//#include "Viscosity_right_side.h"

// You should write multiple files!
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
        int nelts = 100, p = 1 ;
        double h = 2. / (nelts);  //  element size before h=2./(nelts);

        il::Array<double> x{nelts+1};

        il::Array2D<double> xy{nelts+1, 2, 0.0};
        il::Array2D<int> myconn{nelts, 2, 0};
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

        id=hfp2d::dofhandle_dg_full2d( 2, nelts, p,il::io);  // dof handle for DDs

        // some definitions needed for matrix assembly
        il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

        //  SegmentCharacteristic mysege,mysegc;

        il::Array2D<double> K{ndof, ndof };

        std::cout << "Number of elements : " << mesh.nelts() << "\n";
        std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
                  << (nelts) * (p + 1) * 2 << "---" << ndof << "\n";
        std::cout << myconn.size(0) << "\n";

        std::cout << "------\n";
        std::time_t result = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&result));

        il::Timer timer{};
        timer.start();
        K=hfp2d::basic_assembly( mesh, id, p, Ep);  // passing p could be avoided here.
        timer.stop();
        std::cout << "------ " << timer.elapsed() << "  \n";
        std::cout << "---#---\n";
        result = std::time(nullptr);
        std::cout << std::asctime(std::localtime(&result));

        // solve a constant pressurized crack problem...
        il::Array<double> f{ndof, -1.};
        // just opening dds - set shear loads to zero
        for (int i = 0; i < ndof / 2; ++i) {
            f[2 * i] = 0;
        }

        il::Status status;

        il::Array<double> dd = il::linear_solve(K, f, il::io, status);  // lu_decomposition.solve(f);

    std::cout << "------\n";
    std::time_t resultwhole = std::time(nullptr);
    std::cout << std::asctime(std::localtime(&resultwhole));
    il::Timer timerwhole{};
    timerwhole.start();




  //we add here to test the propagation code.
  il::Array2C<double> widthlist;
  il::Array<double> plist;
    il::Array2D<double> plist_2d;
    il::Array<double> l_coh;
    il::Array<double> l_c;
    il::Array<double> energy;
    il::Array<double> energy_g;
    il::Array<double> energy_f;
    il::Array<double> energy_p;
    il::Array<double> energy_coh;
    il::Array<double> energy_j_integral;
    il::Array2C<double> cohlist;
    il::Array2C<double> stresslist;
    il::Array<double> volume_change;
    il::Array2D<double> volume_vary_list;
    il::Array2D<double> elastic_vary_list;
  hfp2d::Material material;
  hfp2d::Initial_condition initial_condition;

    //For Dugdale cohesive law

//    hfp2d::material_condition_col (0.001, 2.,100.0,0.,
//                       0.00001,0.0005,1.,0.,il::io,material, initial_condition);

//  void material_condition(Material &material, Initial_condition &initial_condition,
//                          double wc1, double sigma_t1,double Ep1,double U1,
//                          double pini1,double Q01,double timestep1,double sigma01);
    //wc before=0.0001 changes on the 7th April
//tstepbefore=0.01 changes on the 7th April


    //For linear softening cohesive law
    // sigma_T should be doubled in order to have the same Gc
    // or wc should be doubled

    hfp2d::material_condition_col (0.001*2, 2.,100.0,0.,
                                   0.00001,0.0005,2.,0.,il::io,material, initial_condition);

    //For exponential cohesive law
    //sigma T or wc should be times a constant 6/exp(1.0)/0.95
    //hfp2d::material_condition_col (0.001, 2.*6/exp(1.0)/0.95,100.0,0.,
     //                              0.00001,0.0040,0.1,0.,il::io,material, initial_condition);







    il::Array<double> widthB;
    il::Array<int> mvalue;

    int nstep=1000;
    int break_time=0;

    il::Status status2;

    //fully filled calculation

//    hfp2d::propagation_loop_col(mesh,id,p,material,initial_condition,799,801,nstep,status2,il::io,widthlist,plist,l_coh,l_c,cohlist,mvalue,break_time,stresslist,energy);
//
//    hfp2d::energy_output(widthlist,plist,l_c,l_coh,material,mesh,id,p,2,initial_condition,il::io,energy_f,energy_coh,energy_j_integral);
//
//    volume_change=hfp2d::volume_output(widthlist,2);
//
//    stresslist=hfp2d::deal_with_stress(stresslist,cohlist,plist,initial_condition);


    //partially filled calculation

//    hfp2d::propagation_loop_col_partial(mesh,id,p,material,initial_condition,399,401,nstep,status2,il::io,widthlist,plist,l_coh,l_c,cohlist,mvalue,break_time,stresslist,energy);
//
//    volume_change=hfp2d::volume_output(widthlist,2);
//
//    hfp2d::energy_output_partial(widthlist,plist,l_c,l_coh,cohlist,material,mesh,id,p,2,initial_condition,il::io,energy_f,energy_coh,energy_j_integral);


    //linear softening cohesive law, based on the fully-filled case

//    hfp2d::propagation_loop_linear(mesh,id,p,material,initial_condition,799,801,nstep,status2,il::io,widthlist,plist,l_coh,l_c,cohlist,mvalue,break_time,stresslist,energy);
//
//    hfp2d::energy_output(widthlist,plist,l_c,l_coh,material,mesh,id,p,2,initial_condition,il::io,energy_f,energy_coh,energy_j_integral);
//
//    volume_change=hfp2d::volume_output(widthlist,2);
//
//    stresslist=hfp2d::deal_with_stress(stresslist,cohlist,plist,initial_condition);


    // Considering the viscosity

    // Create the fluid density matrix (rho)
    // Fluid density matrix {{rho_1left, rho_1right},{rho_2left,rho_2right} ..}
    // Remember: continuous linear variation -> rho_2right = rho_1left and so on..

    double CompressFluid=0;
    double Visc=0.001;
    double Density=1.;

    il::Array2D<double> rho{nelts, 2, Density};
    il::Array2D<double> error_matrix;


    // Set the structure members of fluid
    hfp2d::Parameters_fluid fluid_parameters;
    fluid_parameters.compressibility = CompressFluid;
    fluid_parameters.density = rho;
    fluid_parameters.viscosity = Visc;
    hfp2d::propagation_loop_visco(mesh,id,p,material,initial_condition,49,50,nstep,status2,fluid_parameters,il::io,widthlist,plist_2d,l_coh,l_c,cohlist,mvalue,break_time,stresslist,energy,volume_vary_list,elastic_vary_list,error_matrix);

    timerwhole.stop();
    std::cout << "------ " << timerwhole.elapsed() << "  \n";
    std::cout << "---#---\n";
    resultwhole = std::time(nullptr);
    std::cout << std::asctime(std::localtime(&resultwhole));






    if(break_time!=nstep){
        std::cout<<"Oups! At the "<<break_time<< "th time step, the fracture reaches the mesh end point"<<"\n";
    }
    std::cout<<"To draw the curves,nstep="<<break_time<<"\n";




    std::ofstream foutlc;
    foutlc.open("cracklength.txt");
    for(int lca=0;lca<break_time;++lca){
        foutlc<<l_c[lca]<<"\n";
    }
    foutlc.close();


    std::ofstream foutit;
    foutit.open("iteration.txt");
    for(int itera=0;itera<break_time;++itera){
        foutit<<mvalue[itera]<<"\n";
    }
    foutit.close();

    il::Array<double> xlist{2*mesh.nelts(),0.};
    hfp2d::get_xlist_col(xlist,mesh);
  std::ofstream fout;
  fout.open("outputcn1.txt");


  for(int qq=0;qq<break_time;++qq){//qq<widthlist.size(0)
    for(int qqq=0;qqq<widthlist.size(1);++qqq){
      fout<<widthlist(qq,qqq)<<"\t";
    }
    fout<<"\n";
  }
  fout.close();
    

//  std::ofstream fout1;
//  fout1.open("outputpressurecn1.txt");
//
//  for(int mm=0;mm<break_time+1;++mm){//mm<plist.size()
//      fout1<<plist[mm]<<"\n";
//    }
//
//  fout1.close();


//plot only for viscosity case

    std::ofstream fplist;
    fplist.open("outputplist.txt");

    for(int qp=0;qp<break_time+1;++qp){//qq<widthlist.size(0)
        for(int qqp=0;qqp<plist_2d.size(1);++qqp){
            fplist<<plist_2d(qp,qqp)<<"\t";
        }
        fplist<<"\n";
    }
    fplist.close();


    std::ofstream fout2;
    fout2.open("outputlcohcn1.txt");
    for(int cc=0;cc<break_time;++cc){//cc<l_coh.size()
        fout2<<l_coh[cc]<<"\n";
    }
    fout2.close();



//out put for viscosity case
    std::ofstream foutvolume_list;
    foutvolume_list.open("outputvlist.txt");

    for(int vq=0;vq<break_time;++vq){//qq<widthlist.size(0)
        for(int vqq=0;vqq<volume_vary_list.size(1);++vqq){
            foutvolume_list<<volume_vary_list(vq,vqq)<<"\t";
        }
        foutvolume_list<<"\n";
    }
    foutvolume_list.close();


    std::ofstream felas;
    felas.open("outputelas.txt");


    for(int eq=0;eq<break_time;++eq){//qq<widthlist.size(0)
        for(int eqq=0;eqq<elastic_vary_list.size(1);++eqq){
            felas<<elastic_vary_list(eq,eqq)<<"\t";
        }
        felas<<"\n";
    }
    felas.close();

//output only for viso case, check the iteration scheme
    std::ofstream fviscoitera;
    fviscoitera.open("outputviscoitera.txt");
    for(int erv=0;erv<break_time;++erv){//qq<widthlist.size(0)
        for(int erc=0;erc<error_matrix.size(1);++erc){
            fviscoitera<<error_matrix(erv,erc)<<"\t";
        }
        fviscoitera<<"\n";
    }
    fviscoitera.close();





//energy output

//    std::ofstream fenergy;
//    fenergy.open("outputenergyg.txt");
//    for(int ee=0;ee<break_time;++ee){
//        fenergy<<energy_g[ee]<<"\n";
//    }
//    fenergy.close();
//
//    std::ofstream fenergy_p;
//    fenergy_p.open("outputenergyp.txt");
//    for(int ep=0;ep<break_time;++ep){
//        fenergy_p<<energy_p[ep]<<"\n";
//    }
//    fenergy_p.close();


//energy output maybe more important

//    std::ofstream fenergy_f;
//    fenergy_f.open("outputenergyf.txt");
//    for(int ef=0;ef<break_time;++ef){
//        fenergy_f<<energy_f[ef]<<"\n";
//    }
//    fenergy_f.close();
//
//    std::ofstream fenergy_coh;
//    fenergy_coh.open("outputenergycoh.txt");
//    for(int ecoh=0;ecoh<break_time;++ecoh){
//        fenergy_coh<<energy_coh[ecoh]<<"\n";
//    }
//    fenergy_coh.close();
//
//    std::ofstream fenergy_j;
//    fenergy_j.open("outputenergyj.txt");
//    for(int ej=0;ej<break_time+1;++ej){
//        fenergy_j<<energy_j_integral[ej]<<"\n";
//    }
//    fenergy_j.close();
//
//    std::ofstream foutf;
//    foutf.open("outputcohfcn1.txt");
//    for(int cf=0;cf<break_time;++cf){//cf<cohlist.size(0)
//        for(int cff=0;cff<cohlist.size(1);++cff){
//            foutf<<cohlist(cf,cff)<<"\t";
//        }
//        foutf<<"\n";
//    }
//    foutf.close();
//
    std::ofstream foutstress;
    foutstress.open("outputstresscn1.txt");
    for(int cstr=0;cstr<break_time;++cstr){//cf<cohlist.size(0)
        for(int cstre=0;cstre<cohlist.size(1);++cstre){
            foutstress<<stresslist(cstr,cstre)<<"\t";
        }
        foutstress<<"\n";
    }
    foutstress.close();
//
//    std::ofstream fvol;
//    fvol.open("outputvol.txt");
//
//    for(int v=0;v<break_time;++v){//mm<plist.size()
//        fvol<<volume_change[v]<<"\n";
//    }
//
//    fvol.close();






  // Analytical solution at nodes
//  il::Array<double> thex{ndof / 2, 0}, wsol{ndof / 2, 0};

//  int i = 0;
//  for (int e = 0; e < n - 1;
//       ++e) {  // this piece of codes gets 1D mesh of x doubling the nodes of
//               // adjacent elements (for comparison with analytical solution)
//    thex[i] = mesh.Coor(mesh.conn(e, 0), 0);
//    thex[i + 1] = mesh.Coor(mesh.conn(e, 1), 0);
//    i = i + 2;
//  }
//
//  wsol = griffithcrack(thex, 1., 1., 1.);  // call to analytical solution
//
//  // printing out the comparisons for each nodes (left and right values at each
//  // nodes due to the piece-wise constant nature of the solution)...
//  double rel_err;
//  for (int j = 0; j < ndof / 2; ++j) {
//    rel_err = sqrt(pow(dd[j * 2 + 1] - wsol[j], 2)) / wsol[j];
//
//    //    std::cout << "x : " << thex[j] <<"..w anal:" << wsol[j] << " w num: "
//    //    << dd[j*2+1]<<  " rel error: " << rel_err << "\n";
//  }
//
//  return 0;
}


