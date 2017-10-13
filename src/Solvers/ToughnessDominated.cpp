//
// Created by lorenzo on 10/13/17.
//

#include "ToughnessDominated.h"

namespace hfp2d{

double ToughnessDominated(int nelts) {

  int p = 1;
  double h = 2. / (nelts);  //  element size

  // il::Array<double> x{nelts + 1}; // Not needed

  il::Array2D<double> xy{nelts + 1, 2, 0.0};
  il::Array2D<il::int_t> myconn{nelts, 2, 0};
  il::Array2D<il::int_t> id_displ{nelts, 2*(p+1), 0};
  il::Array2D<il::int_t> id_press{nelts, 2, 0};
  il::Array<il::int_t> fracID {nelts,1};
  il::Array<il::int_t> matID {nelts,1};
  il::Array<il::int_t> condID {nelts,1};

  //int ndof = (nelts) * (p + 1) * 2;  // number of dofs
  double Ep = 1.;                    // Plane strain Young's modulus

  //  Array2D M(i, j) -> M(i + 1, j) (Ordre Fortran)
  //  Array2C M(i, j) -> M(i, j + 1) (Ordre C)

  // create a basic 1D mesh ....
  for (il::int_t i = 0; i < xy.size(0); ++i) {
    xy(i, 0) = -1. + i * h;
    xy(i, 1) = 0.;
  };

  for (il::int_t i = 0; i < myconn.size(0); ++i) {
    myconn(i, 0) = i;
    myconn(i, 1) = i + 1;
  };

  for (il::int_t i=0; i < nelts; i++){
    for (il::int_t j=0; j < 2*(p+1); j++){
      id_displ(i,j)=i * 2 * (p+1) + j;
    }
  }

  for (il::int_t i=0; i < nelts; i++){
    id_press(i,0)=i;
    id_press(i,1)=i+1;
  }

  // Mesh initialization
  hfp2d::Mesh mesh(p,xy,myconn,id_displ,id_press,fracID,matID,condID);
  il::int_t ndof = mesh.numDisplDofs();

  // Elastic properties initialization
  hfp2d::ElasticProperties myelas(1.0e7, 0.);
  std::cout << "EP :" << myelas.Ep() << "\n";


  //il::Array2D<int> id = hfp2d::dofhandle_dp(2, nelts, p, il::io);  // dof handle for DDs

  // some definitions needed for matrix assembly
  //il::Array2D<double> xe{2, 2, 0}, xec{2, 2, 0};

  //  SegmentData mysege,mysegc;

  std::cout << "Number of elements : " << mesh.nelts() << "\n";
  //std::cout << "Number of dofs :" << id.size(0) * id.size(1) << "---"
  //          << (nelts) * (p + 1) * 2 << "---" << mesh.numDisplDofs() << "\n";
  //std::cout << myconn.size(0) << "\n";

  std::cout << "------ Assembly: \t";
//  std::time_t result = std::time(nullptr);
//  std::cout << std::asctime(std::localtime(&result));
  il::Timer timer{};
  timer.start();

  il::Array2D<double> K{ndof, ndof};

  K = hfp2d::basic_assembly_new(mesh, myelas,
                                hfp2d::normal_shear_stress_kernel_dp1_dd,
                                0.);  // passing p could be avoided here

  il::Array<double> F(ndof,0.0);

  timer.stop();
  std::cout << timer.elapsed() << " s" << "  \n";

  std::cout << "------ Extract: \t";
  timer.reset();
  timer.start();

  // Set the initial active set
  il::Array<bool> activeSet(mesh.numElems(),false);

  //// SET INITIAL SOURCE POSITION
  // Set the initial location of source (middle point)
  // todo: use a search function to determine which elements are active due to initial source
  il::int_t sourceElem = nelts / 2;
  if(nelts%2==0) { // even number of elements
    activeSet[sourceElem-1]=true;
    activeSet[sourceElem]=true;
  } else {
    activeSet[sourceElem]=true;
  }

  // Count how many active elements are there
  il::int_t numActive=0;
  for(il::int_t i=0; i<mesh.numElems(); i++){
    if(activeSet[i]) { numActive++ ;}
  }

  // Save the active elements
  il::int_t counter = 0;
  il::Array<il::int_t> listActive (numActive);
  for(il::int_t i=0; i<mesh.numElems(); i++){
    if(activeSet[i]){
      listActive[counter] = i;
      counter++;
    }
  }

  // Create K matrix and force vector of the size of the active set + 1 (pressure)
  il::int_t blockSize = mesh.numDisplDofsPerElem();
  il::int_t actDispl = numActive*blockSize;
  il::Array2D<double> Kact(actDispl+1, actDispl+1);
  il::Array<double> Fact(actDispl+1);

  // Fill the matrix with the values of active elements
  counter = 0; // counter of active elements != total number of elements

  for(il::int_t i=0; i<numActive; i++){
    for(il::int_t j=0; j<numActive; j++){

      // the block to copy starts from xStart,yStart
      il::int_t xStart=i*blockSize;
      il::int_t yStart=j*blockSize;

      // save the dof values of the corresponding active element
      // in the blocks of the counter
      for(il::int_t k=0; k<blockSize; k++){
        for(il::int_t l=0; l<blockSize; l++){

          il::int_t dof1 = mesh.dofDispl(listActive[i],k);
          il::int_t dof2 = mesh.dofDispl(listActive[j],l);

          double extValue = K(dof1,dof2);
          Kact(xStart+k,yStart+l) = extValue;

        }
        //save F in Fact
        Fact[xStart+k]=F[mesh.dofDispl(listActive[i],k)];

      }
    }
  }


  // The additional line is full of 1, the additional column is full of -1
  // The additional diagonal term is zero
  for(il::int_t i=0; i < actDispl; i++){
    Kact(i, actDispl)=-1.0;
  }
  for(il::int_t i=0; i < actDispl; i++) {
    Kact(actDispl, i) = 1.0;
  }
  Kact(actDispl,actDispl)=0.0;
  Fact[actDispl]=0.1;

  timer.stop();
  std::cout << timer.elapsed() << " s" << "  \n";

  il::Status status;
//timing reset
  std::cout << "------ Solution:\t";
  timer.reset();
  timer.start();

//  il::LU<il::Array2D<double> > lu_decomposition(K, il::io, status);
//  il::Array<double> Xsol= lu_decomposition.solve(std::move(F));

  il::Array<double> Xsol = il::linearSolve(Kact, Fact, il::io, status);

  status.ok();


  timer.stop();
  std::cout << timer.elapsed() << " s \n";
  std::cout << "---#---\n";

  std::cout << "Opening values: \n" << std::endl;
  for(il::int_t i=0; i<actDispl; i++){
    std::cout << i << "\t" << Xsol[i] << std::endl;
  }

  //// Reconstruct big vector

  il::Array<double> Xfull(mesh.numDisplDofs(),0.0);
  for(il::int_t i=0; i<numActive; i++){
    for(il::int_t j=0; j<blockSize; j++) {
      il::int_t numElem = listActive[i];
      Xfull[mesh.dofDispl(numElem, j)] = Xsol[i * blockSize + j];
    }
  }

  std::cout << "Opening values full: \n" << std::endl;
  for(il::int_t i=0; i<Xfull.size(); i++){
    std::cout << i << "\t" << Xfull[i] << std::endl;
  }

  //// MULTIPLE TIME STEPS
//  for(il::int_t timeStep=0; timeStep<10; timeStep++){
//
//    Fact[actDispl] = timeStep*0.001;
//
//    il::Array<double> Xsol = il::linearSolve(Kact, Fact, il::io, status);
//    status.ok();
//
//    std::cout << "Opening values at timestep " << timeStep << " :" << std::endl;
//    for(il::int_t i=0; i<actDispl; i++){
//      std::cout << i << "\t" << Xsol[i] << std::endl;
//    }
//
//    il::Array<double> stressAtColl = il::dot(Kact,Xsol);
//    std::cout << "Stress values at timestep " << timeStep << " :" << std::endl;
//    for(il::int_t i=0; i<actDispl; i++){
//      std::cout << i << "\t" << stressAtColl[i] << std::endl;
//    }
//
//  }

//  // solve a constant pressurized crack problem...
//  il::Array<double> f{ndof, -1.};
//  // just opening dds - set shear loads to zero
//  for (il::int_t i = 0; i < ndof / 2; ++i) {
//    f[2 * i] = 0;
//  }
//
//  v
//  // use a direct solver
//  il::Array<double> dd = il::linearSolve(K, f, il::io, status);
//  //
//  status.ok();

  std::cout << "End of the analysis" << std::endl;



  return Xsol[actDispl];

}



}
