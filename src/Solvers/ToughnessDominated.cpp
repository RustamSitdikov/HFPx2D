//
// Created by lorenzo on 10/13/17.
//

#include <iomanip>
#include "ToughnessDominated.h"

namespace hfp2d
{

// Convention of signs
//      A dw - dp = - A w + p - sigma_0 - sigma_coh( w + dw )

double
ToughnessDominated(int nelts)
{

    int p = 1;
    double h = 2. / (nelts); //  element size
    il::Array2D<double> xy{nelts + 1, 2, 0.0};
    il::Array2D<il::int_t> myconn{nelts, 2, 0};
    il::Array2D<il::int_t> id_displ{nelts, 2 * (p + 1), 0};
    il::Array2D<il::int_t> id_press{nelts, 2, 0};
    il::Array<il::int_t> fracID{nelts, 1};
    il::Array<il::int_t> matID{nelts, 1};
    il::Array<il::int_t> condID{nelts, 1};

    // create a basic 1D mesh ....
    for (il::int_t i = 0; i < xy.size(0); ++i)
    {
        xy(i, 0) = -1. + i * h;
        xy(i, 1) = 0.;
    };

    for (il::int_t i = 0; i < myconn.size(0); ++i)
    {
        myconn(i, 0) = i;
        myconn(i, 1) = i + 1;
    };

    for (il::int_t i = 0; i < nelts; i++)
    {
        for (il::int_t j = 0; j < 2 * (p + 1); j++)
        {
            id_displ(i, j) = i * 2 * (p + 1) + j;
        }
    }

    for (il::int_t i = 0; i < nelts; i++)
    {
        id_press(i, 0) = i;
        id_press(i, 1) = i + 1;
    }

    // Mesh initialization
    hfp2d::Mesh mesh(p, xy, myconn, id_displ, id_press, fracID, matID, condID);
    const il::int_t totalNumDD = mesh.numDisplDofs();
    const il::int_t DDxElem = mesh.numDisplDofsPerElem(); // dof per element

    // Elastic properties initialization
    // double Ep = 1.;                    // Plane strain Young's modulus
    hfp2d::ElasticProperties myelas(1.0e7, 0.);
    // std::cout << "EP :" << myelas.Ep() << "\n";

    // Stress distribution and other quantities
    double epsiP0 = 1.0e-4;         // variation in pore pressure

    double sigmaS = 0.0;            // in situ stress, shear
    double sigmaW = 1.0e6;          // in situ stress, opening

    double initPress = 1.0e6 + epsiP0; // initial pore pressure
    double injectionRate = 1.0e-1;
    double failStress = 1.0e12;
    double maxOpening = 1.0e-3;

    il::int_t finalTimeStep = 30;
    double deltaTime = 1.0e-4;

    const double tolX1 = 1.0e-4;
    const double tolX2 = 1.0e-8;
    const double tolFX = 1.0e-4;

    il::int_t NLiter = 0;
    il::int_t globalIter = 0;
    bool NLSolConv = false;
    bool actSetConv = false;
    const il::int_t NLiterMax = 10;
    const il::int_t globalIterMax = 10;

    SolidEvolution linearCZM(failStress, maxOpening, totalNumDD);

    /////////////////////////////////////////////////////////////////////////////

    //  std::cout << "Number of elements : " << mesh.nelts() << "\n";
    //  std::cout << "------ Assembly: \t";
    //  il::Timer timer{};
    //  timer.start();

    //// CREATION OF DD MATRIX K
    il::Array2D<double> globalK_DD{totalNumDD, totalNumDD};
    globalK_DD = hfp2d::basic_assembly_new(
        mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    //  timer.stop();
    //  std::cout << timer.elapsed() << " s" << "  \n";

    //// CREATION OF DD VECTOR F
    // Set initial distribution of stresses (constant) and initial pressure
    il::Array<double> inSituStr(totalNumDD);

    for (il::int_t i = 0; i < mesh.numElems(); i++)
    {

        inSituStr[mesh.dofDispl(i, 0)] =
            sigmaS; // 1st collocation pt. of elem. i, shear stress
        inSituStr[mesh.dofDispl(i, 1)] =
            sigmaW; // 1st collocation pt. of elem. i, opening stress
        inSituStr[mesh.dofDispl(i, 2)] =
            sigmaS; // 2nd collocation pt. of elem. i, shear stress
        inSituStr[mesh.dofDispl(i, 3)] =
            sigmaW; // 2nd collocation pt. of elem. i, opening stress

    }


    //// CREATION OF GLOBAL SOLUTION VECTOR
    il::Array<double> globalSol(totalNumDD + 1, 0.0); // full solution vector
    il::Array<double> globalDDs(totalNumDD, 0.0);      // DD part of the solution
    double press = 0.0;   // pore pressure solution (both local and global)

    // and their "last time step" counterparts
    il::Array<double> globalSol_old(totalNumDD + 1, 0.0);
    il::Array<double> globalDDs_old(totalNumDD, 0.0);
    double press_old = 0.0;


    //// CREATION OF "ACTIVE LIST" AND RESPECTIVE MATRICES/VECTORS
    il::Array<il::int_t> activeList;
    activeList.reserve(mesh.numElems()); // maximum capacity of the active list
    // is the number of elements

    // Initialize stiffness matrix of active elements
    il::Array2D<double> Kact;
    Kact.reserve(totalNumDD + 1, totalNumDD + 1); // reserve space for full problem
    // but use only the active part to solveDD

    // Initialize force vector of active elements
    il::Array<double> Fact;
    Fact.reserve(totalNumDD + 1); // again similar as before, here we reserve all
    // but use the active part


    //// SET INITIAL SOURCE POSITION AND INITIAL ACTIVATED ELEMENTS
    // Set the initial location of source (middle point, initial personal choice)
    il::int_t sourceElem = mesh.numElems() / 2;
    if (mesh.numElems() % 2 == 0)
    { // even number of elements

        activeList.resize(4);
        activeList[0]=sourceElem-2; //.append(sourceElem - 1);
        activeList[1]=sourceElem-1; //.append(sourceElem - 1);
        activeList[2]=sourceElem;   //.append(sourceElem);
        activeList[3]=sourceElem+1; //.append(sourceElem);

    }
    else
    { // odd number of elements

        activeList.resize(3);
        activeList[0]=sourceElem-1; //.append(sourceElem);
        activeList[1]=sourceElem; //.append(sourceElem);
        activeList[2]=sourceElem+1; //.append(sourceElem);

    }

    il::int_t numActElems = activeList.size();  // number of active elements
    il::int_t numActDispl = numActElems * DDxElem; // # of active dofs

    il::int_t numActElems_old = activeList.size();  // number of active elements
    il::int_t numActDispl_old = numActElems * DDxElem; // # of active dofs


    //// EXTRACTION OF ACTIVE K MATRIX AND F VECTOR
    //  std::cout << "------ Extract: \t";
    //  timer.reset();
    //  timer.start();

    // Resize the matrix Kact and vector Fact to the active set of DD
    Kact.resize(numActDispl, numActDispl);
    Fact.resize(numActDispl);

    // Fill the matrix with the values of active elements
    for (il::int_t i = 0; i < numActElems; i++)
    {
        il::int_t xStart = i * DDxElem;

        for (il::int_t j = 0; j < numActElems; j++)
        {

            // the block to copy starts from xStart,yStart
            il::int_t yStart = j * DDxElem;

            // save the dof values of the corresponding active element
            // in the blocks of the counter
            for (il::int_t k = 0; k < DDxElem; k++)
            {
                for (il::int_t l = 0; l < DDxElem; l++)
                {

                    il::int_t dof1 = mesh.dofDispl(activeList[i], k);
                    il::int_t dof2 = mesh.dofDispl(activeList[j], l);

                    Kact(xStart + k, yStart + l) = globalK_DD(dof1, dof2);
                }
            }
        }

        // save the values of inSituStr in Fact
        Fact[xStart] = inSituStr[mesh.dofDispl(activeList[i], 0)];
        Fact[xStart + 1] =
            inSituStr[mesh.dofDispl(activeList[i], 1)] - initPress;
        Fact[xStart + 2] = inSituStr[mesh.dofDispl(activeList[i], 2)];
        Fact[xStart + 3] =
            inSituStr[mesh.dofDispl(activeList[i], 3)] - initPress;

    }

    // Pressure terms in Kact and Fact are not needed for THE INITIALIZATION STEP!
    //  timer.stop();
    //  std::cout << timer.elapsed() << " s" << "  \n";

    //  il::LU<il::Array2D<double> > lu_decomposition(globalK_DD, il::io, status);
    //  il::Array<double> initialDDSol= lu_decomposition.solve(std::move(inSituStr));

    //// SOLVE INITIAL OPENING PROFILE, PRESSURE AND VOLUME (active and global)
    //  std::cout << "------ Solution:\t";
    //  timer.reset();
    //  timer.start();

    // initial opening
    il::Status status;
    const il::Array<double> initialDDSol =
        il::linearSolve(Kact, Fact, il::io, status);
    status.ok();
    //  timer.stop();
    //  std::cout << timer.elapsed() << " s \n";
    //  std::cout << "---#---\n";

    // initial pore pressure
    press = initPress;

    // initial volume
    double initVolume = 0.0;
    for (il::int_t i = 0; i < numActElems; i++)
    {
        initVolume = initVolume + (initialDDSol[mesh.dofDispl(i, 1)] +
            initialDDSol[mesh.dofDispl(i, 3)]) * h / 2.0;
    }

    double injectedVol = initVolume;
    double injectedVol_old = initVolume;

    // Reconstruct global vectors: from initialDDSol to globalSol and globalDDs
    // Save the DDs
    for (il::int_t i = 0; i < numActElems; i++)
    {
        for (il::int_t j = 0; j < DDxElem; j++)
        {

            globalSol[mesh.dofDispl(activeList[i], j)] =
                initialDDSol[i * DDxElem + j];
            globalDDs[mesh.dofDispl(activeList[i], j)] =
                initialDDSol[i * DDxElem + j];

        }
    }

    // Save the pressure
    // globalSol has length totalNumDD+1 and the pressure is its last value
    globalSol[totalNumDD] = initPress;

    // save "converged" values
    globalSol_old = globalSol;
    globalDDs_old = globalDDs;
    press_old = press;


    ////////////////////////////////////////////////////////////////////////////

    //// INITIALIZATION OF VECTORS FOR ITERATIVE SOLUTION
    // total solution of active elements at this time step
    il::Array<double> actSol;
    actSol.reserve(totalNumDD + 1);
    // total solution of active elements at last time step
    il::Array<double> actSol_old;
    actSol_old.reserve(totalNumDD + 1);
    // iterative total solution of active elements
    il::Array<double> deltaActSol;
    deltaActSol.reserve(totalNumDD + 1);

    il::Array<double> deltaActSol_old;
    deltaActSol_old.reserve(totalNumDD + 1);


    // resize to the active nodes and save initialization step as old
    actSol.resize(numActDispl + 1);
    actSol = initialDDSol;

    actSol_old.resize(numActDispl + 1);
    actSol_old = initialDDSol;


    // initialize global vector of stresses at collocation points
    il::Array<double> globalStressColl(totalNumDD, 0.0);

    // save first active list as old
    il::Array<il::int_t> activeList_old = activeList;


    // global vector of DD at collocation points
    il::Array<double> DDatColl(totalNumDD,0.0);
    // vector of active DD at collocation points (its max size is the totalNumDD)
    il::Array<double> actDDatColl;
    actDDatColl.reserve(totalNumDD);

    //// INITIALIZATION OF ADDITIONAL MATRICES
    // matrix to compute collocation values starting from nodal values
    il::Array2D<double> fetc = from_edge_to_col_dg_full2d_new(mesh);

    // Set max **HISTORICAL** DD in initial active elements to wc
    // note: not current DD but only initial DD!
    // it is done so the cohesive force is zero
    linearCZM.setInitialCrack(activeList,mesh);


    // printing Kact to check if sliding and opening are coupled (and they
    // are not!!)
    //std::cout << std::setprecision(6) ;
//    for(il::int_t i=0; i<numActDispl; i++){
//        for(il::int_t j=0; j<numActDispl; j++){
//            std::cout << std::scientific << std::setprecision(4) << Kact(i,j)
//                      << "\t" ;
//        }
//        std::cout << std::endl;
//    }
//
//    std::cin.get();

    ////////////////////////////////////////////////////////////////////////////

    const il::String horizontalLine =
        "---------------------------------------------------------------";
    std::cout << horizontalLine << std::endl << std::endl;
    std::cout << std::setprecision(6)
              << "    Number of elements: \t" << nelts << std::endl
              << "    Number of time steps: \t" << finalTimeStep << std::endl
              << "    Time step size: \t" << deltaTime << std::endl
              << "    Injection rate: \t" << injectionRate << std::endl
              << "    Volume per time step: \t" << injectionRate*deltaTime
              << std::endl << std::endl;
    std::cout << horizontalLine << std::endl << std::endl;

    //// MULTIPLE TIME STEPS
    for (il::int_t timeStep = 0; timeStep < finalTimeStep; timeStep++)
    {

        std::cout << "   Timestep: \t" << timeStep << std::endl;

        NLSolConv = false;
        actSetConv = false;
        NLiter = 0;
        globalIter = 0;

        /////////////////////////   GLOBAL LOOP   /////////////////////////
        //// Non linear system and activated set of elements must converge
        //globalIter=0;
        while((!NLSolConv || !actSetConv) && (globalIter < globalIterMax))
        {
            NLSolConv = false;
            actSetConv = false;

            // increase iteration counter
            globalIter++;

            std::cout << "   global iteration: \t" << globalIter << std::endl;

            //std::cin.get(); // pause for enter

            //// GATHER - Reconstruct solution on active set
            /// with previous step values
            // indeed, it is assumed that activeList_temp != activeList
            actSol.resize(numActDispl + 1);
            actSol_old.resize(numActDispl + 1);
            deltaActSol.resize(numActDispl + 1);
            deltaActSol_old.resize(numActDispl+1);

            for (il::int_t i = 0; i < numActElems; i++)
            {
                il::int_t xStart = i * DDxElem;
                for (il::int_t j = 0; j < DDxElem; j++)
                {
                    // Note: last converged solution vector values are copied,
                    // for those that were not activated a zero is copied.
                    actSol[xStart + j] =
                        globalSol_old[mesh.dofDispl(activeList[i], j)];

                    deltaActSol[xStart + j]=0.0;
                }
            }
            actSol[numActDispl] = globalSol[totalNumDD]; // for pressure solution
            actSol_old = actSol; // reset actSol_old to old size/values
            deltaActSol_old = deltaActSol;

            //// INITIALIZATION OF Kact AND Fact FOR THE SOLUTION at this
            /// time step
            Kact.resize(numActDispl + 1, numActDispl + 1);
            Fact.resize(numActDispl + 1);

            // Fill the matrix with the values of active elements
            for (il::int_t i = 0; i < numActElems; i++)
            {
                // the block to copy starts from xStart,yStart
                il::int_t xStart = i * DDxElem;

                for (il::int_t j = 0; j < numActElems; j++)
                {
                    il::int_t yStart = j * DDxElem;

                    // save the dof values of the corresponding active element
                    // in the blocks of the counter
                    for (il::int_t k = 0; k < DDxElem; k++)
                    {
                        for (il::int_t l = 0; l < DDxElem; l++)
                        {
                            il::int_t
                                dof1 = mesh.dofDispl(activeList[i], k);
                            il::int_t
                                dof2 = mesh.dofDispl(activeList[j], l);

                            Kact(xStart + k, yStart + l) = globalK_DD(dof1, dof2);
                        }
                    }
                }
            }

            // filling last colum and row of Kact
            for (il::int_t i = 0; i < numActDispl; i = i + 2)
            {
                Kact(i, numActDispl) = 0.0;
                Kact(i + 1, numActDispl) = -1.0;
            }

            for (il::int_t i = 0; i < numActDispl; i = i + 2)
            {
                Kact(numActDispl, i) = 0.0;
                Kact(numActDispl, i + 1) = h / 2.0; //mesh.eltsize(activeList[i]) / 2.0;
            }
            Kact(numActDispl, numActDispl) = 0.0;

            // LOAD Fact
            // 1. Compute displacements at collocation points
            // Calculating DDs at **ALL** collocation points
            DDatColl = il::dot(fetc, globalDDs);

            // Resize the DDs at **ACTIVE** collocation points
            actDDatColl.resize(numActDispl);
            // reload the active displacement discontinuities
            for (il::int_t i = 0; i < numActElems; i++)
            {
                for (il::int_t j = 0; j < DDxElem; j++)
                {

                    actDDatColl[i * DDxElem + j] =
                        DDatColl[mesh.dofDispl(activeList[i], j)];

                }
            }

            // 2. computation of pressure at collocation points from nodes is
            // not required since it is constant

            // 3. update the force vector with the new cohesive stress values
            // computed from the collocation points opening
            for (il::int_t i = 0; i < numActElems; i++)
            {

                il::int_t coll1 = mesh.dofDispl(activeList[i], 1);
                il::int_t coll2 = mesh.dofDispl(activeList[i], 3);

                double opening1 = DDatColl[coll1];
                double opening2 = DDatColl[coll2];

                double cohStress1 =
                    linearCZM.tractionSeparationLaw(opening1, coll1);
                double cohStress2 =
                    linearCZM.tractionSeparationLaw(opening2, coll2);

                // starting location for element in active lists
                il::int_t xStart = i * DDxElem;

                // save the values of inSituStr in Fact
                Fact[xStart] = inSituStr[mesh.dofDispl(activeList[i], 0)];
                Fact[xStart + 1] =
                    inSituStr[mesh.dofDispl(activeList[i], 1)]
                        - press_old
                        + cohStress1;
                Fact[xStart + 2] =
                    inSituStr[mesh.dofDispl(activeList[i], 2)];
                Fact[xStart + 3] =
                    inSituStr[mesh.dofDispl(activeList[i], 3)]
                        - press_old
                        + cohStress2;

            }

            // the initial crack has already an initial volume that we must
            // take into account in the first iteration.
            // Then it will be always Q*deltaT, the additional inj per time step
            Fact[numActDispl] = injectionRate * deltaTime;

            // subtract the part of previous solution, i.e. compute the
            // residual with - 1.0 Kact W + 1.0 Fact
            auto vec1=il::dot(Kact,actSol_old);
            auto vec2=il::dot(Kact,actSol);
            auto vec3=il::dot(Kact,deltaActSol);
            //il::blas(-1.0, Kact, actSol_old, +1.0, il::io, Fact);
            for(il::int_t i=0; i<numActDispl; i++){
                Fact[i] = Fact[i] - vec1[i];
            }
            injectedVol_old=vec1[numActDispl];

            il::Array<double> sol_DD(numActDispl);
            il::Array2D<double> K_DD(numActDispl,numActDispl);
            for (il::int_t i=0; i < numActDispl; i++){
                for(il::int_t j=0; j < numActDispl; j++){
                 K_DD(i,j) = Kact(i,j);
                }
                sol_DD[i]=actSol_old[i];
            }
            auto res_dd = il::dot(K_DD, sol_DD);

            //// SET INITIAL RESIDUAL NORMS
            // --- computing residual of iteration zero using Kact, deltaW = 0
            // means just computing the norm of Fact in new step, with openings
            // from last time step (or active set) and in particular new Q Deltat
//            double normR0 = il::norm(Fact, il::Norm::L2);
            double normR0_DD = normSplit(Fact, 0, numActDispl);
            double normR0_PP = normSplit(Fact, numActDispl, numActDispl + 1);

//            double normSol = il::norm(actSol_old, il::Norm::L2);
            double normDD = normSplit(actSol_old, 0, numActDispl);
            double normPP = normSplit(actSol_old, numActDispl, numActDispl + 1);

            //// LOOP 1 - ITERATIVE SOLUTION OF NON LINEAR SYSTEM OF EQUATIONS
            // NB: the active set is constant in this loop
            // NB: Kact is constant (in shape and values) in this loop
            // NB: pressure is constant everywhere, at collocation points too
            //NLiter = 0;
            while ((!NLSolConv) && (NLiter < NLiterMax))
            {
                // increase iteration counter
                NLiter++;

                //// SOLVE for the correction deltaActSol
                deltaActSol_old = deltaActSol;
                deltaActSol = il::linearSolve(Kact, Fact, il::io, status);
                status.ok();

                // and update solution
                for (il::int_t i = 0; i < numActDispl + 1; i++)
                {
                    actSol[i] = actSol_old[i] + deltaActSol[i];
//                        0.85*(actSol[i] + deltaActSol[i]) + 0.15*actSol[i];
                }

                //// UPDATE GLOBAL SOLUTION
                // first reset solution to zero
//                for (il::int_t i = 0; i < totalNumDD; i++)
//                {
//                    globalSol[i] = 0.0;
//                    globalDDs[i] = 0.0;
//                }

                // save DDs
                for (il::int_t i = 0; i < numActElems; i++)
                {
                    il::int_t ithElem = activeList[i];

                    for (il::int_t j = 0; j < DDxElem; j++)
                    {
                        globalSol[mesh.dofDispl(ithElem, j)] =
                            actSol[i * DDxElem + j];
                        globalDDs[mesh.dofDispl(ithElem, j)] =
                            actSol[i * DDxElem + j];
                    }
                }
                // save pressure
                globalSol[totalNumDD] = actSol[numActDispl];
                press = actSol[numActDispl];



                //// LOAD updated Fact
                // 1. Compute displacements at collocation points
                // Calculating DDs at **ALL** collocation points
                DDatColl = il::dot(fetc, globalDDs);

                // Resize the DDs at **ACTIVE** collocation points
                actDDatColl.resize(numActDispl);
                // reload the active displacement discontinuities
                for (il::int_t i = 0; i < numActElems; i++)
                {
                    for (il::int_t j = 0; j < DDxElem; j++)
                    {

                        actDDatColl[i * DDxElem + j] =
                            DDatColl[mesh.dofDispl(activeList[i], j)];

                    }
                }

                // 2. computation of pressure at collocation points from nodes is
                // not required since it is constant

                // 3. update the force vector with the new cohesive stress values
                // computed from the collocation points opening
                for (il::int_t i = 0; i < numActElems; i++)
                {
                    il::int_t coll1 = mesh.dofDispl(activeList[i], 1);
                    il::int_t coll2 = mesh.dofDispl(activeList[i], 3);

                    double opening1 = DDatColl[coll1];
                    double opening2 = DDatColl[coll2];

                    double cohStress1 =
                        linearCZM.tractionSeparationLaw(opening1, coll1);
                    double cohStress2 =
                        linearCZM.tractionSeparationLaw(opening2, coll2);

                    // starting location for element in active lists
                    il::int_t xStart = i * DDxElem;

                    // save the values of inSituStr in Fact
                    Fact[xStart] = inSituStr[mesh.dofDispl(activeList[i], 0)];
                    Fact[xStart + 1] =
                        inSituStr[mesh.dofDispl(activeList[i], 1)]
                            - press_old
                            + cohStress1;
                    Fact[xStart + 2] =
                        inSituStr[mesh.dofDispl(activeList[i], 2)];
                    Fact[xStart + 3] =
                        inSituStr[mesh.dofDispl(activeList[i], 3)]
                            - press_old
                            + cohStress2;
                }

                Fact[numActDispl]= injectionRate * deltaTime;
                // subtract the part of previous solution, i.e. compute the
                // residual with - 1.0 Kact W + 1.0 Fact
                //il::blas(-1.0, Kact, actSol_old, +1.0, il::io, Fact);
                vec1=il::dot(Kact,actSol_old);
                for(il::int_t i=0; i<numActDispl; i++){
                    Fact[i] = Fact[i] - vec1[i];
                }
                injectedVol=vec1[numActDispl];

                //// CHECK convergence
                // create current residual
                il::Array<double> R = Fact;

                // compute 1.0 Kact dw - 1.0 R
                il::blas(1.0, Kact, deltaActSol, -1.0, il::io, R);

//                double normR = il::norm(R, il::Norm::L2);

//                double normDeltaSol = il::norm(deltaActSol, il::Norm::L2);
//                double normSol = il::norm(actSol, il::Norm::L2);

                il::Array<double> deltadeltaActSol(numActDispl+1);
                for(il::int_t i=0; i<numActDispl+1; i++){
                    deltadeltaActSol[i] = deltaActSol[i]-deltaActSol_old[i];
                }

                // computation of norms of difference of correction
                double normDeltaDD = normSplit(deltadeltaActSol,
                                               0, numActDispl);


                double normDeltaPP = normSplit(deltadeltaActSol,
                                               numActDispl, numActDispl+1);

                double ratioDeltaDD = normDeltaDD/(tolX1 * normDD + tolX2);
                double ratioDeltaPP = normDeltaPP/(tolX1 * normDD + tolX2);


                // Computation of norms and ratios of residuals
                double normResDD = normSplit(R, 0, numActDispl);
                double normResPP = normSplit(R, numActDispl, numActDispl + 1);

                double ratioResDD = normResDD / (normR0_DD * tolFX);
                double ratioResPP = normResPP / (normR0_PP * tolFX);

                NLSolConv = (ratioDeltaDD < 1.0) &&
                            (ratioDeltaPP < 1.0) &&
                            (ratioResDD < 1.0) &&
                            (ratioResPP < 1.0);

                std::cout << std::setprecision(6) << std::scientific
                          << "NL #" << std::left << std::setw(5) <<   NLiter
                          << "   DD: " << std::left << std::setw(10)
                          << ratioDeltaDD
                          << "   PP: " << std::left << std::setw(10)
                          << ratioDeltaPP
                          << "   ResDD: " << std::left << std::setw(10)
                          << ratioResDD
                          << "   ResPP: " << std::left << std::setw(10)
                          << ratioResPP
                          << std::endl;
                //std::cout << "   NL iteration: \t" << NLiter << std::endl;


            }


            //// SCATTER - Save solution to global vector
            // Update global vector

            // first reset solution to zero
//            for (il::int_t i = 0; i < totalNumDD; i++)
//            {
//                globalSol[i] = 0.0;
//                globalDDs[i] = 0.0;
//            }

            // save DDs
            for (il::int_t i = 0; i < numActElems; i++)
            {
                il::int_t ithElem = activeList[i];

                for (il::int_t j = 0; j < DDxElem; j++)
                {
                    globalSol[mesh.dofDispl(ithElem, j)] =
                        actSol[i * DDxElem + j];
                    globalDDs[mesh.dofDispl(ithElem, j)] =
                        actSol[i * DDxElem + j];
                }
            }
            // save pressure
            globalSol[totalNumDD] = actSol[numActDispl];
            press = actSol[numActDispl];

            //// LOOP 2 - CHECK SET OF ACTIVE ELEMENTS
            /// (restart from list of old time step)
            il::int_t numActDispl_temp = numActDispl_old;
            il::int_t numActElems_temp = numActElems_old;
            il::Array<il::int_t> activeList_temp = activeList_old;

            // compute stress at collocation points (using new global values)
            globalStressColl = il::dot(globalK_DD, globalDDs);

            // recheck the number of active elements
            // to do that, let us recheck that both of the collocation points are active
            // Also, we are not checking again the elements that are already active,
            // so we do not start with numActElems=0 but with the previous value of
            // numActElems
            for (il::int_t i = 0; i < mesh.numElems(); i++)
            {

                if (!isElemAlreadyActive(activeList_temp, i))
                {
                    // if the element was not active, let's check again now

                    // here we should check that the STRESS AT BOTH COLLOCATION POINTS
                    // will activate the element/cohesive zone
                    // node 1, component x (sliding),          component y (opening)
                    //         stressFull[mesh.dofDispl(i,0)] stressFull[mesh.dofDispl(i,1)]
                    // node 2, component x (sliding),          component y (opening)
                    //         stressFull[mesh.dofDispl(i,2)] stressFull[mesh.dofDispl(i,3)]
                    if ((globalStressColl[mesh.dofDispl(i, 1)]
                        > linearCZM.getMaxStress(i)) &&
                        (globalStressColl[mesh.dofDispl(i, 3)]
                            > linearCZM.getMaxStress(i)))
                    {

                        activeList_temp.append(i); // add the element to the list
                        // of active ones
                    }

                }
            }

            // new number of active elements
            numActElems_temp = activeList_temp.size();
            // new size of system of equation (for displacements)
            numActDispl_temp = numActElems_temp * DDxElem;

            //// CHECK convergence on set
            actSetConv = (numActElems_temp == numActElems);

            std::cout << "   Old set: \t" << numActElems
                      << "   New set: \t" << numActElems_temp << std::endl;
            // Save the newly computed set as the working one
            activeList = activeList_temp;
            numActElems = numActElems_temp;
            numActDispl = numActDispl_temp;

        }

        for (il::int_t i = 0; i < numActElems; i++)
        {
            initVolume = initVolume + (initialDDSol[mesh.dofDispl(i, 1)] +
                initialDDSol[mesh.dofDispl(i, 3)]) * h / 2.0;
        }

        globalDDs_old = globalDDs;
        globalSol_old = globalSol;
        press_old = press;

        actSol_old=actSol;

        //// OUTPUT
        std::cout << "injected volume = " << injectedVol <<std::endl;
        std::cout << "---End of time step\n" << std::endl;
//            for (il::int_t i = 0; i < Xfull.size(); i++)
//            {
//                std::cout << i << "\t" << Xfull[i] << std::endl;
//            }
        // print to file

        //


    }

    // TODO: make a better output, write down when it is converged, when it
    // is not, when we are updating the set, when the set continues to be
    // fixed, how much volume is injected (versus the one that is computed) and
    // how much pressure is there. Also add the value of the time step and
    // the length of the crack (or tip position..)

std::cout << "End of the analysis" <<
std::endl;

return initialDDSol[numActDispl];
}

/////////////////////////////////////////////////////////////////////////////

bool isElemAlreadyActive(il::Array<il::int_t> activeList, il::int_t element){

    return (std::find(activeList.begin(), activeList.end(), element) != activeList.end());

};

bool checkActOpening(double strShear, double strOpening){

    double strThreshold = 1.0e6;

    return (strOpening > strThreshold);

};

il::StaticArray<double ,2> tractionSeparation(il::int_t i,double u_x, double u_y){

    il::StaticArray<double, 2> cohForces;

    return cohForces;
};

double normDD(il::Array<double> R, il::int_t dd_dofs){

    il::Array<double> DDres(dd_dofs);
    double theNorm;

    for(il::int_t i=0; i<dd_dofs; i++){
        DDres[i] = R[i];
    }

    theNorm = il::norm(DDres,il::Norm::L2);

    return theNorm;
};

double normSplit(il::Array<double> R, il::int_t begin, il::int_t end){

    IL_EXPECT_FAST(begin >= 0); // not strictly necessary but useful to catch
    IL_EXPECT_FAST(end >= 0);   // possible bugs
    IL_EXPECT_FAST(begin-end < R.size());

    il::int_t newSize = end-begin;

    il::Array<double> newR(newSize);
    double theNorm;

    for(il::int_t i=0; i<newSize; i++){
        newR[i] = R[begin+i];
    }

    theNorm = il::norm(newR,il::Norm::L2);

    return theNorm;
};




}
