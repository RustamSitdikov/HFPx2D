//
// Created by lorenzo on 10/13/17.
//


#include "ToughnessDominated.h"

// TODO: check with active sets which can close or open (use activeList_old)
// TODO: make it modular
// TODO: make the solidEvolution as general
// TODO: the fracture front can move forward or backward
// TODO: use the solution class
// TODO: align with the mesh class of Brice
// TODO: data from the files
// TODO: homogenization/scaling of the matrices and vector to avoid bad precond.
// TODO: set console prints and proper informations in the console and in output

// TODO: other solvers
// TODO: check the case of summing delta_w on w for the iterative solution?
// TODO: quasi newton/newton method?
// TODO: acceleration by "precomputing" the next time step?
// TODO: acceleration by "precomputing" the new set based on inj. and last step?
// TODO: save crack length as output
// TODO: parallelization of loops (at least)
// TODO: create solver with lagrangian multipliers for the positive opening
// TODO: remove all comments, remove all unused variables
// TODO: reorganize the initialization of the matrices, vectors and parameters


namespace hfp2d
{

// Convention of signs
//      A dw - dp = - A w + p - sigma_0 - sigma_coh( w + dw )

double
ToughnessDominated(int nelts)
{

    int p = 1;
    double totLength = 2.0;
    double h = totLength / (nelts); //  element size
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
        xy(i, 0) = -totLength/2.0 + i * h;
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
    const il::int_t totalNumDD = mesh.numDDDofs();
    const il::int_t DDperElem = mesh.numDDDofsPerElem(); // dof per element

    ////////////////////////////////////////////////////////////

    hfp2d::ElasticProperties myelas(100, 0.); //myelas(20.0e9, 0.);
    double epsiP0 = 0.00001;         // variation in pore pressure

    double sigmaS = 0.0;            // in situ stress, shear
    double sigmaW = 0.0;          // in situ stress, opening

    double initPress = 0.0 + epsiP0; // initial pressure
    double injectionRate = 0.0005; //0.0005
    double failStress = 2.0;
    double maxOpening = 0.002; //failStress/myelas.Ep(); //1.0;

    il::int_t finalTimeStep = 500;
    double deltaTime = 0.01; // secs

    const double tolX1 = 1.0e-4;
    const double tolX2 = 1.0e-8;
    const double tolFX = 1.0e-4;
    const double relaxParam=0.85;

    il::int_t NonLinIter = 0;
    il::int_t FracFrontIter = 0;
    bool NonLinConv = false;
    bool actSetConv = false;
    const il::int_t NonLiniterMax = 200;
    const il::int_t globalIterMax = 200;

    std::string trg_dir{"/home/lorenzo/CLionProjects/HFPx2D_Master_refactor"
                            "/resultsBenchmark/"};

    // Initialize solid evolution, i.e. CZM
    SolidEvolution linearCZM(failStress, maxOpening, totalNumDD);

    // Initialize parameters of the simulation
    simulationParams simParams = initSimParams(deltaTime,
        deltaTime, globalIterMax, NonLiniterMax, tolX1, tolX2, relaxParam);

    // Start timing of the analysis
      il::Timer timer{};
      timer.start();

    //// CREATION OF MATRIX K ON DDs
    il::Array2D<double> globalK_DD{totalNumDD, totalNumDD};
    globalK_DD = hfp2d::basic_assembly_new(
        mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    //// CREATION OF VECTOR F ON DDs
    // Set initial distribution of stresses (constant) and initial pressure
    il::Array<double> inSituStr(totalNumDD);

    for (il::int_t i = 0; i < mesh.numElems(); i++)
    {

        inSituStr[mesh.dofDD(i, 0)] =
            sigmaS; // 1st collocation pt. of elem. i, shear stress
        inSituStr[mesh.dofDD(i, 1)] =
            sigmaW; // 1st collocation pt. of elem. i, opening stress
        inSituStr[mesh.dofDD(i, 2)] =
            sigmaS; // 2nd collocation pt. of elem. i, shear stress
        inSituStr[mesh.dofDD(i, 3)] =
            sigmaW; // 2nd collocation pt. of elem. i, opening stress

    }

    /////////////////////////////////////////////////////////////////////////////

    //// CREATION OF GLOBAL SOLUTION VECTOR (for this and last time step)
    il::Array<double> globalDDs(totalNumDD, 0.0);      // DD part of the solution
    double press = 0.0;   // pore pressure solution (both local and global)

    // and their "last time step" counterparts
    il::Array<double> globalDDs_n(totalNumDD, 0.0);
    double press_n = 0.0;


    //// CREATION OF "ACTIVE LIST" AND RESPECTIVE MATRICES/VECTORS
    il::Array<il::int_t> activeList;
    activeList.reserve(mesh.numElems()); // maximum capacity of the active list
    // is the number of elements

    // Initialize stiffness matrix of active elements
    il::Array2D<double> Kact;
    Kact.reserve(totalNumDD + 1, totalNumDD + 1); // reserve space for full problem
    // but use only the active part to solveDD

    // Initialize force vector of active elements
    il::Array<double> Fact_k;
    Fact_k.reserve(totalNumDD + 1); // again similar as before, here we reserve all
    // but use the active part

    //// SET VECTORS OF DD AT COLLOCATION POINT AND FROM NODE TO COLLOCATION MAT
    // global vector of DD at collocation points
    il::Array<double> DDatColl(totalNumDD,0.0);

    // matrix to compute collocation values starting from nodal values
    il::Array2D<double> fetc = from_edge_to_col_dg_full2d_new(mesh);

    // initialize global vector of stresses at collocation points
    il::Array<double> globalStressColl(totalNumDD, 0.0);
    ////////////////////////////////////////////////////////////////////////////

    //// INITIALIZATION OF VECTORS FOR ITERATIVE SOLUTION
    // total solution of active elements at this time step
    il::Array<double> actSol_n;
    actSol_n.reserve(totalNumDD + 1);
    // total solution of active elements at last time step
    il::Array<double> actSol_n_1;
    actSol_n_1.reserve(totalNumDD + 1);
    // iterative total solution of active elements
    il::Array<double> deltaActSol_k;
    deltaActSol_k.reserve(totalNumDD + 1);
    // iterative total solution of active elements at LAST ITERATION
    il::Array<double> deltaActSol_k_1;
    deltaActSol_k_1.reserve(totalNumDD + 1);

    double maxAperture;

    double normDD;
    double normPP;
    double normDeltaDD;
    double normDeltaPP;
    double ratioDeltaDD;
    double ratioDeltaPP;

    double normR0_DD;
    double normR0_PP;
    double normResDD;
    double normResPP;
    double ratioResDD;
    double ratioResPP;

    /////////////////////////////////////////////////////////////////////////////
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
    std::cout.flush();


    /////////////////////////////////////////////////////////////////////////////
    //// TIME STEP ZERO
    /// SET INITIAL SOURCE POSITION AND INITIAL ACTIVATED ELEMENTS
    // Set the initial location of source (middle point, initial personal choice)
    il::int_t sourceElem = mesh.numElems() / 2;
    if (mesh.numElems() % 2 == 0)
    { // even number of elements

        activeList.resize(2);
        activeList[0]=sourceElem-1;
        activeList[1]=sourceElem;
//        activeList.resize(4);
//        activeList[0]=sourceElem-2;
//        activeList[1]=sourceElem-1;
//        activeList[2]=sourceElem;
//        activeList[3]=sourceElem+1;

    }
    else
    { // odd number of elements

        activeList.resize(3);
        activeList[0]=sourceElem-1;
        activeList[1]=sourceElem;
        activeList[2]=sourceElem+1;

    }
    // TODO: generalize with function taking the coordinate of the source
    // and/or inital list (or number) of active elements

    il::int_t numActElems = activeList.size();  // number of active elements
    il::int_t numActDispl = numActElems * DDperElem; // # of active dofs

    // Set max **HISTORICAL** DD in initial active elements to wc
    // note: not current DD but only initial DD!
    // it is done so the cohesive force is zero
    linearCZM.setInitialCrack(activeList,mesh);

    //// EXTRACTION OF ACTIVE K MATRIX AND F VECTOR
    // Resize the matrix Kact to the active set of DDs
    Kact.resize(numActDispl, numActDispl);
    Kact = extractMatFromList(globalK_DD, activeList, mesh);

    // Resize the vector Fact to the active set of DDs
    Fact_k.resize(numActDispl);

    // TODO: to substitute this with a generic function which loads Fact
    // using the initial active set and globalDDs equal to zero, the cohesive
    // force would be zero
    il::Array<double> tempInSitu = extractVecFromList(inSituStr, activeList,
                                                    mesh);

    for(il::int_t i=0; i<numActElems; i++){

        il::int_t xStart = i * DDperElem;

        Fact_k[xStart    ] = -tempInSitu[xStart];
        Fact_k[xStart + 1] = initPress-tempInSitu[xStart+1];
        Fact_k[xStart + 2] = -tempInSitu[xStart+2];
        Fact_k[xStart + 3] = initPress-tempInSitu[xStart+3];

    }


    //// SOLVE INITIAL OPENING PROFILE, PRESSURE AND VOLUME (active and global)

    // initial opening
    il::Status status;
    const il::Array<double> initialDDSol =
        il::linearSolve(Kact, Fact_k, il::io, status);
    status.ok();

    // initial pore pressure
    press = initPress;

    // Reconstruct global vectors: from initialDDSol to globalDDs
    globalDDs = insertVecInList(initialDDSol, activeList, mesh);

    // Save initial solution as time step n-1 and n
    SolutionK SolutionAtN_1(0.0, 0.0, mesh,
                            globalDDs,
                            initPress,
                            Fact_k,
                            activeList,
                            1, 1, 0.0, 0.0, 0.0, 0.0);

    //saveTimeStep
    SolutionK SolutionAtN;
    SolutionAtN=SolutionAtN_1;

    // TODO: make a save output here for the initialization step

    // resize to the active nodes and save initialization step as old
    actSol_n = SolutionAtN.DDvalues();
    actSol_n.resize(numActDispl + 1);
    actSol_n[numActDispl]=SolutionAtN.pressure();
    // active solution values of last time steps
    actSol_n_1 = actSol_n;

    il::int_t numActDispl_temp = numActDispl;
    il::int_t numActElems_temp = numActElems;
    il::Array<il::int_t> activeList_temp;
    activeList_temp.reserve(mesh.numElems());
    activeList_temp = activeList;

    ////////////////////////////////////////////////////////////////////////////

    //// MULTIPLE TIME STEPS
    // TODO: substitute the for on the time step for the more generic while,
    // which is required when we are dealing with adaptive time steps and
    // analysis which requires to reach a certain time
    for (il::int_t timeStep = 1; timeStep <= finalTimeStep; timeStep++)
    {

        std::cout << "\n\n   Timestep: \t" << timeStep << std::endl;

        // compute here double time = old_timestep + delta_time and
        // double source = delta_time * injection rate
        double time = SolutionAtN_1.time() + deltaTime;
        double source = injectionRate * deltaTime;

        NonLinConv = false;
        actSetConv = false;
        NonLinIter = 0;
        FracFrontIter = 0;

        /////////////////////////   GLOBAL LOOP   /////////////////////////
        //// Non linear system and activated set of elements must converge
        while((!NonLinConv || !actSetConv) && (FracFrontIter < globalIterMax))
        {
            NonLinConv = false;
            actSetConv = false;

            // increase iteration counter
            FracFrontIter++;

            std::cout << "   fracture front iteration: \t" << FracFrontIter <<
                                                                    std::endl;

            //// GATHER - Reconstruct solution on active set
            /// with previous step values
            // indeed, it is assumed that activeList_temp != activeList
            // Save DD from last time step
            actSol_n.resize(numActDispl);
            actSol_n = extractVecFromList(SolutionAtN_1.DDvalues(),activeList,
                                          mesh);
            // add the pressure as well from last time step
            actSol_n.resize(numActDispl+1);
            actSol_n[numActDispl] = SolutionAtN_1.pressure();

            // save the last time step solution with this set of active elements
            actSol_n_1.resize(numActDispl + 1);
            actSol_n_1 = actSol_n; // reset actSol_old to old size/values

            // resize the iterative variables
            deltaActSol_k.resize(numActDispl + 1);  // K step
            deltaActSol_k_1.resize(numActDispl + 1); // K-1 step

            for(il::int_t i=0;i<numActDispl+1; i++)
            {
                deltaActSol_k[i] = 0.0;
                deltaActSol_k[i] = 0.0;
            }


            //// INITIALIZATION OF Kact AND Fact FOR THE SOLUTION at this
            /// time step
            prepareK_tough(activeList, globalK_DD, mesh, il::io, Kact);

            prepareF_tough(source, SolutionAtN_1, fetc, globalK_DD, activeList,
                           mesh, inSituStr, il::io,Fact_k, linearCZM, globalDDs);

            //// SET INITIAL RESIDUAL NORMS
            // --- computing residual of iteration zero using Kact, deltaW = 0
            // means just computing the norm of Fact in new step, with openings
            // from last time step (or active set) and in particular new Q Deltat
            normDD = normSplit(actSol_n_1, 0, numActDispl);
            normPP = normSplit(actSol_n_1, numActDispl, numActDispl + 1);

            normR0_DD = normSplit(Fact_k, 0, numActDispl);
            normR0_PP = normSplit(Fact_k, numActDispl, numActDispl + 1);


            //// LOOP 1 - ITERATIVE SOLUTION OF NON LINEAR SYSTEM OF EQUATIONS
            // NB: the active set is constant in this loop
            // NB: Kact is constant (in shape and values) in this loop
            // NB: pressure is constant everywhere, at collocation points too
            NonLinIter = 0;
            while ((!NonLinConv) && (NonLinIter < NonLiniterMax))
            {
                // increase iteration counter
                NonLinIter++;

                //// SOLVE for the correction deltaActSol
                deltaActSol_k_1 = deltaActSol_k;
                deltaActSol_k = il::linearSolve(Kact, Fact_k, il::io, status);
                status.ok();

                // and update solution
                for (il::int_t i = 0; i < numActDispl+1; i++)
                {
                    //actSol_n[i] = actSol_old[i] + deltaActSol[i];
                    actSol_n[i] = actSol_n_1[i]
                        + (    simParams.relaxParam * deltaActSol_k[i]
                        + (1.0-simParams.relaxParam)* deltaActSol_k_1[i]);
                }

                //// UPDATE GLOBAL SOLUTION
                for(il::int_t i=0;i<mesh.numDDDofs(); i++)
                {
                    globalDDs[i] = 0.0;
                }
                globalDDs = insertVecInList(actSol_n, activeList, mesh);
                // save pressure
                press = actSol_n[numActDispl];

                //// LOAD updated Fact
                prepareF_tough(source, SolutionAtN_1, fetc,
                               globalK_DD, activeList, mesh,inSituStr,
                               il::io, Fact_k, linearCZM, globalDDs);

                //// CHECK convergence
                // create current residual
                il::Array<double> R = Fact_k;

                // compute 1.0 Kact dw - 1.0 R
                il::blas(1.0, Kact, deltaActSol_k, -1.0, il::io, R);

                il::Array<double> deltadeltaActSol(numActDispl+1);
                for(il::int_t i=0; i<numActDispl+1; i++){
                    deltadeltaActSol[i] = deltaActSol_k[i]-deltaActSol_k_1[i];
                }

                // computation of norms of difference of correction
                normDeltaDD = normSplit(deltadeltaActSol, 0, numActDispl);
                normDeltaPP = normSplit(deltadeltaActSol, numActDispl, numActDispl+1);

                ratioDeltaDD = normDeltaDD/(tolX1 * normDD + tolX2);
                ratioDeltaPP = normDeltaPP/(tolX1 * normPP + tolX2);

                // Computation of norms and ratios of residuals
                normResDD = normSplit(R, 0, numActDispl);
                normResPP = normSplit(R, numActDispl, numActDispl + 1);

                ratioResDD = normResDD / (tolX1 * normR0_DD + tolX2);
                ratioResPP = normResPP / (tolX1 * normR0_PP + tolX2);

                NonLinConv = (ratioDeltaDD < 1.0) &&
                            (ratioDeltaPP < 1.0) &&
                            (ratioResDD < 1.0) &&
                            (ratioResPP < 1.0);
//
                std::cout << std::setprecision(6) << std::scientific
                          << "NL #" << std::left << std::setw(5) <<   NonLinIter
                          << "   DD: " << std::left << std::setw(10)
                          << ratioDeltaDD
                          << "   PP: " << std::left << std::setw(10)
                          << ratioDeltaPP
                          << "   ResDD: " << std::left << std::setw(10)
                          << ratioResDD
                          << "   ResPP: " << std::left << std::setw(10)
                          << ratioResPP
                          << std::endl;
                std::cout.flush();
                //std::cout << "   NL iteration: \t" << NonLinIter << std::endl;


            } // end of tangent system solver (exits with globalDDs and press)


            //// SCATTER - Save solution to global vector
//            // DDs and pressure should be already saved
//            for (il::int_t i = 0; i < numActElems; i++)
//            {
//                il::int_t ithElem = activeList[i];
//
//                for (il::int_t j = 0; j < DDperElem; j++)
//                {
//                    globalDDs[mesh.dofDD(ithElem, j)] =
//                        actSol_n[i * DDperElem + j];
//                }
//            }
//            press = actSol_n[numActDispl];

            //// LOOP 2 - CHECK SET OF ACTIVE ELEMENTS
            /// (restart from list of old time step)
            numActDispl_temp = numActDispl;
            numActElems_temp = numActElems;
            activeList_temp = activeList;

            // compute stress at collocation points (using new global values)
            globalStressColl = il::dot(globalK_DD,globalDDs);

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
                    //         stressFull[mesh.dofDispl(i,0)] stressFull[mesh.dofDD(i,1)]
                    // node 2, component x (sliding),          component y (opening)
                    //         stressFull[mesh.dofDispl(i,2)] stressFull[mesh.dofDD(i,3)]
                    // TODO: substitute this "if" with a yielding function
                    // it should be generic to take into account elements which
                    // have been already opened and maybe needs to be "closed"

                    double checkStress1 = -1.0 * (
                        simParams.relaxParam*globalStressColl[mesh.dofDD(i,1)]
                            +(1.0-simParams.relaxParam)
                                *SolutionAtN.stressAtColl(mesh.dofDD(i,1)));

                    double checkStress2 = -1.0 * (
                        simParams.relaxParam*globalStressColl[mesh.dofDD(i,3)]
                            +(1.0-simParams.relaxParam)
                                *SolutionAtN.stressAtColl(mesh.dofDD(i,3)));

//                    if ((-globalStressColl[mesh.dofDD(i, 1)]
//                        +inSituStr[mesh.dofDD(i, 1)]
//                        > linearCZM.getMaxStress(i)) &&
//                        (-globalStressColl[mesh.dofDD(i, 3)]
//                            +inSituStr[mesh.dofDD(i, 3)]
//                            > linearCZM.getMaxStress(i)))
                    if ((checkStress1 + inSituStr[mesh.dofDD(i, 1)]
                        > linearCZM.getMaxStress(i)) &&
                        (checkStress2 + inSituStr[mesh.dofDD(i, 3)]
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
            numActDispl_temp = numActElems_temp * DDperElem;

            //// In the case we want to use the active set of previous time
            /// step, it is possible to run the following two comparisons:
            /// 1. the number of added elements is the same
            /// 2. the elements of one list are permutations of the other list
            /// notice that the first point is required for the second point
            /// to function properly with function std::is_permutation

            //// CHECK convergence on set
            actSetConv = (numActElems_temp == numActElems);

            std::cout << "   Old set: \t" << numActElems
                      << "   New set: \t" << numActElems_temp << std::endl;
//            // Save the actual set as the old one
//            if(!actSetConv)
//            {
//                activeList_old.resize(numActElems);
//                activeList_old = activeList;
//                numActElems_old = numActElems;
//                numActDispl_old = numActDispl;
//            }

            // Save the newly computed set as the working one
            activeList = activeList_temp;
            numActElems = numActElems_temp;
            numActDispl = numActDispl_temp;

            maxAperture=*std::max_element(globalDDs.begin(), globalDDs.end());

            std::cout << "Max aperture: " << maxAperture << std::endl;
        } // End of fracture front loop

        //// SAVE HERE THE DATA

        SolutionAtN.saveTimeStep(timeStep*deltaTime,
                                    deltaTime,
                                    mesh, globalDDs, press, globalStressColl,
                                    activeList, NonLinIter, FracFrontIter,
                                    ratioDeltaDD,ratioDeltaPP,ratioResDD,
                                   ratioResPP);

        SolutionAtN_1 = SolutionAtN;

        globalDDs_n = globalDDs;
        press_n = press;

        actSol_n_1=actSol_n;

        // Computation of crack length
        double crackLength = h * numActElems;

        //if()
        //activeList_old.resize(numActElems);
        //activeList_old = activeList;
        //numActElems_old = numActElems;
        //numActDispl_old = numActDispl;

        //// OUTPUT
        std::string of_name; // = std::string{"Test"} + std::to_string
            //(timeStep*deltaTime) +
            //std::string{".txt"};

        char filename[50];

        std::sprintf(filename,"Test%010.6f.txt",(timeStep*deltaTime));

        of_name = std::string(filename);

        std::string f_path = trg_dir + of_name;
        const char *format1 = "%.16g %10";
        const char *format6 = "%6.0d";
        const char *format2 = "N. iteration: %i\n\n";
        const char *format3 = "Current time step: %.5g\n\n";
        const char *format4 = "Current time:\n%2.5g\n\n";
        //const char *format5 = "Slippage length:\n%2.5g";

        FILE *of = std::fopen(f_path.c_str(), "w");
        std::fprintf(of, format2, NonLinIter, 1);
        std::fprintf(of, format3, deltaTime, 1);
        std::fprintf(of, format4, deltaTime*timeStep);
        //std::fprintf(of, format5, SolutionAtTj.slippagezone);

        std::fputs("\n\n******* Pressure *******\n", of);
//        std::fprintf(of, format1, press);
        std::fprintf(of, format1, SolutionAtN.pressure());

        std::fprintf(of, "\n\n******* Total sliding *******\n");
        for (int j = 0; j < totalNumDD; j=j+2) {
            std::fprintf(of, format1, SolutionAtN.DDvalues(j));
//            std::fprintf(of, format1, globalDDs[j]);
        }

        std::fprintf(of, "\n\n******* Total opening *******\n");
        for (int j = 1; j < totalNumDD; j=j+2) {
            std::fprintf(of, format1, SolutionAtN.DDvalues(j));
//            std::fprintf(of, format1, globalDDs[j]);
        }

        std::fprintf(of, "\n\n******* Stress XX *******\n");
        for (int j = 0; j < totalNumDD; j=j+2) {
//            std::fprintf(of, format1, globalStressColl[j]);
            std::fprintf(of, format1, SolutionAtN.stressAtColl(j));
        }

        std::fprintf(of, "\n\n******* Stress YY *******\n");
        for (int j = 1; j < totalNumDD; j=j+2) {
//            std::fprintf(of, format1, globalStressColl[j]);
            std::fprintf(of, format1, SolutionAtN.stressAtColl(j));
        }

        std::fprintf(of, "\n\n******* Crack length *******\n");
//        for (int j = 1; j < totalNumDD; j=j+2) {
            std::fprintf(of, format1, crackLength);
//        }


        // TODO> add the volume computation

        std::fprintf(of, "\n\n######################################\n\n");

        std::fclose(of);

        if(numActElems == nelts){
            std::cerr << "no more elements!!" << std::endl;
            abort();
        }


//        std::cout << "injected volume = " << injectedVol <<std::endl;
//        std::cout << "---End of time step\n" << std::endl;
//            for (il::int_t i = 0; i < Xfull.size(); i++)
//            {
//                std::cout << i << "\t" << Xfull[i] << std::endl;
//            }
        // print to file

        //

//        std::ofstream outputDD("outputDD.txt", std::ios_base::app |
//            std::ios_base::out);
//        for(il::int_t i=1; i<totalNumDD; i=i+2){
//            outputDD << globalDDs[i] << "\n";
//        }
//        outputDD.close();
//
//        std::ofstream foutlc;
//        foutlc.open("cracklength.txt");
//        for(int lca=0;lca<break_time;++lca){
//            foutlc<<l_c[lca]<<"\n";
//        }
//        foutlc.close();
//
//
//        std::ofstream foutit;
//        foutit.open("iteration.txt");
//        for(int itera=0;itera<break_time;++itera){
//            foutit<<mvalue[itera]<<"\n";
//        }
//        foutit.close();
//
//        il::Array<double> xlist{2*mesh.nelts(),0.};
//        hfp2d::get_xlist_col(xlist,mesh);
//        std::ofstream fout;
//        fout.open("outputcn1.txt");
//
//
//        for(int qq=0;qq<break_time;++qq){//qq<widthlist.size(0)
//            for(int qqq=0;qqq<widthlist.size(1);++qqq){
//                fout<<widthlist(qq,qqq)<<"\t";
//            }
//            fout<<"\n";
//        }
//        fout.close();


    }

    // TODO: make a better output, write down when it is converged, when it
    // is not, when we are updating the set, when the set continues to be
    // fixed, how much volume is injected (versus the one that is computed) and
    // how much pressure is there. Also add the value of the time step and
    // the length of the crack (or tip position..)

std::cout << "End of the analysis" << std::endl;

      timer.stop();
      std::cout << timer.elapsed() << " s" << "  \n";

return maxAperture;
}

/////////////////////////////////////////////////////////////////////////////

bool isElemAlreadyActive(il::Array<il::int_t> activeList, il::int_t element){

    return (std::find(activeList.begin(), activeList.end(), element) != activeList.end());

};

bool checkActOpening(double strShear, double strOpening){

    double strThreshold = 1.0e6;

    return (strOpening > strThreshold);

}

il::StaticArray<double ,2> tractionSeparation(il::int_t i,double u_x, double u_y){

    il::StaticArray<double, 2> cohForces;

    return cohForces;
}

double normDD(il::Array<double> R, il::int_t dd_dofs){

    il::Array<double> DDres(dd_dofs);
    double theNorm;

    for(il::int_t i=0; i<dd_dofs; i++){
        DDres[i] = R[i];
    }

    theNorm = il::norm(DDres,il::Norm::L2);

    return theNorm;
}

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
}



il::Array<double> extractVecFromList(const il::Array<double> &fullVec,
                                     const il::Array<il::int_t> &elemList,
                                     const Mesh &mesh)
{

    // compute the size of the vector
    const il::int_t DDperElem = mesh.numDDDofsPerElem();
    const il::int_t numActElem = elemList.size();
    const il::int_t numActDispl = DDperElem * numActElem;

    // create final vector
    il::Array<double> F(numActDispl);

    for (il::int_t i = 0; i < numActElem; i++)
    {
        il::int_t xStart = i * DDperElem;

        // save the values of the fullVec in F
        F[xStart]     = fullVec[mesh.dofDD(elemList[i], 0)];
        F[xStart + 1] = fullVec[mesh.dofDD(elemList[i], 1)];
        F[xStart + 2] = fullVec[mesh.dofDD(elemList[i], 2)];
        F[xStart + 3] = fullVec[mesh.dofDD(elemList[i], 3)];

    }

    return F;

}

il::Array<double> insertVecInList(const il::Array<double> &smallVec,
                                  const il::Array<il::int_t> &elemList,
                                  const Mesh &mesh)
{
    // compute the size of the vector
    const il::int_t DDperElem = mesh.numDDDofsPerElem();
    const il::int_t numActElems = elemList.size();
    const il::int_t numActDispl = DDperElem * numActElems;

    il::Array<double> fullVec(mesh.numDDDofs(),0.0);

    for (il::int_t i = 0; i < numActElems; i++)
    {
        il::int_t ithElem = elemList[i];

        for (il::int_t j = 0; j < DDperElem; j++)
        {
            fullVec[mesh.dofDD(ithElem, j)] =
                smallVec[i * DDperElem + j];
        }
    }

    return fullVec;
};


il::Array<double> extractFfromList(il::Array<double> &inSituStr,
                                   double pressure,
                                   il::Array<il::int_t> &elemList,
                                   Mesh &mesh)
{

    const il::int_t DDperElem = mesh.numDDDofsPerElem();
    const il::int_t numActElem = elemList.size();
    const il::int_t numActDispl = DDperElem * numActElem;

    il::Array<double> F(numActDispl);

    for (il::int_t i = 0; i < numActElem; i++)
    {
        il::int_t xStart = i * DDperElem;

        // RHS is pressure - inSituStr by convention

        // save the values of inSituStr in Fact
        F[xStart]     =          - inSituStr[mesh.dofDD(elemList[i], 0)];
        F[xStart + 1] = pressure - inSituStr[mesh.dofDD(elemList[i], 1)];
        F[xStart + 2] =          - inSituStr[mesh.dofDD(elemList[i], 2)];
        F[xStart + 3] = pressure - inSituStr[mesh.dofDD(elemList[i], 3)];

    }

    return F;

}

//void gatherSolution(il::Array<double> globalSol, )


il::Array2D<double> extractMatFromList(const il::Array2D<double> &fullMat,
                                       const il::Array<il::int_t> &elemList,
                                       const Mesh &mesh)
{

    il::int_t numElem = elemList.size();
    il::int_t DDperElem = mesh.numDDDofsPerElem();
    il::int_t matSize = numElem * DDperElem;

    il::Array2D<double> subMat{matSize, matSize};
//
    il::int_t xStart, yStart, dof1, dof2;

    for (il::int_t i = 0; i < numElem; i++)
    {
        // the block to copy starts from xStart,yStart
        xStart = i * DDperElem;

        for (il::int_t j = 0; j < numElem; j++)
        {
            yStart = j * DDperElem;

            // save the dof values of the corresponding active element
            // in the blocks of the counter
            for (il::int_t k = 0; k < DDperElem; k++)
            {
                for (il::int_t l = 0; l < DDperElem; l++)
                {

                    dof1 = mesh.dofDD(elemList[i], k);

                    dof2 = mesh.dofDD(elemList[j], l);

                    subMat(xStart + k, yStart + l) = fullMat(dof1, dof2);
                }
            }
        }
    }

    return subMat;
}

il::Array<double> takeVecFromList(il::Array<double> const &fullVec,
                                  il::Array<il::int_t> const &elemList,
                                  Mesh const &mesh)
{
    il::int_t numElem = elemList.size();
    il::int_t DDperElem = mesh.numDDDofsPerElem();
    il::int_t matSize = numElem * DDperElem;

    il::Array<double> subVec{matSize};

// Fill the matrix with the values of active elements
    for (il::int_t i = 0; i < numElem; i++)
    {
        il::int_t xStart = i * DDperElem;

        // save the values of inSituStr in Fact
        subVec[xStart]     = fullVec[mesh.dofDD(elemList[i], 0)];
        subVec[xStart + 1] = fullVec[mesh.dofDD(elemList[i], 1)];
        subVec[xStart + 2] = fullVec[mesh.dofDD(elemList[i], 2)];
        subVec[xStart + 3] = fullVec[mesh.dofDD(elemList[i], 3)];

    }

    return subVec;

}

il::Array2D<double> takeMatFromOldToNewList(il::Array2D<double> const &fullMat,
                                            il::Array<il::int_t> const &oldList,
                                            il::Array<il::int_t> const &newList,
                                            Mesh const &mesh)
{
    // TODO: check that the function works correctly

    il::int_t oldNumElem = oldList.size();
    il::int_t newNumElem = newList.size();
    il::int_t DDperElem = mesh.numDDDofsPerElem();
    il::int_t matSize = newNumElem * DDperElem;

    il::Array2D<double> subMat{matSize, matSize};
//
    il::int_t xStart, yStart, dof1, dof2;

    for (il::int_t i = oldNumElem; i < newNumElem; i++)
    {
        // the block to copy starts from xStart,yStart
        xStart = i * DDperElem;

        for (il::int_t j = 0; j < newNumElem; j++)
        {
            yStart = j * DDperElem;

            // save the dof values of the corresponding active element
            // in the blocks of the counter
            for (il::int_t k = 0; k < DDperElem; k++)
            {
                for (il::int_t l = 0; l < DDperElem; l++)
                {

                    dof1 = mesh.dofDD(newList[i], k);

                    dof2 = mesh.dofDD(newList[j], l);

                    subMat(xStart + k, yStart + l) = fullMat(dof1, dof2);
                }
            }
        }
    }

    return subMat;
};


il::Array<double> takeVecFromOldToNewList(il::Array<double> const &fullVec,
                                          il::Array<il::int_t> const &oldList,
                                          il::Array<il::int_t> const &newList,
                                          Mesh const &mesh)
{
    // TODO: check that the function works correctly

    il::int_t oldNumElem = oldList.size();
    il::int_t newNumElem = newList.size();
    il::int_t DDperElem = mesh.numDDDofsPerElem();
    il::int_t matSize = newNumElem * DDperElem;

    il::Array<double> subVec{matSize};

// Fill the matrix with the values of active elements
    for (il::int_t i = oldNumElem; i < newNumElem; i++)
    {
        il::int_t xStart = i * DDperElem;

        // save the values of inSituStr in Fact
        subVec[xStart]     = fullVec[mesh.dofDD(newList[i], 0)];
        subVec[xStart + 1] = fullVec[mesh.dofDD(newList[i], 1)];
        subVec[xStart + 2] = fullVec[mesh.dofDD(newList[i], 2)];
        subVec[xStart + 3] = fullVec[mesh.dofDD(newList[i], 3)];

    }

    return subVec;
};

void prepareK_tough(const il::Array<il::int_t> &activeList,
                    const il::Array2D<double> &globalK_DD,
                    Mesh &mesh,
                    il::io_t,
                    il::Array2D<double> &Kact)
{
    const il::int_t numActElems = activeList.size();
    const il::int_t DDperElem = mesh.numDDDofsPerElem();
    const il::int_t numActDispl = numActElems * DDperElem;

    Kact.resize(numActDispl, numActDispl);
    Kact = extractMatFromList(globalK_DD, activeList, mesh);

    // filling last colum and row of Kact
    Kact.resize(numActDispl + 1, numActDispl + 1);

    for (il::int_t i = 0; i < numActElems; i++)
    {
        il::int_t xStart = i*DDperElem;
        Kact(xStart    , numActDispl) =  0.0;
        Kact(xStart + 1, numActDispl) = -1.0;
        Kact(xStart + 2, numActDispl) =  0.0;
        Kact(xStart + 3, numActDispl) = -1.0;
    }

    for (il::int_t i = 0; i < numActElems; i++)
    {
        il::int_t xStart = i*DDperElem;
        Kact(numActDispl, xStart    ) = 0.0;
        Kact(numActDispl, xStart + 1) = mesh.elt_size(activeList[i]) / 2.0;
        Kact(numActDispl, xStart + 2) = 0.0;
        Kact(numActDispl, xStart + 3) = mesh.elt_size(activeList[i]) / 2.0;
    }
    Kact(numActDispl, numActDispl) = 0.0;

}

void prepareF_tough(const double source,
                    const SolutionK &solutionAtN,
                    const il::Array2D<double> &fetc,
                    const il::Array2D<double> &globalK_DD,
                    const il::Array<il::int_t> &activeList,
                    const Mesh &mesh,
                    const il::Array<double> &inSituStress,
                    il::io_t,
                    il::Array<double> &Fact_k,
                    SolidEvolution &CZM,
                    il::Array<double> &globalDDs)
{
    const il::int_t numActElems = activeList.size();
    const il::int_t DDperElem = mesh.numDDDofsPerElem();
    const il::int_t numActDispl = numActElems * DDperElem;

    Fact_k.resize(numActDispl);

    // LOAD Fact
    // 1. Compute displacements at collocation points
    // Calculating DDs at **ALL** collocation points
    il::Array<double> DDatColl = il::dot(fetc, globalDDs); // we
    // use old
    // values
    // because we are still in the initialization of the timestep

    // 2. computation of pressure at collocation points from nodes is
    // not required since it is constant

    // 3. update the force vector with the new cohesive stress values
    // computed from the collocation points opening
    for (il::int_t i = 0; i < numActElems; i++)
    {

        il::int_t coll1 = mesh.dofDD(activeList[i], 1);
        il::int_t coll2 = mesh.dofDD(activeList[i], 3);

        double opening1 = DDatColl[coll1];
        double opening2 = DDatColl[coll2];

        double cohStress1 =
            CZM.tractionSeparationLaw(opening1, coll1);
        double cohStress2 =
            CZM.tractionSeparationLaw(opening2, coll2);

        // starting location for element in active lists
        il::int_t xStart = i * DDperElem;

        // save the values of inSituStr in Fact
        Fact_k[xStart]     = - inSituStress[mesh.dofDD(activeList[i], 0)];

        Fact_k[xStart + 1] = - inSituStress[mesh.dofDD(activeList[i], 1)]
                             - cohStress1;
        Fact_k[xStart + 2] = - inSituStress[mesh.dofDD(activeList[i], 2)];

        Fact_k[xStart + 3] = - inSituStress[mesh.dofDD(activeList[i], 3)]
                             - cohStress2;

    }

    // TODO: check that it is working properly
    // recover the elastic matrix with correct number of active elements from
    // last step
    //
    // or
    //
    // compute the full vector of stresses K
    il::Array2D<double> K = extractMatFromList(globalK_DD,
                                               activeList, //solutionAtN.activeList(),
                                               mesh); //solutionAtN.mesh());

    il::Array<double> ddvals = extractVecFromList(solutionAtN.DDvalues(),
                                                  activeList,
                                                  mesh);

    // subtract the contribution due to DD from previous solution,
    // compute [-K] w_n and subtract it from Fact
//    auto oldContribution = il::dot(K,ddvals);
//    for(il::int_t i=0; i < numActElems; i++)
//    {
//        il::int_t xStart = i * DDperElem;
//
//        Fact_k[xStart    ] -= oldContribution[mesh.dofDD(activeList[i], 0)];
//        Fact_k[xStart + 1] -= oldContribution[mesh.dofDD(activeList[i], 1)];
//        Fact_k[xStart + 2] -= oldContribution[mesh.dofDD(activeList[i], 2)];
//        Fact_k[xStart + 3] -= oldContribution[mesh.dofDD(activeList[i], 3)];
//    }


    il::blas(-1.0, K, ddvals,1.0, il::io, Fact_k);

    // also we have to add pressure contribution in the active nodes
    // of last time step
    for(il::int_t i=0; i < solutionAtN.activeList().size(); i++)
    {
        il::int_t xStart = i * DDperElem;

        Fact_k[xStart + 1] += solutionAtN.pressure();
        Fact_k[xStart + 3] += solutionAtN.pressure();
    }

    // Source value is Q * DeltaT
    Fact_k.resize(numActDispl + 1);
    Fact_k[numActDispl] = source;

} // end prepareF_tough


//bool nonLinSolver(SolutionK &SolutionAtN, SolutionK &SolutionAtN_1){


simulationParams initSimParams(double minDeltaTime,
                               double maxDeltaTime,
                               il::int_t fracfrontMaxIter,
                               il::int_t nonlinMaxIter,
                               double tolX1,
                               double tolX2,
                               double relaxParam){

    simulationParams simParams = {};

    simParams.minDeltat = minDeltaTime;
    simParams.maxDeltat = maxDeltaTime;
    simParams.ffMaxIter = fracfrontMaxIter;
    simParams.nlMaxIter = nonlinMaxIter;
    simParams.tolX1 = tolX1;
    simParams.tolX2 = tolX2;
    simParams.relaxParam = relaxParam;

    return simParams;
};

/*void initialSolution(Mesh &mesh,
                     double &initPress,
                     il::Array<double> &inSituStr,
                     il::Array2D<double> &globalK_DD,
                     il::io_t,
                     SolidEvolution &linearCZM,
                     il::Array<il::int_t> &activeList,
                     il::Array2D<double> &Kact,
                     il::Array<double> &Fact,
                     il::Array<double> &initialDDSol){

    /// SET INITIAL SOURCE POSITION AND INITIAL ACTIVATED ELEMENTS
    // Set the initial location of source (middle point, initial personal choice)
    il::int_t sourceElem = mesh.numElems() / 2;
    if (mesh.numElems() % 2 == 0)
    { // even number of elements

        activeList.resize(4);
        activeList[0]=sourceElem-1;
        activeList[1]=sourceElem;
        activeList[2]=sourceElem+1;
        activeList[3]=sourceElem+2;

    }
    else
    { // odd number of elements

        activeList.resize(3);
        activeList[0]=sourceElem-1;
        activeList[1]=sourceElem;
        activeList[2]=sourceElem+1;

    }
    // TODO: generalize with function taking the coordinate of the source
    // and/or inital list (or number) of active elements

    il::int_t numActElems = activeList.size();  // number of active elements
    il::int_t DDperElem = mesh.numDDDofsPerElem(); // number of DD dofs per elem
    il::int_t numActDispl = numActElems * DDperElem; // # of active dofs

    // Set max **HISTORICAL** DD in initial active elements to wc
    // note: not current DD but only initial DD!
    // it is done so the cohesive force is zero
    linearCZM.setInitialCrack(activeList,mesh);

    //// EXTRACTION OF ACTIVE K MATRIX AND F VECTOR
    // Resize the matrix Kact to the active set of DDs
    Kact.resize(numActDispl, numActDispl);
    Kact = extractMatFromList(globalK_DD, activeList, mesh);

    // Resize the vector Fact to the active set of DDs
    Fact.resize(numActDispl);

    // TODO: to substitute this with a generic function which loads Fact
    // using the initial active set and globalDDs equal to zero, the cohesive
    // force would be zero
    il::Array<double> tempInSitu = extractVecFromList(inSituStr, activeList,
                                                      mesh);

    for(il::int_t i=0; i<numActElems; i++){

        il::int_t xStart = i * DDperElem;

        Fact[xStart    ] = -tempInSitu[xStart];
        Fact[xStart + 1] = initPress-tempInSitu[xStart+1];
        Fact[xStart + 2] = -tempInSitu[xStart+2];
        Fact[xStart + 3] = initPress-tempInSitu[xStart+3];

    }


    //// SOLVE INITIAL OPENING PROFILE, PRESSURE AND VOLUME (active and global)

    // initial opening
    il::Status status;
    const il::Array<double> initialDDSol =
        il::linearSolve(Kact, Fact, il::io, status);
    status.ok();

}*/




} // end namespace




