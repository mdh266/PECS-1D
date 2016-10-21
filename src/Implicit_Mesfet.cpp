/* Main.cpp is a code for solving the Drift Diffusion Poisson Equations
 *  Copyright 2013 Moderate Perfomance Computing Company
 *
 *
 * "Computing at the full speed of Freedom!"
 *
 */

#include "../include/main.hpp"
#include "../include/ddpPrintState.hpp"
#include "../include/ddpComputeDeltaT.hpp"
#include "../include/carrier.hpp"
#include "../include/poisson.hpp"
#include "../include/dopingprofile.hpp"
#include "../include/testUtilities.hpp"


int
main()
{
//  grvy_timer_reset();
//  grvy_timer_begin("Main");

    /////////////////////////////////////////////////////////////////////////////
    // Set Up
    /////////////////////////////////////////////////////////////////////////////

    // Boiler plate stuff
    ddpDomain_type domain;
    ddpGrid_type grid;
    ddpProblemInfo_type problem;
    ddpCarrierConstants_type carrierConstants;
    std::vector<double> timeStamps;
    ddpDenseVector_type ElecFieldDof;
    ddpDenseVector_type EPTS;

    // Read in inputs about domain, grid and problem.
    readInput(domain, grid, problem, carrierConstants, "ddp-input");

    // Doping Profile function pointer
    double (*DP)(const double & x) = &NPlus_N_NPlus;

    // Make the grid for the carrier.
    ddpMakeUniformGrid(problem, domain, grid);
//  ddpPrintGrid(grid);


    // Create the electron object
    ddpCarrier_type electrons;
    electrons.initialize(grid, problem, carrierConstants, "electrons");

    // Set the BC first!
    electrons.setLeft_DirBC(DP);
    electrons.setRight_DirBC(DP);

    // Assemble matrices and vectors for LDG fluxes
    // NOTE: Depends on BC Type
    electrons.assembleLDGMatrices(grid, problem);

    // Aseemble matrices and vectors for Timestepping type
    electrons.setSolver(grid,problem);

    // Now set initial conditions
    electrons.setInitialConditions(grid, DP);


    // create the doping profile object
    ddpDopingProfile_type DopingProfile;
    DopingProfile.setDopingProfile(grid, electrons, DP);

    // create the poisson object
    ddpPoisson_type Poisson;
    Poisson.initialize(grid, problem, carrierConstants);

    // set the boundary conditions
    Poisson.setBias(problem);

    /////////////////////////////////////////////////////////////////////////////
    // Time Stepping
    /////////////////////////////////////////////////////////////////////////////

    //Prepare time stamps
    ddpMakeTimeStamps(problem, timeStamps);
    double tCurrent = timeStamps.at(0);
    double tNext;
    double tDelta;

    // Get and print Initial conditions
    unsigned int timeStampLabel = 0;
    Poisson.solveSystem(electrons,
                        DopingProfile,
                        tCurrent,
                        problem);

    ddpPrintState(problem,
                  electrons,
                  Poisson,
                  timeStampLabel);



    tDelta = 0.01;

    electrons.AssembleLDGSystem(problem,
                                grid,
                                tDelta);

//	cout << "tDelta = " << tDelta << endl;


    for(timeStampLabel = 1; timeStampLabel < timeStamps.size(); ++timeStampLabel)
    {
        tNext = timeStamps.at(timeStampLabel);
        progressBar(timeStampLabel, timeStamps.size(), 50);

        // Get the update Electric field DOFs
        while(tCurrent < tNext)
        {
            // Get Current state of electrons and solve for Poisson

            Poisson.solveSystem(electrons,
                                DopingProfile,
                                tCurrent,
                                problem);

            // Get current state of Poisson and solve for electrons
            ElecFieldDof = Poisson.PoissonState.elecDof;

            electrons.UpdateExplicitDriftTerm(ElecFieldDof);

            electrons.AssembleRHS(tCurrent,
                                  tDelta);

            electrons.BackwardEuler();

            tCurrent += tDelta;
        } // inner while until nex time stamp


        // Output results
        ddpPrintState(problem,
                      electrons,
                      Poisson,
                      timeStampLabel);


        // Update time
    } // outer time stepping while

//  grvy_timer_finalize();
//   grvy_timer_summarize();

    return 0;
} // end main
