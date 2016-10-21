////////////////////////////////////////////////////////////////////////////
// main_interface.cpp
// By Michael Harmon
//
// This code is for solving the Drift Diffusion Poisson Equations
// to simulate a semiconductor-electrolyte devices with a reactive
// interfaceU ues a local discontinous Galerkin with forward Euler to solve
// the transport equations, while it uses a mixed CG finite element
// discretization to solve for Poisson's equation with the constraint
// that the displacement electric field be continuous across the interface.
///////////////////////////////////////////////////////////////////////////

#include "includes.hpp"
#include "main.hpp"
#include "carrier.hpp"
#include "poisson.hpp"
#include "dopingprofile.hpp"
#include "ddpPrintState.hpp"
#include "ddpComputeDeltaT.hpp"
#include "../include/ddpTimeStepping.hpp"

int
main()
{
    grvy_timer_init("GRVY TIMING");
    grvy_timer_reset();
    grvy_timer_begin("Main");

/////////////////////////////////////////////////////////////////////////
// Set Up
/////////////////////////////////////////////////////////////////////////

    grvy_timer_begin("Boiler Plate");

    // Read in inputs about domain, grid and problem.
    ddpDomain_type total_domain, semiconductor_domain, electrolyte_domain;
    ddpGrid_type semiconductor_grid, electrolyte_grid;
    ddpProblemInfo_type problem;
    ddpCarrierConstants_type carrierConstants;

    readInput(total_domain, semiconductor_grid, problem, carrierConstants, "ddp-input.ini");
    readInput(total_domain, electrolyte_grid, problem, carrierConstants, "ddp-input.ini");

    // make the left subdomain
    semiconductor_domain.LeftEndPoint = total_domain.LeftEndPoint;
    semiconductor_domain.RightEndPoint = 0.0;

    // make the right subdomain
    electrolyte_domain.LeftEndPoint = 0.0;
    electrolyte_domain.RightEndPoint = total_domain.RightEndPoint;


    // Make the semiconductor grid
    ddpMakeRightRefinedBoundaryGrid(problem, semiconductor_domain, semiconductor_grid);


    // make the electrolyte grid
    ddpMakeLeftRefinedBoundaryGrid(problem, electrolyte_domain, electrolyte_grid);
//  ddpMakeRefinedBoundaryGrid(problem, domain, grid);

    //Prepare time stamps
    std::vector<double> timeStamps;
    ddpMakeTimeStamps(problem, timeStamps);

    // Make the carriers
    ddpCarrier_type electrons,
                    holes,
                    reductants,
                    oxidants;

    // Initialize all physical constants and allocate memory
    electrons.initialize(semiconductor_grid, problem, carrierConstants, "electrons");
    holes.initialize(semiconductor_grid, problem, carrierConstants, "holes");
    reductants.initialize(electrolyte_grid, problem, carrierConstants, "reductants");
    oxidants.initialize(electrolyte_grid, problem, carrierConstants, "oxidants");

    // These are the initial condition function pointers
    double (* ElectronsIC)(const double & x) = &ddpPEC_DonorDopingProfile;
    double (* HolesIC)(const double & x) = &ddpPEC_AcceptorDopingProfile;
    double (* ReductantsIC)(const double & x) = &ddpPEC_ReductantsIC;
    double (* OxidantsIC)(const double & x) = &ddpPEC_OxidantsIC;

    // set the Ohmic contacts at Left BC
    electrons.setLeft_DirBC(ElectronsIC);
    holes.setLeft_DirBC(HolesIC);

    // set the Redox values at the Right BC
    reductants.setRight_DirBC(ReductantsIC);
    oxidants.setRight_DirBC(OxidantsIC);

    // Set the boundary conditions at the interface to be of Robin type
    electrons.setRight_RobinBC();
    holes.setRight_RobinBC();
    reductants.setLeft_RobincBC();
    oxidants.setLeft_RobincBC();

    // Assemble the LDG matrices for this problem
    electrons.assembleLDGMatrices(semiconductor_grid, problem);
    holes.assembleLDGMatrices(semiconductor_grid, problem);
    reductants.assembleLDGMatrices(electrolyte_grid, problem);
    oxidants.assembleLDGMatrices(electrolyte_grid, problem);

    // Assemble the LDG matrices for this problem
    electrons.setSolver(semiconductor_grid, problem);
    holes.setSolver(semiconductor_grid, problem);
    reductants.setSolver(electrolyte_grid, problem);
    oxidants.setSolver(electrolyte_grid, problem);

    // Set the equilibrium values
    electrons.setEquilibriumDensity(ElectronsIC);
    holes.setEquilibriumDensity(HolesIC);

    // Make the doping profile
    ddpDopingProfile_type  DopingProfile;
    DopingProfile.setDopingProfile(semiconductor_grid,
                                   electrons,
                                   holes,
                                   ElectronsIC,
                                   HolesIC);

    // Initialize all the physical constants, allocate memory,
    // and build matrices
    ddpPoisson_type Poisson;
    Poisson.initializeWithInterface(semiconductor_grid, electrolyte_grid,
                                    problem, carrierConstants);

    // Set the boundary conditions for Poisson's equation
    Poisson.setBias(problem);

    /////////////////////////////////////////////////////////////////////////////
    // Initial Conditions
    /////////////////////////////////////////////////////////////////////////////



    // Set the initial conditions of each carrier
    electrons.setInitialConditions(semiconductor_grid, ElectronsIC);
    holes.setInitialConditions(semiconductor_grid, HolesIC);
    reductants.setInitialConditions(electrolyte_grid, ReductantsIC);
    oxidants.setInitialConditions(electrolyte_grid, OxidantsIC);

    // Time variables
    double tCurrent = timeStamps.at(0);
    double tNext;
    double tDeltaProposed;
    double tDelta;
    unsigned int timeStampLabel = 0;
    progressBar(timeStampLabel, timeStamps.size(), 50);
    ddpDenseVector_type EPTS;
    ddpTimeStepping_type TimeStepper(electrons,
                                     holes,
                                     reductants,
                                     oxidants,
                                     Poisson,
                                     DopingProfile,
                                     problem,
                                     tCurrent);

    // Print the initial state
    ddpPrintState(problem,
                  electrons,
                  holes,
                  reductants,
                  oxidants,
                  Poisson,
                  timeStampLabel,
                  tCurrent);

    ///////////////////////////////////////////////////////////////////////
    // 	Compute DeltaT from CFL Condition from initial elec field
    ////////////////////////////////////////////////////////////////////////
    Poisson.solveSystem(electrons,
                        holes,
                        reductants,
                        oxidants,
                        DopingProfile,
                        tCurrent,
                        problem);

    Poisson.getElecFieldVals(EPTS); // over [-1, 1]

    // Computes deltaT from the CFL condition
    ddpComputeDeltaT(semiconductor_grid,
                     electrolyte_grid,
                     EPTS,
                     electrons,
                     holes,
                     reductants,
                     oxidants,
                     tDeltaProposed);

    // increase it by increased_time_step_factor
    tDelta = 0.04;


    ///////////////////////////////////////////////////////////////////////
    // Assemble the LDG Systems
    ////////////////////////////////////////////////////////////////////////
    // these correspond to the time independent terms i.e. diffusion terms
    TimeStepper.AssembleLDGMatrices(electrons,
                                    holes,
                                    reductants,
                                    oxidants,
                                    problem,
                                    semiconductor_grid,
                                    electrolyte_grid,
                                    tDelta);

    grvy_timer_end("Boiler Plate");


    /////////////////////////////////////////////////////////////////////////////
    // Time Stepping (Gummel Iteration)
    /////////////////////////////////////////////////////////////////////////////
    for(timeStampLabel=1; timeStampLabel < timeStamps.size(); ++timeStampLabel)
    {
        // tNext = next printing time
        tNext = timeStamps.at(timeStampLabel);
        progressBar(timeStampLabel, timeStamps.size(), 50);


        // Keep doing backward Euler until get the the next printing time.
        while(tCurrent < tNext)
        {
            ///////////////////////////////////////////////////////////////////////
            // Find Potential and Electric field at every time step until the
            // semiconductor converges then every 5th time step update it
            ////////////////////////////////////////////////////////////////////////
            Poisson.solveSystem(electrons,
                                holes,
                                reductants,
                                oxidants,
                                DopingProfile,
                                tCurrent,
                                problem);

            ////////////////////////////////////////////////////////////////////
            // Assemble the RHSes for the LDG systems
            ////////////////////////////////////////////////////////////////////
            TimeStepper.UpdateExplicitDriftTerms(electrons,
                                                 holes,
                                                 reductants,
                                                 oxidants,
                                                 Poisson);


            TimeStepper.AssembleBackwardEulerRHS(electrons,
                                                 holes,
                                                 reductants,
                                                 oxidants,
                                                 semiconductor_grid,
                                                 tCurrent,
                                                 tDelta);

            ////////////////////////////////////////////////////////////////////
            // Do parallel backward Euler time stepping
            ////////////////////////////////////////////////////////////////////
            TimeStepper.BackwardEulerSolve(electrons,
                                           holes,
                                           reductants,
                                           oxidants);

            tCurrent += tDelta;
        }
        grvy_timer_begin("Print State");

        //////////////////////////////////////////////////////////////////////
        // printing the state
        //////////////////////////////////////////////////////////////////////
        ddpPrintState(problem,
                      electrons,
                      holes,
                      reductants,
                      oxidants,
                      Poisson,
                      timeStampLabel,
                      tCurrent);

        grvy_timer_end("Print State");

    } // for  end time stepping

    //////////////////////////////////////////////////////////////////////
    // Print the dofs of the final states incase want to use for rerunning
    //////////////////////////////////////////////////////////////////////
    ddpPrintFinalDofs(electrons,
                      holes,
                      reductants,
                      oxidants,
                      Poisson);

    /////////////////////////////////////////////////////////////////////////
    // Print the run time performance results
    //////////////////////////////////////////////////////////////////////////
    grvy_timer_end("Main");
    grvy_timer_finalize();
    grvy_timer_summarize();

    return 0;
}
