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

#include "../include/includes.hpp"
#include "../include/main.hpp"
#include "../include/carrier.hpp"
#include "../include/poisson.hpp"
#include "../include/dopingprofile.hpp"
#include "../include/ddpPrintState.hpp"
#include "../include/ddpComputeDeltaT.hpp"
#include "../include/ddpTimeStepping.hpp"

int
main()
{
//  grvy_timer_init("GRVY TIMING");

//  grvy_timer_reset();
//  grvy_timer_begin("Main");

/////////////////////////////////////////////////////////////////////////
// Set Up
/////////////////////////////////////////////////////////////////////////

    // Read in inputs about domain, grid and problem.
    ddpDomain_type domain;
    ddpGrid_type grid;
    ddpProblemInfo_type problem;
    ddpCarrierConstants_type carrierConstants;
    ddpDenseVector_type EPTS;
    ddpDenseVector_type ElecFieldDof;


    readInput(domain, grid, problem, carrierConstants, "ddp-input.ini");

    // Make the grid
//  ddpMakeUniformGrid(problem, domain, grid);
    ddpMakeRightRefinedBoundaryGrid(problem, domain, grid);
    //ddpPrintGrid(grid);

    //Prepare time stamps
    std::vector<double> timeStamps;
    ddpMakeTimeStamps(problem, timeStamps);

    // Make the carriers
    ddpCarrier_type electrons,
                    holes;

    // Initialize all physical constants and allocate memory
    electrons.initialize(grid, problem, carrierConstants, "electrons");
    holes.initialize(grid, problem, carrierConstants, "holes");

    // These are the initial condition function pointers
    double (* ElectronsIC)(const double & x) = &ddpPEC_DonorDopingProfile;
    double (* HolesIC)(const double & x) = &ddpPEC_AcceptorDopingProfile;

    // set the Ohmic contacts at Left BC
    electrons.setLeft_DirBC(ElectronsIC);
    holes.setLeft_DirBC(HolesIC);

    // set the Ohmic contacts at Right BC
//  electrons.setRight_DirBC(ElectronsIC);
//  holes.setRight_DirBC(HolesIC);

    // set the Schottky approximation for the right BC
    electrons.setRight_RobinBC();
    holes.setRight_RobinBC();


    // Assemble the LDG matrices for this problem
    electrons.assembleLDGMatrices(grid, problem);
    holes.assembleLDGMatrices(grid, problem);

    // Assemble the LDG matrices for this problem
    electrons.setSolver(grid, problem);
    holes.setSolver(grid, problem);

    // Set the equilibrium values
    electrons.setEquilibriumDensity(ElectronsIC);
    holes.setEquilibriumDensity(HolesIC);
    double electrons_e =	electrons.carrierProps.EquilibriumDensity;
    double holes_e = holes.carrierProps.EquilibriumDensity;

// set the recombination velocities (already rescaled)
    double v_n = electrons.carrierProps.RecombinationVelocity;
    double v_p = holes.carrierProps.RecombinationVelocity;

    cout << "v_n = " << v_n << endl;
    cout << "v_p = " << v_p << endl;


    // initialize the interface densities
    double electrons_interface = electrons.getInterfaceDensity();
    double holes_interface = holes.getInterfaceDensity();

    ddpRecombination_type SRH_Recombination(electrons, holes);

    // Make the doping profile
    ddpDopingProfile_type  DopingProfile;

    DopingProfile.setDopingProfile(grid,
                                   electrons,
                                   holes,
                                   ElectronsIC,
                                   HolesIC);

    // Make the electric field and potential
    ddpPoisson_type Poisson;

    // Initialize all the physical constants, allocate memory,
    // and build matrices
    Poisson.initialize(grid, problem, carrierConstants);

    // Set the boundary conditions for Poisson's equation
    Poisson.setBias(problem);

    /////////////////////////////////////////////////////////////////////////////
    // Initial Conditions
    /////////////////////////////////////////////////////////////////////////////

    if(Original == problem.IVP_Type)
    {
        // Set the initial conditions of each carrier
        electrons.setInitialConditions(grid, ElectronsIC);
        holes.setInitialConditions(grid, HolesIC);
    }
    else// if(Continuation == problem.IVP_Type)
    {
        ddpReadInFinalStates(electrons,
                             holes,
                             Poisson);
    }

//  grvy_timer_begin("Time Stepping");

    // Time variables
    double tCurrent = timeStamps.at(0);
    double tNext;
    double tDeltaProposed;
    double tDelta;
    unsigned int timeStampLabel = 0;
    progressBar(timeStampLabel, timeStamps.size(), 50);

    // find the initial electric field
    Poisson.solveSystem(electrons,
                        holes,
                        DopingProfile,
                        tCurrent,
                        problem);


    // Print the initial state
    ddpPrintState(problem,
                  electrons,
                  holes,
                  Poisson,
                  timeStampLabel,
                  tCurrent);

    Poisson.getElecFieldVals(EPTS); // over [-1,0]

    // Computes deltaT from the CFL condition
    ddpComputeDeltaT(grid,
                     EPTS,
                     electrons,
                     holes,
                     tDeltaProposed);

    cout << tDeltaProposed << endl;

    tDelta = 0.01;

    electrons.AssembleLDGSystem(problem, grid, tDelta);
    holes.AssembleLDGSystem(problem, grid, tDelta);

    /////////////////////////////////////////////////////////////////////////////
    // Time Stepping (Gummel Iteration)
    /////////////////////////////////////////////////////////////////////////////

    // Do large grain time stepping to compute the compute the next "printable"
    // state
    for(timeStampLabel = 1; timeStampLabel < timeStamps.size(); ++timeStampLabel)
    {
        // tNext = next printing time
        tNext = timeStamps.at(timeStampLabel);
        progressBar(timeStampLabel, timeStamps.size(), 50);

//		cout << "tNext = " << tNext << endl;

        // Keep doing forward Euler until get the the next printing time.
        while(tCurrent < tNext)
        {

            ///////////////////////////////////////////////////////////////////////
            // Find Potential and Electric field at next time step
            ////////////////////////////////////////////////////////////////////////

            Poisson.solveSystem(electrons,
                                holes,
                                DopingProfile,
                                tCurrent,
                                problem);

            ///////////////////////////////////////////////////////////////////////
            // Update RHS terms
            ////////////////////////////////////////////////////////////////////////

            // update the drift terms
            ElecFieldDof = Poisson.PoissonState.elecDof;
            electrons.UpdateExplicitDriftTerm(ElecFieldDof);
            holes.UpdateExplicitDriftTerm(ElecFieldDof);


            // update the interface terms
            electrons_interface = electrons.getInterfaceDensity();
            holes_interface = holes.getInterfaceDensity();
            electrons.setRight_RobinValue( v_n * (electrons_interface - electrons_e) );
            holes.setRight_RobinValue( v_p * (holes_interface - holes_e) );

            // update the recombination terms
            SRH_Recombination.updateRecombination(grid, electrons, holes);

            // assemble the right hand sides
            electrons.AssembleRHS(tCurrent, tDelta);
            holes.AssembleRHS(tCurrent, tDelta);

            ////////////////////////////////////////////////////////////////////
            // Do backward Euler solve
            // ////////////////////////////////////////////////////////////////////
            electrons.BackwardEuler();
            holes.BackwardEuler();

            ////////////////////////////////////////////////////////////////////
            // Update time
            ////////////////////////////////////////////////////////////////////
            tCurrent += tDelta;
        }
        //	grvy_timer_begin("Print State");
        // printing the state
        ddpPrintState(problem,
                      electrons,
                      holes,
                      Poisson,
                      timeStampLabel,
                      tCurrent);
//		grvy_timer_end("Print State");

    } // end time stepping

    // print the final state
    ddpPrintState(problem,
                  electrons,
                  holes,
                  Poisson,
                  timeStampLabel,
                  tCurrent);

    // print out the dofs
    ddpPrintFinalDofs(electrons, holes, Poisson);


    /////////////////////////////////////////////////////////////////////////
    // Print the run time performance results
    //////////////////////////////////////////////////////////////////////////
    cout << "\tExiting Time Stepping" << endl;
// grvy_timer_end("Time Stepping");
// grvy_timer_end("Main");
// grvy_timer_finalize();
// grvy_timer_summarize();

    return 0;
}
