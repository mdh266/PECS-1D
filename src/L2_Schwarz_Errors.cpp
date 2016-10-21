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
#include "ddpTimeStepping.hpp"

int
main()
{
//	grvy_timer_init("GRVY TIMING");
//  grvy_timer_reset();
//   grvy_timer_begin("Main");

/////////////////////////////////////////////////////////////////////////
// Set Up
/////////////////////////////////////////////////////////////////////////

//	grvy_timer_begin("Boiler Plate");

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

//	cout << "dx left = " << semiconductior_grid.DeltaxMin << endl;
//	cout << "dx right = " << electrolyte_grid.DeltaxMin << endl;

//  ddpMakeRefinedBoundaryGrid(problem, domain, grid);
//	ddpPrintGrid(semiconductor_grid);
//	ddpPrintGrid(electrolyte_grid);

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

    if(Original == problem.IVP_Type)
    {
        // Set the initial conditions of each carrier
        electrons.setInitialConditions(semiconductor_grid, ElectronsIC);
        holes.setInitialConditions(semiconductor_grid, HolesIC);
        reductants.setInitialConditions(electrolyte_grid, ReductantsIC);
        oxidants.setInitialConditions(electrolyte_grid, OxidantsIC);
    }
    else// if(Continuation == problem.IVP_Type)
    {
        ddpReadInFinalStates(electrons,
                             holes,
                             reductants,
                             oxidants,
                             Poisson);
    }

    // Time variables
    double tol 									  = 1.0e-13;
    bool semiconductor_converged	= false;
    int elec_field_solver_counter = 0;
    int elec_field_solver_period  = 5;
    int converge_counter					= 0;
    int converge_counter_period		= 5000;
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
    tDelta = problem.increase_time_step_factor * tDeltaProposed;
    cout << "tDelta = " << tDelta << endl;
//	cout << "Estimated run time = " << (3.0e-4/60.0)  * (problem.TimeFinal / tDelta)
//			 << " Minutes.\n";


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

//	grvy_timer_end("Boiler Plate");


    /////////////////////////////////////////////////////////////////////////////
    // Time Stepping (Gummel Iteration)
    /////////////////////////////////////////////////////////////////////////////
    for(timeStampLabel=1; timeStampLabel < timeStamps.size(); ++timeStampLabel)
    {
        // tNext = next printing time
        tNext = timeStamps.at(timeStampLabel);
        progressBar(timeStampLabel, timeStamps.size(), 50);

//		cout << "tNext = " << tNext << endl;

        // Keep doing backward Euler until get the the next printing time.
        while(tCurrent < tNext)
        {
            if(!semiconductor_converged)
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

                if((converge_counter % converge_counter_period) == 0)
                {
                    /////////////////////////////////////////////////////////////////////
                    // Check to see if semiconductor converged. If it has recompute
                    // the time step, build and refactorize matrices
                    ////////////////////////////////////////////////////////////////////
//					grvy_timer_begin("Check converged");
                    semiconductor_converged = TimeStepper.check_converged(electrons, holes, tol);
//					grvy_timer_begin("Check converged");

                    if(semiconductor_converged)
                    {
                        //				grvy_timer_begin("Reassemble and refactor");

                        Poisson.getElecFieldVals(EPTS); // over [-1, 1]

                        ddpComputeDeltaT(electrolyte_grid,
                                         EPTS,
                                         reductants,
                                         oxidants,
                                         tDeltaProposed);

                        // increase it by increased_time_step_factor
                        tDelta = 2000 * tDeltaProposed;

//						tDelta = 1.0;
                        cout << "New tDelta = " << tDelta << endl;

                        TimeStepper.AssembleLDGMatrices(electrons,
                                                        holes,
                                                        reductants,
                                                        oxidants,
                                                        problem,
                                                        semiconductor_grid,
                                                        electrolyte_grid,
                                                        tDelta);

                        //				grvy_timer_end("Reassemble and refactor");
                    } // if check_convergedi reassemble
                } // if converge counter
                converge_counter++;

            } // if not_coverged
            else
            {
                ///////////////////////////////////////////////////////////////////////
                // Find Potential and Electric field at every 5th time step update it
                ////////////////////////////////////////////////////////////////////////
                if(0 == (elec_field_solver_counter % elec_field_solver_period) )
                {
                    Poisson.solveSystem(electrons,
                                        holes,
                                        reductants,
                                        oxidants,
                                        DopingProfile,
                                        tCurrent,
                                        problem);

                }
                ////////////////////////////////////////////////////////////////////
                // Assemble the RHSes for the LDG systems
                ////////////////////////////////////////////////////////////////////
                TimeStepper.UpdateExplicitDriftTerms(reductants,
                                                     oxidants,
                                                     Poisson);

                TimeStepper.AssembleBackwardEulerElectrolyteRHS(electrons,
                        holes,
                        reductants,
                        oxidants,
                        tCurrent,
                        tDelta);

                ////////////////////////////////////////////////////////////////////
                // Do parallel backward Euler time stepping
                ////////////////////////////////////////////////////////////////////
                TimeStepper.BackwardEulerSolve(reductants,
                                               oxidants);
            } // end else

            ////////////////////////////////////////////////////////////////////
            // Update current value of time
            ////////////////////////////////////////////////////////////////////

            tCurrent += tDelta;
            elec_field_solver_counter++;

        }// end while < tNext
//	 	grvy_timer_begin("Print State");

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

//		grvy_timer_end("Print State");

    } // for  end time stepping

    //////////////////////////////////////////////////////////////////////
    // Print the dofs of the final states incase want to use for rerunning
    //////////////////////////////////////////////////////////////////////
    /*	ddpPrintFinalDofs(electrons,
    										holes,
    										reductants,
    										oxidants,
    										Poisson);
    */
    int N = electrons.PTSDense.size();
    std::ofstream prt("quad_points.dat");
    for(int i = 0; i < N; i++)
        prt << electrons.PTSDense(i) << "\n";
    prt.close();

    prt.open("quad_weights.dat");
    for(int i = 0; i < N; i++)
        prt << electrons.weightsDense(i) << "\n";
    prt.close();

    cout << "dx = " << semiconductor_grid.DeltaxMin << endl;
    /*
    	ddpDenseVector_type electrons_u = ddpDenseVector_type::Zero(1000);
    	ddpDenseVector_type electrons_q = ddpDenseVector_type::Zero(1000);
    	ddpReadInElectronDOFS(electrons_u, electrons_q);
    	cout << electrons.PTSDense.size() << endl;
    */

    //////////////////////////////////////////////////////////////////////
    // Read in the old density values
    //////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////////////
    // Print the run time performance results
    //////////////////////////////////////////////////////////////////////////
//  grvy_timer_end("Main");
//  grvy_timer_finalize();
//  grvy_timer_summarize();

    return 0;
}
