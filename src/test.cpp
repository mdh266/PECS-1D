#include "../include/test.hpp"

int
ddpTest_type::
initializeTest(const ddpTestChargeCarrier_type & carrierType,
               const ddpGrid_type & inputgrid,
               const ddpProblemInfo_type & inputproblem,
               const ddpCarrierConstants_type & testConstants)
{
    //This just needs to be set before initializing the carriers
    // set the domain, difusivity, mobility, etc. for the testing
    // but does not set anything specifically related tests.
    // use the setTestCarrier or setTestPoisson functions for those
    // purposes.

    testgrid = inputgrid;
    testproblem = inputproblem;

    // initialize the test carrier and solution carrier
    if(carrierType == NegativeCharge)
    {
        testCarrier.initialize(testgrid, testproblem, testConstants, "electrons");
        trueCarrier.initialize(testgrid, testproblem, testConstants, "electrons");
    }
    else if(carrierType == PositiveCharge)
    {
        testCarrier.initialize(testgrid, testproblem, testConstants, "holes");
        trueCarrier.initialize(testgrid, testproblem, testConstants, "holes");
    }
    else
    {
        cout << "No carrier of that type" << endl;
        assert(false);
    }

    testCarrier1.initialize(testgrid, testproblem, testConstants, "electrons");
    testCarrier2.initialize(testgrid, testproblem, testConstants, "holes");


    // initialize poisson
    testPoisson.initialize(testgrid, testproblem, testConstants);

    // Set the bijections
    DGForward = testPoisson.Bijections.DGForward;
    MXForward = testPoisson.Bijections.MXForward;


    // Convert mobility and diffusivity to identity matrix!
    testCarrier1.carrierProps.CurrentScale = 1.0;
    testCarrier2.carrierProps.CurrentScale = 1.0;
    testCarrier.carrierProps.CurrentScale = 1.0;

    //Set it so that coupling is off.
    testPoisson.PoissonProps.coupledOrNot = 0.0;

    // Get the maximum mobility and diffusivity values!
    testCarrier1.carrierProps.MaxMobility = 1.0;
    testCarrier1.carrierProps.MaxDiffusivity = 1.0;

    testCarrier2.carrierProps.MaxMobility = 1.0;
    testCarrier2.carrierProps.MaxDiffusivity = 1.0;

    testCarrier.carrierProps.MaxMobility = 1.0;
    testCarrier.carrierProps.MaxDiffusivity = 1.0;


    // basically sets lamba = 1 in poissons equation
    testPoisson.remakeProperties(testproblem, testgrid);

    // set the vandemondematrices and the datastructures needed
    // to project functions on the DG and CG basis.
    globalVandeMondeDG = testPoisson.VandeMondeMatrices.globalVandeMondeDG;

    globalVandeMondeMX = testPoisson.VandeMondeMatrices.globalVandeMondeMX;

    PTForward = testPoisson.Bijections.PTForward;

    sparseWeights = testPoisson.weightsSparse;

    MassU = testCarrier.carrierProps.MassU;
    A00 = testPoisson.PoissonProps.A00;

    BIJFlagDG = DG;
    BIJFlagMX = Mixed;

    return 0;
}

int
ddpTest_type::
ReinitializeTestWithScaling(const ddpTestChargeCarrier_type & carrierType,
                            const ddpGrid_type & inputgrid,
                            const ddpProblemInfo_type & inputproblem,
                            const ddpCarrierConstants_type & testConstants)
{
    //This just needs to be set before initializing the carriers
    // set the domain, difusivity, mobility, etc. for the testing
    // but does not set anything specifically related tests.
    // use the setTestCarrier or setTestPoisson functions for those
    // purposes.

    testgrid = inputgrid;
    testproblem = inputproblem;

    // initialize the test carrier and solution carrier
    if(carrierType == NegativeCharge)
    {
        testCarrier.initialize(testgrid, testproblem, testConstants, "electrons");
        trueCarrier.initialize(testgrid, testproblem, testConstants, "electrons");
    }
    else if(carrierType == PositiveCharge)
    {
        testCarrier.initialize(testgrid, testproblem, testConstants, "holes");
        trueCarrier.initialize(testgrid, testproblem, testConstants, "holes");
    }
    else
    {
        cout << "No carrier of that type" << endl;
        assert(false);
    }

    testCarrier1.initialize(testgrid, testproblem, testConstants, "electrons");
    testCarrier2.initialize(testgrid, testproblem, testConstants, "holes");

    // initialize poisson
    testPoisson.initialize(testgrid, testproblem, testConstants);

    // Set the bijections
    DGForward = testPoisson.Bijections.DGForward;
    MXForward = testPoisson.Bijections.MXForward;


    // Convert mobility and diffusivity to identity matrix!
    testCarrier1.carrierProps.CurrentScale = 1.0;
    testCarrier2.carrierProps.CurrentScale = 1.0;
//	testCarrier.carrierProps.CurrentScale = 1.0;

    cout << testCarrier.carrierProps.CurrentScale << endl;

    //Set it so that coupling is off.
    testPoisson.PoissonProps.coupledOrNot = 0.0;

    // Get the maximum mobility and diffusivity values!
    testCarrier1.carrierProps.MaxMobility = 1.0;
    testCarrier1.carrierProps.MaxDiffusivity = 1.0;

    testCarrier2.carrierProps.MaxMobility = 1.0;
    testCarrier2.carrierProps.MaxDiffusivity = 1.0;

    testCarrier.carrierProps.MaxMobility = 1.0;
    testCarrier.carrierProps.MaxDiffusivity = 1.0;


    // basically sets lamba = 1 in poissons equation
    testPoisson.remakeProperties(testproblem, testgrid);

    // set the vandemondematrices and the datastructures needed
    // to project functions on the DG and CG basis.
    globalVandeMondeDG = testPoisson.VandeMondeMatrices.globalVandeMondeDG;

    globalVandeMondeMX = testPoisson.VandeMondeMatrices.globalVandeMondeMX;

    PTForward = testPoisson.Bijections.PTForward;

    sparseWeights = testPoisson.weightsSparse;

    MassU = testCarrier.carrierProps.MassU;
    A00 = testPoisson.PoissonProps.A00;

    BIJFlagDG = DG;
    BIJFlagMX = Mixed;

    return 0;
}
int
ddpTest_type::
initializeReactiveFluxTest(const ddpGrid_type & left_grid,
                           const ddpGrid_type & right_grid,
                           const ddpProblemInfo_type & inputproblem,
                           const ddpCarrierConstants_type & testConstants)
{
    //This just needs to be set before initializing the carriers
    // set the domain, difusivity, mobility, etc. for the testing
    // but does not set anything specifically related tests.
    // use the setTestCarrier or setTestPoisson functions for those
    // purposes.

    testgrid = left_grid;
    testgrid2 = right_grid;

    testproblem = inputproblem;

    testCarrier.initialize(testgrid, testproblem, testConstants, "electrons");
    testCarrier1.initialize(testgrid2, testproblem, testConstants, "reductants");
    testCarrier2.initialize(testgrid2, testproblem, testConstants, "oxidants");


    // initialize poisson
    // doesnst matter since no electic field
    testPoisson.initialize(testgrid, testproblem, testConstants);

    // Set the bijections
    DGForward = testPoisson.Bijections.DGForward;
    MXForward = testPoisson.Bijections.MXForward;

    // Convert mobility and diffusivity to identity matrix!
    testCarrier1.carrierProps.CurrentScale = 1.0;
    testCarrier2.carrierProps.CurrentScale = 1.0;
    testCarrier.carrierProps.CurrentScale = 1.0;

    //Set it so that coupling is off.
    testPoisson.PoissonProps.coupledOrNot = 0.0;

    // Get the maximum mobility and diffusivity values!
    testCarrier1.carrierProps.MaxMobility = 1.0;
    testCarrier1.carrierProps.MaxDiffusivity = 1.0;

    testCarrier2.carrierProps.MaxMobility = 1.0;
    testCarrier2.carrierProps.MaxDiffusivity = 1.0;

    testCarrier.carrierProps.MaxMobility = 1.0;
    testCarrier.carrierProps.MaxDiffusivity = 1.0;


    // basically sets lamba = 1 in poissons equation
    testPoisson.remakeProperties(testproblem, testgrid);

    // set the vandemondematrices and the datastructures needed
    // to project functions on the DG and CG basis.
    globalVandeMondeDG = testPoisson.VandeMondeMatrices.globalVandeMondeDG;

    globalVandeMondeMX = testPoisson.VandeMondeMatrices.globalVandeMondeMX;

    PTForward = testPoisson.Bijections.PTForward;

    sparseWeights = testPoisson.weightsSparse;

    MassU = testCarrier1.carrierProps.MassU;
    A00 = testPoisson.PoissonProps.A00;

    BIJFlagDG = DG;
    BIJFlagMX = Mixed;

    return 0;
}


// for poisson
int
ddpTest_type::
ReinitializeTestWithInterface(const ddpGrid_type & left_grid,
                              const ddpGrid_type & right_grid,
                              const ddpProblemInfo_type & inputproblem,
                              const ddpCarrierConstants_type & testConstants)
{

    /*	cout << "A00 = \n" << A00 << endl;
    	cout << "MassU = \n" << MassU << endl;
    	cout << "globalVandeMondeDG = \n " << globalVandeMondeDG << endl;
    	cout << "globalVandeMondeMX= \n " << globalVandeMondeMX << endl;
    */

    testgrid = left_grid;
    testgrid2 = right_grid;
    testproblem = inputproblem;

    // set these to 1.0 so lambda will be 1.0
    testproblem.vacuumPermittivity = 1.0;
    testproblem.thermalVoltage = 1.0;
    testproblem.characteristicLength = 1.0;
    testproblem.electronCharge = 1.0;

    // set the permitivitie so we can get A00 over [-1,1] Without the
    // change in permativity, so we can figure out the DOF of true solution
    testproblem.semiCondRelativePerm = 1.0;
    testproblem.electrolyteRelativePerm = 1.0;
    testPoisson.initializeWithInterface(testgrid, testgrid2, testproblem, testConstants);
    A00 = testPoisson.PoissonProps.A00; // If we didnt have an extra A00
    // when project true solution we would
    // invert true solution dof
    // with a discontinuity when it should
    // have a continuity


    // Now change the epsilson and build everything else
    testproblem.semiCondRelativePerm = M_PI;
    testproblem.electrolyteRelativePerm = 2.0 * M_PI;

//	cout << "REDOING" << endl;

    testPoisson.initializeWithInterface(testgrid, testgrid2, testproblem, testConstants);


    //Set it so that coupling is on.
    testPoisson.PoissonProps.coupledOrNot = 1.0;

    // reset the grid to go from [-1, 1]
//	ddpPrintGrid(testgrid);
    testgrid = testPoisson.Poisson_grid;
//	ddpPrintGrid(testgrid);
//	assert(false);

    // reset A00

    // Set the bijections
    DGForward = testPoisson.Bijections.DGForward;
    MXForward = testPoisson.Bijections.MXForward;

//	cout << "NumELements = " << testgrid.NumElementsWithGhost << endl;

    // set the vandemondematrices and the datastructures needed
    // to project functions on the DG and CG basis.
    globalVandeMondeDG = testPoisson.VandeMondeMatrices.globalVandeMondeDG;
    globalVandeMondeMX = testPoisson.VandeMondeMatrices.globalVandeMondeMX;

    PTForward = testPoisson.Bijections.PTForward;

    sparseWeights = testPoisson.weightsSparse;

    int N = sparseWeights.size();
    // make the pt values for epsilon
    EpsilonPTValues = ddpDenseVector_type(N);
    EpsilonInvPTValues = ddpDenseVector_type(N);
    for(int i = 0; i < N; i++)
    {
        if( i < N/2 )
        {
            EpsilonPTValues(i) = M_PI;
            EpsilonInvPTValues(i) = 1.0/M_PI;
        }
        else
        {
            EpsilonInvPTValues(i) = 0.5/M_PI;
            EpsilonPTValues(i) = 2.0 * M_PI;
        }
    }

    //cout << MassU << endl;

    // REDO MASS U so it is for over [-1, +1]
    ddpBijFlag_type FlagI = DG;
    ddpMakeGeneralMatrix(testproblem, testgrid, sparseWeights,
                         FlagI, DGForward, globalVandeMondeDG,
                         FlagI, DGForward, globalVandeMondeDG,
                         MassU);

    /*
    	cout << "A00 = \n " << A00 << endl;
    	cout << "MassU = \n" << MassU << endl;
    	cout << "DFDof = " << DGForward.size() << endl;

    	cout << "globalVandeMondeDG = \n " << globalVandeMondeDG << endl;
    	cout << "globalVandeMondeMX= \n " << globalVandeMondeMX << endl;

    	assert(false);
    */

    BIJFlagDG = DG;
    BIJFlagMX = Mixed;

    return 0;
}

int ddpTest_type::
setImplicitTestCarrier()
{
    setSteadyStateTestCarrier();
}


int
ddpTest_type::
setSteadyStateTestCarrier()
{
    // This sets all the carrier and poisson datastructures, boundary
    // conditions and initial conditions for testing the DG solver for
    // the drift diffusion equations for carrier transport
    // after the necessary parameters have been specified


    // set the carriers left BC type
    if(CarrierLeftBCType == Dirichlet)
    {
        testCarrier.setLeftTestDirichletBC(CarrierLeftBCFunction);
    }
    else if(CarrierLeftBCType == Robin)
    {
        testCarrier.setLeftTestRobincBC(CarrierLeftBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // set the carriers right BC type
    if(CarrierRightBCType == Dirichlet)
    {
        testCarrier.setRightTestDirichletBC(CarrierRightBCFunction);
    }
    else if(CarrierRightBCType == Robin)
    {
        testCarrier.setRightTestRobinBC(CarrierRightBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // remake (ALL THE FLUXES MATRICES)
    // NOTE: This must be done so that the diffusive fluxes can be remade to pick
    // up Robin boundary conditions if necessary.
    testCarrier.assembleLDGMatrices(testgrid,testproblem);
    testCarrier.setSolver(testgrid, testproblem);


    // set the initial conditions
    testCarrier.setInitialConditions(testgrid, CarrierICFunction);

    // set Phi(0), Phi(pi)
    testPoisson.setDirichletBoundaryConditions(PotentialLeftBCFunction,
            PotentialRightBCFunction);

    // Now set C(x) = 0
    DopingProfile.setDopingProfile(testgrid, testCarrier,DopingProfileFunction);

    return 0;
}

int
ddpTest_type::
setTestDoping()
{
    // This functions sets the dopingprofile object to the specified
    // parameters of the test.

    // Now set C(x)
    DopingProfile.setDopingProfile(testgrid, testCarrier1, testCarrier2,
                                   DopingProfileFunction1, DopingProfileFunction2);

    // Get DOFs for C(x)
    ddpProjectFunction(TrueDopingProfileFunction, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU, TrueDopingDof);

    return 0;
}

bool
ddpTest_type::
runTestDoping()
{
    // makes sure that the dopingprofile has been set correctly
    // by comparin the degrees of freedom for the set one and what
    // it should be.
    return 	Are_Same(DopingProfile.Dof, TrueDopingDof);
}




int
ddpTest_type::
setTestPoisson()
{
    // This sets all the carrier and poisson datastructures, boundary
    // conditions and initial conditions for testing Mixed CG solver for
    // the potential and electric field in poissons equation
    // after the necessary parameters have been specified


    // set the initial conditions for the charge carries
    testCarrier1.setInitialConditions(testgrid, Carrier1ICFunction);
    testCarrier2.setInitialConditions(testgrid, Carrier2ICFunction);

    // set Phi(0), Phi(pi)
    testPoisson.setDirichletBoundaryConditions(PotentialLeftBCFunction,
            PotentialRightBCFunction);

    // Now set C(x)
    DopingProfile.setDopingProfile(testgrid, testCarrier1, testCarrier2,
                                   DopingProfileFunction1, DopingProfileFunction2);

    //   Construct TRUE solution potential
    ddpProjectFunction(TrueSolutionPotential, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       truePoisson.PoissonState.potDof);

    // Construct the TRUE solution for E
    ddpProjectFunction(TrueSolutionElecField, testgrid, MXForward,
                       globalVandeMondeMX, PTForward, BIJFlagMX,
                       sparseWeights, A00,
                       truePoisson.PoissonState.elecDof);
    return 0;
}

bool
ddpTest_type::
runTestPoisson()
{
    // solve Poissons equation for the coupled and uncouple
    // problem.
    testPoisson.solveSystem(testCarrier1, testCarrier2,
                            DopingProfile, 0.0, testproblem);

//	L2_Error(testPoisson, truePoisson);

    // test to see if the calculated solutions and the analytic
    // solutions have the same DOFs
    return 	Dofs_Are_Same(testPoisson, truePoisson);

}

void
ddpTest_type::
runTestPoissonL2Error(ddpDenseVector_type & errors)
{
    // solve Poissons equation for the coupled and uncouple
    // problem.
    testPoisson.solveSystem(testCarrier1, testCarrier2,
                            DopingProfile, 0.0, testproblem);

    L2_Error(testPoisson, TrueSolutionPotential, TrueSolutionElecField, errors);
    // test to see if the calculated solutions and the analytic
    // solutions have the same DOFs

}

bool
ddpTest_type::
runTestPoissonWithInterface()
{
    // solve Poissons equation for the coupled and uncouple
    // problem.
    testPoisson.solveSystem(testCarrier1,
                            testCarrier2,
                            testCarrier1,
                            testCarrier2,
                            DopingProfile,
                            0.0, testproblem);

    // test to see if the calculated solutions and the analytic
    // solutions have the same DOF:s  ddpDenseVector_type
    ddpDenseVector_type
    testElecPTS = globalVandeMondeMX * testPoisson.PoissonState.elecDof,
    testPOTPTS = globalVandeMondeDG * testPoisson.PoissonState.potDof,
    trueElecPTS = globalVandeMondeMX * truePoisson.PoissonState.elecDof,
    truePOTPTS = globalVandeMondeDG * truePoisson.PoissonState.potDof;

    trueElecPTS = trueElecPTS.cwiseProduct(EpsilonInvPTValues);
    testElecPTS = testElecPTS.cwiseProduct(EpsilonInvPTValues);


    bool elecfieldSame = Are_Same(trueElecPTS, testElecPTS);
    bool potsSame = Are_Same(truePOTPTS, testPOTPTS);


    if( (true == potsSame) &&  (true == elecfieldSame))
        return true;
    else
        return false;

}



bool
ddpTest_type::
runImplicitTestCarrier()
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);

    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    double tCurrent = 0.0;
    double Tend = 1.0;
    double DeltaT = 0.001;

    testCarrier.AssembleLDGSystem(testproblem,
                                  testgrid,
                                  DeltaT);


    // Run time stepping from t= 0 to t = 1.0
    while(tCurrent < Tend)
    {
        testCarrier.UpdateExplicitDriftTerm(ElecFieldDof);
        testCarrier.AssembleRHS(tCurrent, DeltaT);
        testCarrier.BackwardEuler();
        tCurrent += DeltaT;
    }

    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.
    return Are_Same(testCarrier, trueCarrier);
}

bool
ddpTest_type::
runImplicitTestCarrierWithScaling()
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.
    double C = 5.0;
    double tau = 2.0;
    double L = M_PI;
    double mu = 1.0;

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);


    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    double tCurrent = 0.0;
    double Tend = 1.0;
    double DeltaT = 0.001;

    testCarrier.AssembleLDGSystem(testproblem,
                                  testgrid,
                                  DeltaT);



    // Run time stepping from t= 0 to t = 1.0
    while(tCurrent < Tend)
    {
        testCarrier.UpdateExplicitDriftTerm(ElecFieldDof);
        testCarrier.AssembleRHS(tCurrent, DeltaT);
        testCarrier.BackwardEuler();
        tCurrent += DeltaT;
    }

    testCarrier.carrierState.uDof *= C ;
    testCarrier.carrierState.qDof *= C*L/tau;

    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.
    return Are_Same(testCarrier, trueCarrier);
}

bool
ddpTest_type::
runSteadyStateShurTestCarrier()
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);

    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    /// solve the linear system
    testCarrier.ShurSolve(ElecFieldDof,
                          testproblem,
                          testgrid,
                          0.0);

    testCarrier.carrierState.qDof = trueCarrier.carrierState.qDof;

    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.
    return Are_Same(testCarrier, trueCarrier);
}

bool
ddpTest_type::
runSteadyStateTestCarrier()
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);

    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    /// solve the linear system
    testCarrier.Solve(ElecFieldDof,
                      testproblem,
                      testgrid,
                      0.0);

    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.
    return Are_Same(testCarrier, trueCarrier);
}
void
ddpTest_type::
runTestCarrierL2Error(ddpDenseVector_type & errors)
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);



    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    /// solve the linear system
    testCarrier.Solve(ElecFieldDof,
                      testproblem,
                      testgrid,
                      0.0);

    L2_Error(testCarrier, TrueSolutionU, TrueSolutionQ, errors);


}
void
ddpTest_type::
runImplicitTestCarrierL2Error(ddpDenseVector_type & errors)
{
    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier, DopingProfile, 0.0, testproblem);

    //   Construct TRUE solution
    ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
                       globalVandeMondeDG,	PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.uDof);

    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);



    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;
    double tCurrent = 0.0;
    double Tend = 1.0;
    double DeltaT = 0.001;

    testCarrier.AssembleLDGSystem(testproblem,
                                  testgrid,
                                  DeltaT);


    // Run time stepping from t= 0 to t = 1.0
    while(tCurrent < Tend)
    {
        testCarrier.UpdateExplicitDriftTerm(ElecFieldDof);
        testCarrier.AssembleRHS(tCurrent, DeltaT);
        testCarrier.BackwardEuler();
        tCurrent += DeltaT;
    }

    L2_Error(testCarrier, TrueSolutionU, TrueSolutionQ, errors);


}


int
ddpTest_type::
setTestReactiveFluxes()
{
    // This sets all the carrier and poisson datastructures, boundary
    // conditions and initial conditions for testing the DG solver for
    // the drift diffusion equations for carrier transport
    // after the necessary parameters have been specified


    // set the carriers left BC type
    if(CarrierLeftBCType == Dirichlet)
    {
        testCarrier.setLeftTestDirichletBC(CarrierLeftBCFunction);
    }
    else if(CarrierLeftBCType == Robin)
    {
        testCarrier.setLeftTestRobincBC(CarrierLeftBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // set the carriers left BC type
    if(Carrier1LeftBCType == Dirichlet)
    {
        testCarrier1.setLeftTestDirichletBC(Carrier1LeftBCFunction);
    }
    else if(Carrier1LeftBCType == Robin)
    {
        testCarrier1.setLeftTestRobincBC(Carrier1LeftBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // NOTE: CARRIER2 AND CARRIER WILL HAVE SAME BOUNDARY CONDITIONS

    if(Carrier2LeftBCType == Dirichlet)
    {
        testCarrier2.setLeftTestDirichletBC(Carrier2LeftBCFunction);
    }
    else if(Carrier2LeftBCType == Robin)
    {
        testCarrier2.setLeftTestRobincBC(Carrier2LeftBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // set the carriers right BC type
    if(CarrierRightBCType == Dirichlet)
    {
        testCarrier.setRightTestDirichletBC(CarrierRightBCFunction);
    }
    else if(CarrierRightBCType == Robin)
    {
        testCarrier.setRightTestRobinBC(CarrierRightBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    if(Carrier1RightBCType == Dirichlet)
    {
        testCarrier1.setRightTestDirichletBC(Carrier1RightBCFunction);
    }
    else if(Carrier1RightBCType == Robin)
    {
        testCarrier1.setRightTestRobinBC(Carrier1RightBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }


    if(Carrier2RightBCType == Dirichlet)
    {
        testCarrier2.setRightTestDirichletBC(Carrier2RightBCFunction);
    }
    else if(Carrier2RightBCType == Robin)
    {
        testCarrier2.setRightTestRobinBC(Carrier2RightBCFunction);
    }
    else
    {
        cout << "No BC of that kind" << endl;
        assert(false);
    }

    // remake (ALL THE FLUXES MATRICES)
    // NOTE: This must be done so that the diffusive fluxes can be remade to pick
    // up Robin boundary conditions if necessary.
    testCarrier.assembleLDGMatrices(testgrid,testproblem);
    testCarrier1.assembleLDGMatrices(testgrid,testproblem);
    testCarrier2.assembleLDGMatrices(testgrid,testproblem);



    // set Phi(0), Phi(pi)
    testPoisson.setDirichletBoundaryConditions(PotentialLeftBCFunction,
            PotentialRightBCFunction);

    // Now set C(x) = 0
    DopingProfile.setDopingProfile(testgrid, testCarrier1, DopingProfileFunction);

    return 0;
}

bool
ddpTest_type::
runImplicitTestReactiveFluxes()
{

    // NEEDS TO HAPPEN BEFORE SETTING INITIAL CONDITIONS
    testCarrier.setSolver(testgrid,testproblem);
    testCarrier1.setSolver(testgrid, testproblem);
    testCarrier2.setSolver(testgrid, testproblem);


    // set the initial conditions
    testCarrier.setInitialConditions(testgrid, CarrierICFunction);
    testCarrier1.setInitialConditions(testgrid, Carrier1ICFunction);
    testCarrier2.setInitialConditions(testgrid, Carrier2ICFunction);


    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier1, DopingProfile, 0.0, testproblem);

    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    testPoisson.getElecFieldVals(EPTS);

    // prepare for timestepping

    double tCurrent = 0.0;
    double Tend = 10.0;
    double ket = 100.0;
    double DeltaT = 0.0001;

    ddpDenseVector_type timeStamps(100);

    double timeStampDeltaT = Tend/100;

    for(int i = 0; i < 100; i++)
        timeStamps(i) = (i+1)*timeStampDeltaT;

    double u1_RightDensity;
    double u2_LeftDensity;
    double u1_RightFlux;

    ddpPrintStateTest(testproblem,
                      testCarrier,
                      testCarrier1,
                      testCarrier2,
                      testPoisson,
                      0);


    // electron
    testCarrier.carrierProps.Material = Semiconductor;

    // reductant
    testCarrier1.carrierProps.Material = Electrolyte;

    // oxidant
    testCarrier2.carrierProps.Material = Electrolyte;


    // for loop over the time stamps
    for(int stamp = 0; stamp < 100; stamp++)
    {
        while(tCurrent < timeStamps(stamp))
        {
            // electron density
            u1_RightDensity = testCarrier.getInterfaceDensity();

            // oxidant density
            u2_LeftDensity = testCarrier2.getInterfaceDensity();

            // electrons implicit data
            testCarrier.AssembleLDGSystem(testproblem,
                                          testgrid,
                                          ket * u2_LeftDensity,
                                          DeltaT);


            // electons explicit data
            testCarrier.setRight_RobinValue(-ket * 0.5 * u2_LeftDensity);

            testCarrier.UpdateExplicitDriftTerm(ElecFieldDof);
            testCarrier.AssembleRHS(tCurrent, DeltaT);
            testCarrier.BackwardEuler();

            // update electron density at interface
            u1_RightDensity = testCarrier.getInterfaceDensity();

            // oxidants implicit data
            testCarrier2.AssembleLDGSystem(testproblem,
                                           testgrid,
                                           ket * (u1_RightDensity - 0.5),
                                           DeltaT);




            // oxidants explicit data
            testCarrier2.setLeft_RobinValue(0.0); // * u2_LeftDensity);

            testCarrier2.UpdateExplicitDriftTerm(ElecFieldDof);
            testCarrier2.AssembleRHS(tCurrent, DeltaT);
            testCarrier2.BackwardEuler();

            tCurrent += DeltaT;
        }

        ddpPrintStateTest(testproblem,
                          testCarrier,
                          testCarrier1,
                          testCarrier2,
                          testPoisson,
                          (stamp+1));
    }
    /*
    	//   Construct TRUE solution
    	ddpProjectFunction(TrueSolutionU, testgrid, DGForward,
    					 	   	 		globalVandeMondeDG,	PTForward, BIJFlagDG,
    	    		 				  sparseWeights, MassU,
    	    		  				trueCarrier.carrierState.uDof);

    	// Construct the TRUE solution for q
     	ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
    	    		   					globalVandeMondeDG, PTForward, BIJFlagDG,
    	    		   					sparseWeights, MassU,
    			   							trueCarrier.carrierState.qDof);
    */
    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.



}
bool
ddpTest_type::
runImplicitExplicitTestReactiveFluxes()
{

    testCarrier.setSolver(testgrid,testproblem);
    testCarrier1.setSolver(testgrid2, testproblem);
    testCarrier2.setSolver(testgrid2, testproblem);

    // set the initial conditions
    testCarrier.setInitialConditions(testgrid, CarrierICFunction);
    testCarrier1.setInitialConditions(testgrid2, Carrier1ICFunction);
    testCarrier2.setInitialConditions(testgrid2, Carrier2ICFunction);

    // solve the Drift diffusion equations with a given electric field
    // up until a time T = 1.0 and then compare the calcuted solutions
    // and analytic solutions at every point in the mesh.

    // Solve for Electric field at time = 0.0 (wont change with time)
    testPoisson.solveSystem(testCarrier1, DopingProfile, 0.0, testproblem);

    // get the electric field
    ElecFieldDof = testPoisson.PoissonState.elecDof;

    testPoisson.getElecFieldVals(EPTS);

    // prepare for timestepping

    double tCurrent = 0.0;
    double Tend = 15.0;
    double ket = 100.0;
    double DeltaT = 0.0001;

    double u1_RightDensity;
    double u2_LeftDensity;
    double u1_RightFlux;

    ddpDenseVector_type timeStamps(100);

    double timeStampDeltaT = Tend/100;

    for(int i = 0; i < 100; i++)
        timeStamps(i) = (i+1)*timeStampDeltaT;

    ddpPrintStateTest(testproblem,
                      testCarrier,
                      testCarrier1,
                      testCarrier2,
                      testPoisson,
                      0);

    testCarrier.AssembleLDGSystem(testproblem,
                                  testgrid,
                                  DeltaT);

    testCarrier1.AssembleLDGSystem(testproblem,
                                   testgrid2,
                                   DeltaT);

    testCarrier2.AssembleLDGSystem(testproblem,
                                   testgrid2,
                                   DeltaT);

    for(int stamp = 0; stamp < 100; stamp++)
    {
        while(tCurrent < timeStamps(stamp))
        {
            // electron density
            u1_RightDensity = testCarrier.getInterfaceDensity();

            // oxidant density
            u2_LeftDensity = testCarrier2.getInterfaceDensity();

            // calculate boundary conditions: kf * \rho_{o}
            u1_RightFlux = ket*(u1_RightDensity - 0.5)*u2_LeftDensity;

            // set the boundary conditions
            testCarrier.setRight_RobinValue(u1_RightFlux);
            testCarrier.UpdateExplicitDriftTerm(ElecFieldDof);
            testCarrier.AssembleRHS(tCurrent, DeltaT);
            testCarrier.BackwardEuler();


            testCarrier1.setLeft_RobinValue(u1_RightFlux);
            testCarrier2.setLeft_RobinValue(-1.0*u1_RightFlux);

            testCarrier1.AssembleRHS(tCurrent, DeltaT);
            testCarrier1.UpdateExplicitDriftTerm(ElecFieldDof);
            testCarrier1.BackwardEuler();

            testCarrier2.UpdateExplicitDriftTerm(ElecFieldDof);
            testCarrier2.AssembleRHS(tCurrent, DeltaT);
            testCarrier2.BackwardEuler();

            tCurrent += DeltaT;
        }


        ddpPrintStateTest(testproblem,
                          testCarrier,
                          testCarrier1,
                          testCarrier2,
                          testPoisson,
                          stamp+1);
    }


    // Construct the TRUE solution for q
    ddpProjectFunction(TrueSolutionQ, testgrid, DGForward,
                       globalVandeMondeDG, PTForward, BIJFlagDG,
                       sparseWeights, MassU,
                       trueCarrier.carrierState.qDof);

    // return the result of testing whether the two solutions
    // give the same point value result.  That is whether our
    // numerical solution has converged in the L_{\infinity} norm.

    bool result1 = Are_Same(trueCarrier.carrierState.qDof,
                            testCarrier.carrierState.qDof);

    bool result2 = Are_Same(trueCarrier.carrierState.qDof,
                            testCarrier2.carrierState.qDof);

    if(result1 && result2)
        return true;
    else
        return false;
}

int
ddpTest_type::
outputTestCarrier(int testNumber)
{

    // print the tests carriers, true carriers, and electric fields
    // solution and derivative of the solutions point values.
    ofstream prt;
    std::string prefix = "testCarrier";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    ss << prefix << testNumber << extension << '\0';
    filename = ss.str();
    char str[filename.size()];

    for(unsigned int i = 0; i < filename.size(); i++)
    {
        str[i] = filename[i];
    }

    prt.open(str);

    ddpDenseVector_type
    UPTS_test = globalVandeMondeDG * testCarrier.carrierState.uDof,
    UPTS_true = globalVandeMondeDG * trueCarrier.carrierState.uDof;

    ddpDenseVector_type
    QPTS_test = globalVandeMondeDG * testCarrier.carrierState.qDof,
    QPTS_true = globalVandeMondeDG * trueCarrier.carrierState.qDof;

    ddpDenseVector_type
    ElecPTS = globalVandeMondeMX * testPoisson.PoissonState.elecDof,
    POTPTS = globalVandeMondeDG * testPoisson.PoissonState.potDof;

    for(int i = 0; i < UPTS_test.size(); i++)
    {
        prt << testCarrier.PTSDense(i) << " "
            << UPTS_test(i) << " "
            << QPTS_test(i) << " "
            << UPTS_true(i) << " "
            << QPTS_true(i) << " "
            << ElecPTS(i) << " "
            << POTPTS(i) << endl;
    }


    prt.close();
    return 0;
}

int
ddpTest_type::
outputTestPoisson(int testNumber)
{

    // outputs the results of the test numbered testNumber
    // for poisson, Will output the potential and electric field
    // of the numericallly computed solution and analytic solution
    ofstream prt;
    std::string prefix = "testPoisson";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    ss << prefix << testNumber << extension << '\0';
    filename = ss.str();
    char str[filename.size()];

    for(unsigned int i = 0; i < filename.size(); i++)
        str[i] = filename[i];

    prt.open(str);

    ddpDenseVector_type
    testElecPTS = globalVandeMondeMX * testPoisson.PoissonState.elecDof,
    testPOTPTS = globalVandeMondeDG * testPoisson.PoissonState.potDof,
    trueElecPTS = globalVandeMondeMX * truePoisson.PoissonState.elecDof,
    truePOTPTS = globalVandeMondeDG * truePoisson.PoissonState.potDof;



    for(int i = 0; i < testElecPTS.size(); i++)
    {
        prt << testPoisson.PTSDense(i) << " "
            << testElecPTS(i) << " "
            << testPOTPTS(i) << " "
            << trueElecPTS(i) << " "
            << truePOTPTS(i) <<  endl;
    }


    prt.close();
    return 0;

}

int
ddpTest_type::
outputTestPoissonWithInterface(int testNumber)
{

    // outputs the results of the test numbered testNumber
    // for poisson, Will output the potential and electric field
    // of the numericallly computed solution and analytic solution
    ofstream prt;
    std::string prefix = "testPoisson";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    ss << prefix << testNumber << extension << '\0';
    filename = ss.str();
    char str[filename.size()];

    for(unsigned int i = 0; i < filename.size(); i++)
        str[i] = filename[i];

    prt.open(str);

    ddpDenseVector_type
    testElecPTS = globalVandeMondeMX * testPoisson.PoissonState.elecDof,
    testPOTPTS = globalVandeMondeDG * testPoisson.PoissonState.potDof,
    trueElecPTS = globalVandeMondeMX * truePoisson.PoissonState.elecDof,
    truePOTPTS = globalVandeMondeDG * truePoisson.PoissonState.potDof;

    trueElecPTS = trueElecPTS.cwiseProduct(EpsilonInvPTValues);
    testElecPTS = testElecPTS.cwiseProduct(EpsilonInvPTValues);

    for(int i = 0; i < testElecPTS.size(); i++)
    {
        prt << testPoisson.PTSDense(i) << " "
            << testElecPTS(i) << " "
            << testPOTPTS(i) << " "
            << trueElecPTS(i) << " "
            << truePOTPTS(i) <<  endl;
    }

    prt.close();
    return 0;

}




