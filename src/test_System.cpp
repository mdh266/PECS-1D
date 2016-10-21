#include "../include/test.hpp"


int ddpTestingMakeInputFile(void);

int main()
{

    // Initial set up stuff
    ddpTestingMakeInputFile();

    ddpDomain_type testdomain1;
    ddpGrid_type testgrid1;
    ddpProblemInfo_type testproblem;
    ddpCarrierConstants_type testConstants;

    readInput(testdomain1, testgrid1, testproblem, testConstants, "ddp-test.ini");

    //Changing end points so that I can come up with manufactored solutions
    testdomain1.LeftEndPoint = -M_PI;
    testdomain1.RightEndPoint = 0;
    ddpMakeUniformGrid(testproblem, testdomain1, testgrid1);

    ddpDomain_type testdomain2;
    ddpGrid_type testgrid2;
    readInput(testdomain2, testgrid2, testproblem, testConstants, "ddp-test.ini");

    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = M_PI;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);


//	testgrid.NumBoundaryElements = 50;
//	testgrid.BoundaryLayerWidth = 0.5;
//	ddpMakeRefinedBoundaryGrid(testproblem, testdomain, testgrid);
    testproblem.BoltzmannConst = 1.0;
    testproblem.temperature = 1.0;
    testproblem.electronCharge = 1.0;
    testproblem.vacuumPermittivity = 1.0;


    // Test For Electrons
    ddpTestChargeCarrier_type carrierType = NegativeCharge;
    ddpTest_type testElectrons;
    testElectrons.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    // Test For Holes
    carrierType = PositiveCharge;
    ddpTest_type testHoles;
    testHoles.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    // Test For Poisson's Equation
    carrierType = PositiveCharge; // doesnt matter
    ddpTest_type testPoisson;
    testPoisson.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    // Test the Reactive fluxes
    ddpTest_type testReactiveFluxes;
    testReactiveFluxes.initializeReactiveFluxTest(testgrid1, testgrid2, testproblem, testConstants);



    bool testPassed;
    int numTestsPassed = 0;
    int numTestsTotal = 0;



///////////////////////////////////////////////////////////////////////////////
// Doping Profile
//////////////////////////////////////////////////////////////////////////////

// Test to make sure Doping profile is being setting up correctly.
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << "\t \t \t \t DOPING  \n" << endl;
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << endl << endl;


    cout << "\n 1.) Testing: C(x) = N_D(x)" << endl << endl;

    // set C = N_D
    testPoisson.DopingProfileFunction1 = &Doping_N_D;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    testPoisson.TrueDopingProfileFunction = &Doping_N_D;

    testPoisson.setTestDoping();

    testPassed = testPoisson.runTestDoping();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    cout << "\n 2.) Testing: C(x) = -N_A(x)" << endl << endl;

    // set C = -N_A
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &Doping_N_A;

    testPoisson.TrueDopingProfileFunction = &True_N_A;

    testPoisson.setTestDoping();

    testPassed = testPoisson.runTestDoping();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    cout << "\n 3.) Testing: C(x) = N_D(x) - N_A(x)" << endl << endl;

    // set C = N_D - N_A
    testPoisson.DopingProfileFunction1 = &Doping_N_D;
    testPoisson.DopingProfileFunction2 = &Doping_N_A;

    testPoisson.TrueDopingProfileFunction = &Doping_PN;

    testPoisson.setTestDoping();

    testPassed = testPoisson.runTestDoping();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

////////////////////////////////////////////////////////////////////////////
// Poisson's Equation
////////////////////////////////////////////////////////////////////////////

// Tests to make sure solving Poissons' Equation correctly.

    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << "\t \t \t \t MIXED METHOD POISSON \n" << endl;
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << endl << endl;

// turn on coupling (include n and p on RHS of Poisson's eq)
    testPoisson.testPoisson.PoissonProps.coupledOrNot = 1.0;

    cout << "1.)  No Doping\n" << endl;
    // set n = p = 0
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = 0
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 1
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &OneFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &OneOverPiLineFunction;
    testPoisson.TrueSolutionElecField = &MinusOneOverPiFunction;


    testPoisson.setTestPoisson();
    testPassed = testPoisson.runTestPoisson();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

//	assert(false);


    // output results
//  testPoisson.outputTestPoisson(1);



    cout << "2.)  Doping With Sin(x)\n" << endl;

    // set n = p = 0
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = sin(x)
    testPoisson.DopingProfileFunction1 = &SineFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &ZeroFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &SineFunction;
    testPoisson.TrueSolutionElecField = &MinusCosineFunction;


    testPoisson.setTestPoisson();

    testPassed = testPoisson.runTestPoisson();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testPoisson.outputTestPoisson(2);


    cout << "3.)  No Doping: n = sin(x)\n" << endl;

    // set n = sin, p = 0
    testPoisson.Carrier1ICFunction = &SineFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = 0
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &ZeroFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &MinusSineFunction;
    testPoisson.TrueSolutionElecField = &CosineFunction;

    testPoisson.setTestPoisson();

    testPassed = testPoisson.runTestPoisson();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
// testPoisson.outputTestPoisson(3);

    cout << "4.)  No Doping: p = sin(x)\n" << endl;

    // set n = 0 p = sin
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &SineFunction;

    // set C = 0
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &ZeroFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &SineFunction;
    testPoisson.TrueSolutionElecField = &MinusCosineFunction;

    testPoisson.setTestPoisson();

    testPassed = testPoisson.runTestPoisson();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testPoisson.outputTestPoisson(4);


    cout << "5.)  No Doping: n, p = sin(x)\n" << endl;

    // set n = sin p = sin
    testPoisson.Carrier1ICFunction = &SineFunction;
    testPoisson.Carrier2ICFunction = &SineFunction;

    // set C = 0
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 1
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &OneFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &OneOverPiLineFunction;
    testPoisson.TrueSolutionElecField = &MinusOneOverPiFunction;

    testPoisson.setTestPoisson();

    testPassed = testPoisson.runTestPoisson();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testPoisson.outputTestPoisson(5);


///////////////////////////////////////////////////////////////////////////////
// Test Poisson With an interface
///////////////////////////////////////////////////////////////////////////////

    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << "\t \t \t  MIXED METHOD POISSON W/ INTERFACE \n" << endl;
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << endl << endl;

    //Changing end points so that I can come up with manufactored solutions
    testdomain1.LeftEndPoint = -1.0;
    testdomain1.RightEndPoint = 0.0;
    ddpMakeUniformGrid(testproblem, testdomain1, testgrid1);


    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = 1.0;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);


    cout << "6.)  Zero RHS: \n\tN_D, NA = 0,\n"
         << "\trho_n, rho_p, rho_r, rho_o = 0 \n" << endl;

    testPoisson.ReinitializeTestWithInterface(testgrid1,
            testgrid2,
            testproblem,
            testConstants);

    // rho_n, rho_r = 0, rho_p, rho_o = 0
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = 0
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &OneFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &PotentialInterfaceZeroRHS;
    testPoisson.TrueSolutionElecField = &DisplElecFieldInterfaceZeroRHS;


    testPoisson.setTestPoisson();
    testPassed = testPoisson.runTestPoissonWithInterface();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }


    // output results
// 	testPoisson.outputTestPoissonWithInterface(6);

    cout << "7.) Nonzero RHS:\n \tN_D = 1, N_A = 0,\n"
         << "\trho_n, rho_r = 0, rho_p, rho_o = 0 \n" << endl;

    // rho_n = \rho_r = 0, \rho_o = \rho_p = 0
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = N_{D} - N_{A} = 1 - 0 =  1
    testPoisson.DopingProfileFunction1 = &OneFunction;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 1.0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &OneFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &PotentialInterfaceNonZeroRHS1;
    testPoisson.TrueSolutionElecField = &DisplElecFieldInterfaceNonZeroRHS1;


    testPoisson.setTestPoisson();
    testPassed = testPoisson.runTestPoissonWithInterface();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
// 	testPoisson.outputTestPoissonWithInterface(7);

    cout << "8.) Nonzero RHS:\n \tN_D = 1, N_A = 0,\n"
         << "\trho_n, rho_r = 0, rho_p, rho_o = 1 \n" << endl;

    // rho_n = \rho_r = 0, \rho_o = \rho_p = 1
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &OneFunction;

    // set C = N_{D} - N_{A} = 0 -1 =
    testPoisson.DopingProfileFunction1 = &ZeroFunction;
    testPoisson.DopingProfileFunction2 = &OneFunction;

    // set Phi(0) = 0, Phi(pi) = 1.0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &OneFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &PotentialInterfaceNonZeroRHS2;
    testPoisson.TrueSolutionElecField = &DisplElecFieldInterfaceNonZeroRHS2;


    testPoisson.setTestPoisson();
    testPassed = testPoisson.runTestPoissonWithInterface();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
// 	testPoisson.outputTestPoissonWithInterface(8);

//////////////////////////////////////////////////////////////////////////////
// LDG STEADY STATE SOLVERS
/////////////////////////////////////////////////////////////////////////////

    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << "\t \t \t LDG STEADY STATE SOLVERS\n" << endl;
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << endl << endl;

    cout << "\t ************************* \n \n";
    cout << "\t   MONOLITHIC SOLVERS \n \n";
    cout << "\t ************************* \n \n" << endl;

    cout << "1.) Poisson w/ Dirichlet \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &OneFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &OneOverPiLineFunction;
    testHoles.TrueSolutionQ = &MinusOneOverPiFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(1);


    cout << "2.) Poisson w/ Neumann \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Robin;
    testHoles.CarrierLeftBCFunction = &OneFunction;
    testHoles.CarrierRightBCFunction = &PiFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &Diffusion_SteadyState_Robin_U;
    testHoles.TrueSolutionQ = &PiFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(2);


    cout << "3.) Drift-Diffusion w/ Dirichlet \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &OneFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &OneFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_SteadyState_Dir_U;
    testHoles.TrueSolutionQ = &HDD_SteadyState_Dir_Q;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(3);

    cout << "4.) Drift-Diffusion w/ Robin \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Robin;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &OneFunction;
    testHoles.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &OneFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_SteadyState_Robin_U;
    testHoles.TrueSolutionQ = &OneFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(4);

//////////////////////////////////////////////////////////////////////////////

    cout << "\t ************************* \n \n";
    cout << "\t  SHUR COMPLEMENT SOLVER\n \n";
    cout << "\t ************************* \n \n" << endl;

    cout << "5.) Poisson w/ Dirichlet \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &OneFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &OneOverPiLineFunction;
    testHoles.TrueSolutionQ = &MinusOneOverPiFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateShurTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(5);


    cout << "6.)Poisson w/ Neumann \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Robin;
    testHoles.CarrierLeftBCFunction = &OneFunction;
    testHoles.CarrierRightBCFunction = &PiFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &Diffusion_SteadyState_Robin_U;
    testHoles.TrueSolutionQ = &PiFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateShurTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(6);


    cout << "7.) Drift-Diffusion w/ Dirichlet \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &OneFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &OneFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_SteadyState_Dir_U;
    testHoles.TrueSolutionQ = &HDD_SteadyState_Dir_Q;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateShurTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(7);


    cout << "8.) Drift-Diffusion w/ Robin \n\n";

    testHoles.CarrierICFunction = &ZeroFunction;

    testHoles.CarrierLeftBCType = Robin;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &OneFunction;
    testHoles.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &OneFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_SteadyState_Robin_U;
    testHoles.TrueSolutionQ = &OneFunction;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testPassed = testHoles.runSteadyStateShurTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(8);



//////////////////////////////////////////////////////////////////////////////
// LDG TRANSIENT SOLVERS
/////////////////////////////////////////////////////////////////////////////

    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << "\t \t \t LDG TRANSIENT SOLVERS\n" << endl;
    cout << "\t \t --------------------------------------------------  \n" << endl;
    cout << endl << endl;

///////////////////////////////////////////////////////////////////////////////
// DIFFUSION EQUATION
///////////////////////////////////////////////////////////////////////////////
    cout << "\t ************************* \n \n";
    cout << "\t Diffusion w/ Implicit Time Stepping\n" << endl;
    cout << "\t ************************* \n \n";
// u_t = u_xx

// E = 0 in drift diffusion equation

    cout << "14.) Dirichlet Homogeneous\n\n";

// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &SineFunction;

    // set n(0,t) = 0, n(pi,t) = 0
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Dirichlet;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction;
    testElectrons.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionHomoU;
    testElectrons.TrueSolutionQ = &DiffusionSolutionHomoQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testElectrons.outputTestCarrier(14);


    cout << "15.) Dirichlet Non-Homogenous.\n\n";

    // set n(x,0) = cos(x)
    testElectrons.CarrierICFunction = &CosineFunction;

    // set n(0,t) = exp(-t), n(pi,t) = -exp(-t)
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Dirichlet;
    testElectrons.CarrierLeftBCFunction = &DiffusionNonHomoLeftBC;
    testElectrons.CarrierRightBCFunction = &DiffusionNonHomoRightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNonHomoU;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNonHomoQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }
    // output results
//  testElectrons.outputTestCarrier(15);


// E = 0 in drift diffusion equation

    cout << "16.)  Homogeneous Neumann\n\n";

// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &Cosine2Function;

    // set n_x(0,t) = 0, n_x(pi,t) = 0
    testElectrons.CarrierLeftBCType = Robin;
    testElectrons.CarrierRightBCType = Robin;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction;
    testElectrons.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNeumannHomoU;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNeumannHomoQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testElectrons.outputTestCarrier(16);

    cout << "17.) Non-Homogeneous Neumann\n\n";

// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &Sine2Function;

    // set n_x(0,t) = 2 exp(-4t) , n_x(pi,t) = 2 exp(-4t)
    testElectrons.CarrierLeftBCType = Robin;
    testElectrons.CarrierRightBCType = Robin;
    testElectrons.CarrierLeftBCFunction = &DiffusionNeumannNonHomoBC;
    testElectrons.CarrierRightBCFunction = &DiffusionNeumannNonHomoBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNeumannNonHomoU;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNeumannNonHomoQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;

    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testElectrons.outputTestCarrier(17);

// E = 0 in drift diffusion equation

    cout << "18.)  Left Homogenous Dirichlet, Right Non-Homogenous Neumann\n\n";

// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &Sine2Function;

    // set n(0,t) = 0, n_x(pi,t) = 2*exp(-4t)
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Robin;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction;
    testElectrons.CarrierRightBCFunction = &DiffusionNeumannNonHomoBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNeumannNonHomoU;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNeumannNonHomoQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testElectrons.outputTestCarrier(18);

    cout << "19.) Scaling, Non-Homogenous Dirichlet.\n\n";
    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = 1.0; //.0;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    testproblem.characteristicLength = M_PI;
    testproblem.characteristicTime = 2.0;
    testproblem.characteristicDensity = 5.0;
    testConstants.electron_Mobility = 1.0;
    testproblem.thermalVoltage = 1.0;

    carrierType = NegativeCharge;
    testElectrons.ReinitializeTestWithScaling(carrierType, testgrid2,
            testproblem, testConstants);

    // set n(x,0) = cos(x)
    testElectrons.CarrierICFunction = &CosinePiFunction;

    // set n(0,t) = exp(-t), n(pi,t) = -exp(-t)
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Dirichlet;
    testElectrons.CarrierLeftBCFunction = &DiffusionNonHomoLeftBC_Scaling;
    testElectrons.CarrierRightBCFunction = &DiffusionNonHomoRightBC_Scaling;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNonHomoU_Scaling;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNonHomoQ_Scaling;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrierWithScaling();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }
    // output results
//  testElectrons.outputTestCarrier(19);

    cout << "20.) Scaling, Homogenous Neumann.\n\n";
//	testdomain2.LeftEndPoint = 0.0;
//  testdomain2.RightEndPoint = M_PI; //.0;
//  ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    testproblem.characteristicLength = M_PI;
    testproblem.characteristicTime = 2.0;
    testproblem.characteristicDensity = 5.0;
    testConstants.electron_Mobility = 1.0;
    testproblem.thermalVoltage = 1.0;

    carrierType = NegativeCharge;
    testElectrons.ReinitializeTestWithScaling(carrierType, testgrid2,
            testproblem, testConstants);
// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &Cosine2PiFunction;

    // set n_x(0,t) = 2 exp(-4t) , n_x(pi,t) = 2 exp(-4t)
    testElectrons.CarrierLeftBCType = Robin;
    testElectrons.CarrierRightBCType = Robin;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction; //DiffusionNeumannNonHomoBC_Scaling;
    testElectrons.CarrierRightBCFunction = &ZeroFunction; //DiffusionNeumannNonHomoBC_Scaling;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionSolutionNeumannHomoU_Scaling;
    testElectrons.TrueSolutionQ = &DiffusionSolutionNeumannHomoQ_Scaling;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrierWithScaling();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;

    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }


    // output results
//  testElectrons.outputTestCarrier(20);


    cout << "20.) Scaling, Non-Homogenous Neumann.\n\n";
    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = 1.0; //.0;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    testproblem.characteristicLength = M_PI;
    testproblem.characteristicTime = 2.0;
    testproblem.characteristicDensity = 5.0;
    testConstants.electron_Mobility = 1.0;
    testproblem.thermalVoltage = 1.0;

    carrierType = NegativeCharge;
    testElectrons.ReinitializeTestWithScaling(carrierType, testgrid2,
            testproblem, testConstants);
// set n(x,0) = sin(x)
    testElectrons.CarrierICFunction = &DiffusionWithScalingInitial;

    // set n_x(0,t) = 2 exp(-4t) , n_x(pi,t) = 2 exp(-4t)
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Robin;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction;
    testElectrons.CarrierRightBCFunction = &DiffusionWithScalingAndNeumannRightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &ZeroFunction;

    // Now set C(x) = 0 => E = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &DiffusionWithScalingU;
    testElectrons.TrueSolutionQ = &DiffusionWithScalingQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrierWithScaling();

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;

    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }


    // output results
//  testElectrons.outputTestCarrier(21);

///////////////////////////////////////////////////////////////////////////////
// Drift Diffusion Problem
///////////////////////////////////////////////////////////////////////////////

    cout << "\t ************************* \n \n";
    cout << "\t Drift-Diffusion w/ Implicit Time Stepping\n" << endl;
    cout << "\t ************************* \n \n";
    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = M_PI;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    carrierType = NegativeCharge;
    testElectrons.initializeTest(carrierType, testgrid2, testproblem, testConstants);



    cout << "27.) Electrons: Dirichlet Homogeneous\n\n";

// set n(x,0) = exp(-0.5x) sin(x)
    testElectrons.CarrierICFunction = &EDD_Dirichlet_Homo_IC;

    // set n(0,t) = 0, n(pi,t) = 0
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Dirichlet;
    testElectrons.CarrierLeftBCFunction = &ZeroFunction;
    testElectrons.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &EDD_Dirichlet_Homo_SolutionU;
    testElectrons.TrueSolutionQ = &EDD_Dirichlet_Homo_SolutionnQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }
    // output results
//  testElectrons.outputTestCarrier(27);


    cout << "28.)  Holes: Dirichlet Homogeneous\n\n";

// set p(x,0) = exp(0.5x) sin(x)
    testHoles.CarrierICFunction = &HDD_Dirichlet_Homo_IC;

    // set p(0,t) = 0, p(pi,t) = 0
    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Dirichlet_Homo_SolutionU;
    testHoles.TrueSolutionQ = &HDD_Dirichlet_Homo_SolutionnQ;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
///  testHoles.outputTestCarrier(28);

    cout << "29.) Electrons: Dirichlet Non-Homogeneous\n\n";

// set n(x,0) = exp(-0.5x) cos(x)
    testElectrons.CarrierICFunction = &EDD_Dirichlet_Non_Homo_IC;

    // set n(0,t) = exp(-1.25t) , n(pi,t) = -e(-0.5pi)exp(-1.25t)
    testElectrons.CarrierLeftBCType = Dirichlet;
    testElectrons.CarrierRightBCType = Dirichlet;
    testElectrons.CarrierLeftBCFunction = &EDD_Dirichlet_Non_Homo_LeftBC;
    testElectrons.CarrierRightBCFunction = &EDD_Dirichlet_Non_Homo_RightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testElectrons.PotentialLeftBCFunction = &ZeroFunction;
    testElectrons.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testElectrons.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testElectrons.TrueSolutionU = &EDD_Dirichlet_Non_Homo_SolutionU;
    testElectrons.TrueSolutionQ = &EDD_Dirichlet_Non_Homo_SolutionnQ;

    // set up the test
    testElectrons.setImplicitTestCarrier();

    // run the test
    testPassed = testElectrons.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
// testElectrons.outputTestCarrier(29);

    cout << "30.)  Holes: Dirichlet Non-Homogeneous\n\n";

    // set p(x,0) = exp(0.5x)cos(x)
    testHoles.CarrierICFunction = &HDD_Dirichlet_Non_Homo_IC;

    // set p(0,t) = exp(-1.25t) , n(pi,t) = -e(0.5pi)exp(-1.25t)
    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &HDD_Dirichlet_Non_Homo_LeftBC;
    testHoles.CarrierRightBCFunction = &HDD_Dirichlet_Non_Homo_RightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Dirichlet_Non_Homo_SolutionU;
    testHoles.TrueSolutionQ = &HDD_Dirichlet_Non_Homo_SolutionnQ;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
// testHoles.outputTestCarrier(30);

///////////////////////////////////////////////////////////////////////////////
// Drift Diffusion Problem With MIXED Boundary Conditions
///////////////////////////////////////////////////////////////////////////////
// Refer to documentation for BC and IC or testUtilities.cpp

    cout << " Left Neumann Non-Homogeneous\n";
    cout << " Right Dirichlet Homogeneous \n\n";

    cout << "31.) Holes: \n\n";

    testHoles.CarrierICFunction = HDD_Dirichlet_Homo_IC;

    testHoles.CarrierLeftBCType = Robin;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &HDD_Non_Homo_Neumann_Homo_Dir_LeftBC;
    testHoles.CarrierRightBCFunction = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Dirichlet_Homo_SolutionU;
    testHoles.TrueSolutionQ = &HDD_Dirichlet_Homo_SolutionnQ;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();
    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }
    // output results
    //testHoles.outputTestCarrier(31);

    cout << " Left Robin Non-Homogeneous \n";
    cout << " Right Dirichlet Non-Homogeneous\n\n";

    cout << "32.) Holes: \n\n";

    testHoles.CarrierICFunction = HDD_Dirichlet_Non_Homo_IC;

    testHoles.CarrierLeftBCType = Robin;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &HDD_Non_Homo_Robin_Non_Homo_Dir_LeftBC;
    testHoles.CarrierRightBCFunction = &HDD_Dirichlet_Non_Homo_RightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Dirichlet_Non_Homo_SolutionU;
    testHoles.TrueSolutionQ = &HDD_Dirichlet_Non_Homo_SolutionnQ;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();

    // see if it passed
    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
//  testHoles.outputTestCarrier(32);

    cout << " Left Dirichlet Homogeneous\n";
    cout << " Right Neumann Non-Homogeneous\n\n";

    cout << "33.) Holes: \n\n";

    testHoles.CarrierICFunction = HDD_Dirichlet_Homo_IC;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Robin;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &HDD_Homo_Dir_Non_Homo_Neumann_RightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Dirichlet_Homo_SolutionU;
    testHoles.TrueSolutionQ = &HDD_Dirichlet_Homo_SolutionnQ;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();

    // see if it passed

    if(testPassed)
    {
        cout << "\tPassed! \n\n";
        numTestsPassed++;
        numTestsTotal++;
    }
    else
    {
        cout << "\tFailed! \n\n";
        numTestsTotal++;
    }

    // output results
    //testHoles.outputTestCarrier(33);


    cout << " Left Dirichlet Homogeneous\n";
    cout << " Right Robin Non-Homogeneous\n\n";

    cout << "34.) Holes: \n\n";

    testHoles.CarrierICFunction = &HDD_Right_Robin_IC;

    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Robin;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;
    testHoles.CarrierRightBCFunction = &HDD_Homo_Dir_Non_Homo_Robin_RightBC;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &MinusPiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

    // Set the True solutions
    testHoles.TrueSolutionU = &HDD_Right_Robin_Solution_U;
    testHoles.TrueSolutionQ = &HDD_Right_Robin_Solution_Q;

    // set up the test
    testHoles.setImplicitTestCarrier();

    // run the test
    testPassed = testHoles.runImplicitTestCarrier();

///////////////////////////////////////////////////////////////////////////////
// REACTIVE INTERFACE FLUXES DIFFUSION EQUATION
///////////////////////////////////////////////////////////////////////////////
// u_t = u_xx
    /*
    // Solves for each carrier to make sure working properly.
    cout << "\t \t DIFFUSION EQUATION W/ REACTIVE INTERFACES \n \n ";

    // E = 0 in drift diffusion equation

    cout << "CHECK RESULTS" << endl;

     //Set the initial conditions
    testReactiveFluxes.CarrierICFunction = &OneFunction;
    testReactiveFluxes.Carrier1ICFunction = &OneFunction;//PiFunction;
     testReactiveFluxes.Carrier2ICFunction = &OneFunction; // Two

    // Set the boundary conditions
     testReactiveFluxes.CarrierLeftBCType = Dirichlet;
     testReactiveFluxes.CarrierRightBCType = Robin;
     testReactiveFluxes.CarrierLeftBCFunction = &OneFunction;
     testReactiveFluxes.CarrierRightBCFunction = &ZeroFunction;

    testReactiveFluxes.Carrier1LeftBCType = Robin;
    	testReactiveFluxes.Carrier1RightBCType =  Dirichlet;
     testReactiveFluxes.Carrier1LeftBCFunction = &ZeroFunction;
     testReactiveFluxes.Carrier1RightBCFunction = &OneFunction; // Pi

    testReactiveFluxes.Carrier2LeftBCType = Robin;
    	testReactiveFluxes.Carrier2RightBCType =  Dirichlet;
     testReactiveFluxes.Carrier2LeftBCFunction = &ZeroFunction;
     testReactiveFluxes.Carrier2RightBCFunction = &OneFunction; // Two

     // set Phi(0) = 0, Phi(pi) = 0
     testReactiveFluxes.PotentialLeftBCFunction = &ZeroFunction;
     testReactiveFluxes.PotentialRightBCFunction = &ZeroFunction;

     // Now set C(x) = 0 => E = 0
     testReactiveFluxes.DopingProfileFunction = &ZeroFunction;

    testReactiveFluxes.TrueSolutionQ = &ReactiveFluxesQ;

     // set up the test
     testReactiveFluxes.setTestReactiveFluxes();

     // run the test
     testPassed = testReactiveFluxes.runImplicitExplicitTestReactiveFluxes();
     //testPassed = testReactiveFluxes.runImplicitTestReactiveFluxes();

     // see if it passed
     if(testPassed)
     {
    	cout << "\tPassed! \n\n";
    	numTestsPassed++;
    	numTestsTotal++;
     }
     else
     {
    	cout << "\tFailed! \n\n";
    	numTestsTotal++;
     }

     // output results
    //  testHoles.outputTestCarrier(34);
    */


    cout << "\n\nTesting Complete:" << endl;

    cout << "\n \t\t\tA total of " << numTestsPassed << "/" << numTestsTotal
         << " tests Passed.\n\n";


    return 0;
}



int ddpTestingMakeInputFile(void)
{
    ofstream inputfile;
    inputfile.open("ddp-test.ini");
    inputfile <<
              "[computational]\n \
		timeStepFactor = 1.0 \n \
		maxOrderMX = 1 \n \
		maxOrderDG = 1 \n \
		numElements = 150 \n \
		numBoundaryElements = 0 \n \
		boundaryLayerWidth = 0.0 \n \
		numTimeStamps = 0 \n \
		GaussLegendreNumPoints = 16\n \
		\n \n \
		[physical] \n \
		xLeftEndPoint = 0.0\n \
		xRightEndPoint = 1.0\n \
		characteristicLength	= 1.0\n \
		IVP_Type = Original \n \
		timeInitial = 0.0\n \
	  timeFinal = 1.00\n \
		characteristicTime 	= 1.0\n \
		appliedBias = 0.0 \n \
		builtInBias	= 0.0 \n \
		semiCondRelativePerm = 1.0\n \
		electrolyteRelativePerm = 1.0\n \
		photonFlux = 0.0\n\
		absorptionCoeff = 0.0\n\
		characteristicDensity = 1.0\n \
		intrinsicDensity = 1.0\n \
		\n \n \
		[illuminationStatus] \n \
		illumination = Off \n \
		\n \n \
		[couplingStatus]\n \
		couplingToPoisson = Off \n \
		\n \n \
  	[electrons] \n \
		Mobility = 1.0 \n \
		ChargeSign = Negative \n \
		recombinationTime = 1.0 \n \
		TransferRate = 0.0 \n \
		RecombinationVelocity = 0.0 \n \
		\n \n \
		[holes] \n \
		Mobility = 1.0 \n \
		ChargeSign = Positive \n \
		recombinationTime = 1.0 \n \
		TransferRate = 0.0 \n \
		RecombinationVelocity = 0.0 \n \
		\n \n \
  	[reductants] \n \
		Mobility = 0.0 \n \
		ChargeSign = Negative \n \
  	[oxidants] \n \
		Mobility = 0.0 \n \
		ChargeSign = Positive\n";
    inputfile.close();
    return 0;
}



