#include "test.hpp"

int ddpTestingMakeInputFile(int order, int numEl);
int Compute_Poisson_Errors(ddpDenseVector_type & errors, int runs);
int Compute_Carrier_Errors(ddpDenseVector_type & errors, int runs);

int main()
{
    enum Solver {MX, LDG};

    // SOLVER TYPE
    Solver solver_type = LDG;

    int N = 10;
    int start = 2;
    int order = 0;
    int numElem = 2;
    int runs = 0;
    ddpDenseVector_type const_errors(3);
    ddpDenseVector_type linear_errors(3);
    ddpDenseVector_type quad_errors(3);

    ofstream prt;
    if(MX == solver_type)
    {
        prt.open("Mixed_Errors.dat");

        for(int i = start; i < N; i++)
        {
            numElem = (int)std::pow(2.0,i);

            // consts
            order = 0;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Poisson_Errors(const_errors, runs);

            // linear
            order = 1;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Poisson_Errors(linear_errors, runs);

            // quads
            order = 2;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Poisson_Errors(quad_errors, runs);

            prt << (const_errors(2)) << "\t"
                << (const_errors(0)) << "\t"
                << (const_errors(1)) << "\t"
                << (linear_errors(0)) << "\t"
                << (linear_errors(1)) << "\t"
                << (quad_errors(0)) << "\t"
                << (quad_errors(1)) << endl;
        } // end for
        prt.close();
    } // end if
    else
    {
        prt.open("LDG_Errors.dat");

        for(int i = start; i < N; i++)
        {
            numElem = (int)std::pow(2.0,i);

            // const
            order = 0;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Carrier_Errors(const_errors, runs);

            // linear
            order = 1;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Carrier_Errors(linear_errors, runs);

            // quads
            order = 2;
            ddpTestingMakeInputFile(order, numElem);
            Compute_Carrier_Errors(quad_errors, runs);

            prt << (const_errors(2)) << "\t"
                << (const_errors(0)) << "\t"
                << (const_errors(1)) << "\t"
                << (linear_errors(0)) << "\t"
                << (linear_errors(1)) << "\t"
                << (quad_errors(0)) << "\t"
                << (quad_errors(1)) << endl;
        } // end for
        prt.close();
    } // end if




    return 0;
}

int Compute_Poisson_Errors(ddpDenseVector_type & errors, int runs)
{
    ddpProblemInfo_type testproblem;
    ddpCarrierConstants_type testConstants;


    ddpDomain_type testdomain2;
    ddpGrid_type testgrid2;
    readInput(testdomain2, testgrid2, testproblem, testConstants, "ddp-test.ini");

    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = M_PI;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);


    // Test For Poisson's Equation
    ddpTestChargeCarrier_type carrierType = PositiveCharge; // doesnt matter
    ddpTest_type testPoisson;
    testPoisson.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);
    //ddpPrintGrid(testgrid2);

    testPoisson.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    // set n = p = 0
    testPoisson.Carrier1ICFunction = &ZeroFunction;
    testPoisson.Carrier2ICFunction = &ZeroFunction;

    // set C = sin(x)
    testPoisson.DopingProfileFunction1 = &SineFunction; //Poisson_L2_RHS;
    testPoisson.DopingProfileFunction2 = &ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testPoisson.PotentialLeftBCFunction = &ZeroFunction;
    testPoisson.PotentialRightBCFunction = &ZeroFunction;

    /// set true solutions
    testPoisson.TrueSolutionPotential = &SineFunction; //Poisson_L2_Potential;
    testPoisson.TrueSolutionElecField = &MinusCosineFunction; //Poisson_L2_ElecField;


    testPoisson.setTestPoisson();

    testPoisson.runTestPoissonL2Error(errors);

    errors(2) = testgrid2.DeltaxMin;

    testPoisson.outputTestPoisson(runs);

    return 0;
}

int Compute_Carrier_Errors(ddpDenseVector_type & errors, int runs)
{
    ddpProblemInfo_type testproblem;
    ddpCarrierConstants_type testConstants;

    ddpDomain_type testdomain2;
    ddpGrid_type testgrid2;
    readInput(testdomain2, testgrid2, testproblem, testConstants, "ddp-test.ini");

    testdomain2.LeftEndPoint = 0.0;
    testdomain2.RightEndPoint = M_PI;
    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    // Test For Poisson's Equation
    ddpTestChargeCarrier_type carrierType = PositiveCharge; // doesnt matter
    ddpTest_type testHoles;
    testHoles.initializeTest(carrierType, testgrid2, testproblem, testConstants);

    ddpMakeUniformGrid(testproblem, testdomain2, testgrid2);

    testHoles.initializeTest(carrierType, testgrid2, testproblem, testConstants);

// transients
    testHoles.CarrierICFunction = &HDD_Dirichlet_Homo_IC;
    testHoles.CarrierLeftBCType = Dirichlet;
    testHoles.CarrierRightBCType = Dirichlet;
    testHoles.CarrierLeftBCFunction = &ZeroFunction;//Non_Homo_Neumann_Homo_Dir_LeftBC;
    testHoles.CarrierRightBCFunction = &OneFunction; //ZeroFunction;

    // set Phi(0) = 0, Phi(pi) = 0
    testHoles.PotentialLeftBCFunction = &ZeroFunction;
    testHoles.PotentialRightBCFunction = &OneFunction; //PiFunction;

    // Now set C(x) = 0
    testHoles.DopingProfileFunction = &ZeroFunction;

//	testHoles.testCarrier.setTestGenerationRHS(testgrid2, &DD_L2_RHS); //&SineFunction);

//  testHoles.TrueSolutionU = &HDD_Dirichlet_Homo_SolutionU;
//  testHoles.TrueSolutionQ = &HDD_Dirichlet_Homo_SolutionnQ;

    testHoles.TrueSolutionU = &HDD_SteadyState_Dir_U; //Poisson_L2_Potential;
    testHoles.TrueSolutionQ = &HDD_SteadyState_Dir_Q; //Poisson_L2_ElecField;

    // set up the test
    testHoles.setSteadyStateTestCarrier();

    // run the test
    testHoles.runTestCarrierL2Error(errors);

    // set up the test
//  testHoles.setImplicitTestCarrier();

    // run the test
//  testHoles.runImplicitTestCarrierL2Error(errors);


    errors(2) = testgrid2.DeltaxMin;

    testHoles.outputTestCarrier(runs);

    /*
    	// steady state
    	testHoles.CarrierICFunction = &ZeroFunction;
    	testHoles.CarrierLeftBCType = Dirichlet;
      testHoles.CarrierRightBCType = Robin;
      testHoles.CarrierLeftBCFunction = &ZeroFunction;
      testHoles.CarrierRightBCFunction = &DD_L2_BC;

    	// set Phi(0) = 0, Phi(pi) = 0
      testHoles.PotentialLeftBCFunction = &ZeroFunction;
     	testHoles.PotentialRightBCFunction = &ZeroFunction;

      // Now set C(x) = 0
      testHoles.DopingProfileFunction = &ZeroFunction;

    	testHoles.testCarrier.setTestGenerationRHS(testgrid2, &DD_L2_RHS);

     	testHoles.TrueSolutionU = &DD_L2_Density;
      testHoles.TrueSolutionQ = &DD_L2_Current;


      // set up the test
      testHoles.setSteadyStateTestCarrier();

      // run the test
      testHoles.runTestCarrierL2Error(errors);

    	errors(2) = testgrid2.DeltaxMin;
    	testHoles.outputTestCarrier(runs);
    */

    return 0;
}


int ddpTestingMakeInputFile(int order, int num_Elem)
{
    ofstream inputfile;
    inputfile.open("ddp-test.ini");
    inputfile <<
              "[computational]\n \
		timeStepFactor = 1.0 \n \
		maxOrderMX = " << order << " \n \
		maxOrderDG = " << order << " \n \
		numElements = " << num_Elem << " \n \
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

