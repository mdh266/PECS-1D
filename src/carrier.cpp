#include "../include/carrier.hpp"

///////////////////////////////////////////////////////////////////////////////
//  Initialization funcitons
///////////////////////////////////////////////////////////////////////////////

int
ddpCarrier_type::
initialize(const ddpGrid_type & grid,
           const ddpProblemInfo_type & problem,
           const ddpCarrierConstants_type & carrierConstants,
           char const * nameOfCarrier
          )
{

    //////////////////////////////////////////////////////////////////////
    // Store all the physical constants depending on the sign of the
    // chosen charge.
    /////////////////////////////////////////////////////////////////////


    std::string CarrierName(nameOfCarrier);
    std::string electronName("electrons");
    std::string holeName("holes");
    std::string reductantName("reductants");
    std::string oxidantName("oxidants");


    if(CarrierName.compare(electronName)  == 0)
    {
        // get electrons semiconductor mobility
        carrierProps.Mobility = 	carrierConstants.electron_Mobility;

        // get electrons semiconductor Diffusivity by einsteins relation
        carrierProps.Diffusivity = 	carrierProps.Mobility *
                                    problem.thermalVoltage;

        // get other properties of carrier and scale them.
        carrierProps.ChargeSign = carrierConstants.electron_ChargeSign;

        // get recombination time
        carrierProps.RecombinationTime =
            carrierConstants.electron_RecombinationTime /
            problem.characteristicTime;

        // get the transfer rate for the electrons and scale it
        carrierProps.TransferRate =
            (carrierConstants.electron_TransferRate *
             problem.characteristicTime *
             problem.characteristicDensity /
             problem.characteristicLength);


        // get the transfer rate for the electrons and scale it
        carrierProps.RecombinationVelocity =
            (carrierConstants.electron_RecombinationVelocity *
             problem.characteristicTime /
             problem.characteristicLength);


        // set the material type
        carrierProps.Material = Semiconductor;

        // set and scale the absorption coefficent
        GenerationFunction.AbsorpCoeff =  problem.Absorption_Coeff;

        // set and scale the photon flux
        GenerationFunction.PhotonFlux = (problem.Photon_Flux
                                         * problem.characteristicTime
                                        ) / problem.characteristicDensity;


//		cout << "Scaled Photon Flux = " << GenerationFunction.PhotonFlux << endl;

        // set and scale the absorption coefficient for the expoential
        GenerationFunction.ScaledAbsorpCoeff = problem.Absorption_Coeff
                                               * problem.characteristicLength;

//		cout << "Scaled Absorption = " << GenerationFunction.ScaledAbsorpCoeff << endl;

        // Set the intrinsic denisty
        carrierProps.IntrisincDensity = problem.intrinsicDensity
                                        / problem.characteristicDensity;
    }
    else if(CarrierName.compare(holeName) == 0 )
    {
        // get holes semiconductor mobility
        carrierProps.Mobility = 	carrierConstants.hole_Mobility;

        // get holes semiconductor Diffusivity by einsteins relation
        carrierProps.Diffusivity = carrierProps.Mobility *
                                   problem.thermalVoltage;

        // get the rest of the carrier properties
        carrierProps.ChargeSign = carrierConstants.hole_ChargeSign;

        // get the recombination time
        carrierProps.RecombinationTime =
            carrierConstants.hole_RecombinationTime /
            problem.characteristicTime;

        // GEt the hole transfer rate and scale it
        carrierProps.TransferRate =
            (carrierConstants.hole_TransferRate *
             problem.characteristicTime *
             problem.characteristicDensity /
             problem.characteristicLength);

        carrierProps.RecombinationVelocity =
            (carrierConstants.hole_RecombinationVelocity *
             problem.characteristicTime) /
            problem.characteristicLength;


        // set the material type
        carrierProps.Material = Semiconductor;

        // set and scale the absorption coefficent
        GenerationFunction.AbsorpCoeff =  problem.Absorption_Coeff;

        // set and scale the photon flux
        GenerationFunction.PhotonFlux = (problem.Photon_Flux
                                         * problem.characteristicTime)
                                        / problem.characteristicDensity;

        // set and scale the absorption coefficient for the expoential
        GenerationFunction.ScaledAbsorpCoeff = problem.Absorption_Coeff
                                               * problem.characteristicLength;

        // Set the intrinsic denisty
        carrierProps.IntrisincDensity = problem.intrinsicDensity
                                        / problem.characteristicDensity;

    }
    else if(CarrierName.compare(reductantName) == 0 )
    {

        // get reductants mobility
        carrierProps.Mobility = 	carrierConstants.reductant_Mobility;

        // get electrons semiconductor Diffusivity by einsteins relation
        carrierProps.Diffusivity = 	carrierProps.Mobility *
                                    problem.thermalVoltage;

        // get other properties of carrier and scale them.
        carrierProps.ChargeSign = carrierConstants.reductant_ChargeSign;

        // set the material type
        carrierProps.Material = Electrolyte;

        // set and scale the absorption coefficent
        GenerationFunction.AbsorpCoeff = 0.0;

        // set and scale the photon flux
        GenerationFunction.PhotonFlux = 0.0;

        // Set the intrinsic denisty
        carrierProps.IntrisincDensity = 0.0;

    }
    else if(CarrierName.compare(oxidantName) == 0 )
    {
        // get holes semiconductor mobility
        carrierProps.Mobility = 	carrierConstants.oxidant_Mobility;

        // get holes semiconductor Diffusivity by einsteins relation
        carrierProps.Diffusivity = carrierProps.Mobility *
                                   problem.thermalVoltage;

        // get the rest of the carrier properties
        carrierProps.ChargeSign = carrierConstants.oxidant_ChargeSign;

        // set the material type
        carrierProps.Material = Electrolyte;

        // set and scale the absorption coefficent
        GenerationFunction.AbsorpCoeff = 0.0;

        // set and scale the photon flux
        GenerationFunction.PhotonFlux = 0.0;

        // Set the intrinsic denisty
        carrierProps.IntrisincDensity = 0.0;
    }
    else
    {
        cout << "\n \nERROR: No carrier of that name. You may use: \n \
			electrons \n \
			holes\n" << endl;

        assert(false);
    }


    // variables to scale and unscale the variables
    carrierProps.Scale4Mobility = (problem.thermalVoltage *
                                   problem.characteristicTime /
                                   (problem.characteristicLength *
                                    problem.characteristicLength));

    carrierProps.Scale4Diffusivity = (problem.characteristicTime /
                                      (problem.characteristicLength *
                                       problem.characteristicLength));


    carrierProps.RecomboScale = 1.0; //problem.characteristicTime;

    // Using Einstein relationship the diffusivity we just have
    // J = \omega ( E u - du/dx)
    // where \oemga = Scale4Mobility * Mobility
    // 							= Scale4Diffusivity * Diffusivity

    double Scale1 = 	carrierProps.Mobility	* carrierProps.Scale4Mobility;
    double Scale2 = 	carrierProps.Diffusivity * carrierProps.Scale4Diffusivity;

    // checks to  make sure they are the same.. rounding error issue
    assert(abs(Scale1 - Scale2) < 1.0e-7);

    carrierProps.CurrentScale = carrierProps.Mobility *
                                carrierProps.Scale4Mobility;

    carrierProps.LengthScale = problem.characteristicLength;
    // set the constants for charge sign this is used to determine the force
    // on the electric field on the charge and the sign the right hand side
    // of poissons equation when fully coupled.
    if(carrierProps.ChargeSign == Negative)
    {
        carrierProps.Sign4Force 	= -1.0;
        carrierProps.Sign4Poisson = 1.0;
    }
    else
    {
        carrierProps.Sign4Force = 1.0;
        carrierProps.Sign4Poisson = -1.0;
    }

    // TODO:  These could be carrier independent and passed around
    // as Data structures IF the memory becomes and isuse or having
    // to build them each for a carrier becomes an issue timewise

    // make the quadrature weights and points vectors
    ddpMakeWeightsAndPoints(grid,
                            problem,
                            weightsSparse,
                            PTSSparse,
                            PTSDense,
                            weightsDense);


    // Set the left end point for the boundary condition.  This
    // needs the domain to set it and thats why it is done here.
    // The domain could be [0,1] or [-1, 1]
    carrierProps.BCLeft.xLeftEndPoint = (grid.ElementList[0]).Left;


    // Make the bijections from the local DOF to global DOF
    ddpMakeAllBijections(grid, Bijections);


    // Make the global VandeMonde Matrice and Flux Matrices
    // DOF -> PT Values
    ddpMakeVandeMondeMatrices(grid,
                              problem,
                              Bijections,
                              VandeMondeMatrices,
                              DGFluxMatrices);


    // Make all the carrier Matrices
    // Matrices related to DOFs that WONT change in time.
    ddpMakeCarrierProperies(problem,
                            grid,
                            weightsSparse,
                            Bijections,
                            VandeMondeMatrices,
                            DGFluxMatrices,
                            carrierProps);

    // Get the number of DOF for DG problem and the number of
    // cell boundary points
    int DGDof = Bijections.DGForward.size();
    int numEndPoints = grid.NumElementsNoGhost + 1;

    // Set the RHS for the Generation
    setGenerationRHS(grid, problem.IlluminationStatus);


    carrierProps.MaxMobility = carrierProps.CurrentScale;
    carrierProps.MaxDiffusivity = carrierProps.CurrentScale;


    return 0;
}


int
ddpCarrier_type::
setSolver(ddpGrid_type const & grid,
          ddpProblemInfo_type const & problem)
{
    // Get the number of DOF for DG problem and the number of
    // cell boundary points
    int DGDof = Bijections.DGForward.size();
    int numEndPoints = grid.NumElementsNoGhost + 1;
    int numElements = numEndPoints-1;
    int Rows = carrierProps.TotalFluxFromBCRHS.rows();

    // Allocate memory for carrierState Dofs and uNext
    carrierState.uDof = ddpDenseVector_type::Zero(DGDof);
    carrierState.qDof = ddpDenseVector_type::Zero(DGDof);
    carrierState.qDof_prev = ddpDenseVector_type::Zero(DGDof);

    // allocate memory for FOR RHS Vectors of dx/dt = Ax + b
    carrierState.BC_Dir_Input = ddpDenseVector_type::Zero(2);
    carrierState.BC_Robin_Input = ddpDenseVector_type::Zero(2);
    carrierState.RHSFromRecombination = ddpDenseVector_type::Zero(DGDof);
    carrierProps.Dir_RHS = ddpSparseMatrix_type(2*Rows,2);
    carrierProps.Robin_RHS = ddpSparseMatrix_type(2*Rows,2);
    carrierState.LDG_RHS = ddpDenseVector_type(2*Rows);
    carrierProps.Penalty_RHS = ddpSparseMatrix_type(2*Rows,2);
    carrierState.BigSolution = ddpDenseVector_type(2*Rows);
    carrierState.OldBigSolution = ddpDenseVector_type(2*Rows);

    // allocate memory for work space vectors
    carrierState.ElecPTS = ddpDenseVector_type::Zero(PTSDense.size());
    carrierState.ElecPTSTimesWeights = ddpDenseVector_type::Zero(PTSDense.size());
    carrierState.DriftTerm = ddpDenseVector_type::Zero(2*Rows);

    // allocate memory for the work space matrices
    carrierState.ElecTimesWeightsMatrix =
        ddpDiagonalMatrix_type(weightsDense.size() );

    carrierState.PE = ddpSparseMatrix_type(Rows, Rows);

    // construct LDG matrices
    carrierProps.B2 = -carrierProps.StiffQFromU
                      + carrierProps.FluxRightQFromU
                      - carrierProps.FluxLeftQFromU;

    carrierProps.B1 = -carrierProps.StiffUFromQ
                      + carrierProps.FluxRightUFromQ
                      - carrierProps.FluxLeftUFromQ;

    // TODO : JUST FOR MONOLITHIC CASE

    // construct BC Sparse matrices for Dirichlet data
    for(int k = 0; k < carrierProps.TotalFluxFromBCRHS.outerSize(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.TotalFluxFromBCRHS,k);
                it; ++it)
        {
            carrierProps.Dir_RHS.insert(it.row(), k) = it.value();
        }
    }

    // construct BC Sparse matrices for Robin data
    for(int k = 0; k < carrierProps.TotalFluxFromBCRHS.outerSize(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.TotalFluxFromBCRHS,k);
                it; ++it)
        {
            carrierProps.Robin_RHS.insert(Rows + it.row(), k) = it.value();
        }
    }

    // construct Penalty BC Sparse matrices for Dirichlet data
    for(int k = 0; k < carrierProps.TotalFluxFromBCPenalty.outerSize(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.TotalFluxFromBCPenalty,k);
                it; ++it)
        {
            carrierProps.Penalty_RHS.insert(Rows + it.row(), k) = it.value();
        }
    }


    carrierProps.MassUToSystem = ddpSparseMatrix_type(2*Rows, Rows);
    // construct matrix to convert MassU * u^{n} to right size
    for(int k = 0; k < carrierProps.MassU.rows(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.MassU,k);
                it; ++it)
        {
            carrierProps.MassUToSystem.insert(Rows + it.row(), k) = it.value();
        }
    }
    //cout << "MassUToSystem = \n" << carrierProps.MassUToSystem << endl;

    carrierProps.RecomboToSystem = ddpSparseMatrix_type(2*Rows,Rows);
    for(int k = 0; k < carrierProps.MassU.rows(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.MassU,k);
                it; ++it)
        {
            carrierProps.RecomboToSystem.insert(Rows + it.row(), k) = 1.0;
        }
    }

    carrierProps.DriftToSystem = ddpSparseMatrix_type(2*Rows,Rows);
    for(int k = 0; k < carrierProps.MassU.rows(); k++)
    {
        for(ddpSparseMatrix_type::InnerIterator it(carrierProps.MassU,k);
                it; ++it)
        {
            carrierProps.DriftToSystem.insert(it.row(), k) = 1.0;
        }
    }
    /*
    	cout << "DIR_RHS = \n " << carrierProps.Dir_RHS << endl;
    	cout << "Penalty_RHS = \n" << carrierProps.Penalty_RHS << endl;
    	cout << "ROBIN_RHS = \n " << carrierProps.Robin_RHS << endl;
    	cout << "MassU = \n" << carrierProps.MassU << endl;
    	cout << "MassUToSystem = \n" << carrierProps.MassUToSystem << endl;
    	cout << "RecomboToSystem = \n" << carrierProps.RecomboToSystem << endl;
    	cout << "DrftToSystem = \n" << carrierProps.DriftToSystem << endl;

    	assert(false);
    */

    return 0;
}

int
ddpCarrier_type::
setTestGenerationRHS(ddpGrid_type const & grid,
                     double (*GenerationFunction)(const double & x) )
{

    // does int_{I} G f(x) dx
    ddpBijFlag_type BijFlag = DG;
    ddpProjectRHSFunction(GenerationFunction,
                          grid,
                          Bijections.DGForward,
                          VandeMondeMatrices.globalVandeMondeDG,
                          Bijections.PTForward,
                          BijFlag,
                          weightsSparse,
                          carrierState.RHSFromGeneration);


    int ROWS = carrierState.RHSFromGeneration.size();
    carrierState.Generation_RHS = ddpDenseVector_type::Zero(2*ROWS);


    for(int k = 0; k < ROWS; k++)
    {
        carrierState.Generation_RHS(ROWS + k) = carrierState.RHSFromGeneration(k);
    }

    return 0;
}

int
ddpCarrier_type::
setGenerationRHS(ddpGrid_type const & grid,
                 const ddpIlluminationStatus_type & IlluminationStatus)
{
    // Create the RHS vector for generation for this specific carrier
    int DGDof = Bijections.DGForward.size();
    carrierState.RHSFromGeneration = ddpDenseVector_type::Zero(DGDof);

    ddpDenseVector_type temp_Generation = ddpDenseVector_type::Zero(DGDof);

    if(Illuminated == IlluminationStatus)
    {
        ddpBijFlag_type BijFlag = DG;

        //	cout << "illuminated" << endl;
        // Project this carrier's initial condition function onto the DG basis
        // and store result in this carriers.UDof

        // does M^{-1} * int_{I} v G(x) dx
        ddpProjectFunction(&ddpGenerationFunction_type::GenFun,
                           GenerationFunction,
                           grid,
                           Bijections.DGForward,
                           VandeMondeMatrices.globalVandeMondeDG,
                           Bijections.PTForward,
                           BijFlag,
                           weightsSparse,
                           carrierProps.MassU,
                           temp_Generation);
    }

    // makes it so its M * M^{-1} * int_{I}  v G(x) dx
    carrierState.RHSFromGeneration = carrierProps.MassU *  temp_Generation;

    int ROWS = carrierState.RHSFromGeneration.size();
    carrierState.Generation_RHS = ddpDenseVector_type::Zero(2*ROWS);

    for(int k = 0; k < ROWS; k++)
    {
        carrierState.Generation_RHS(ROWS + k) = carrierState.RHSFromGeneration(k);
    }
//	cout << "Absoprtion Coeff = " << GenerationFunction.AbsorpCoeff << endl;
//	cout << "Solar Flux = " << GenerationFunction.PhotonFlux << endl;
    //:"	cout << carrierState.Generation_RHS << endl;
    return 0;

}


int
ddpCarrier_type::
assembleLDGMatrices(ddpGrid_type const & grid,
                    ddpProblemInfo_type const & problem)
{
    // Create flux matrices to pick up appropriate
    // boundary conditions
    // Theses matrices will have different values at the
    // boundaryies and interface depending on whether what kind of
    // what kind of boundary problem we are solving.

    // This assures that the direction of the LDG Fluxes
    // is appropriate for the boundary conditions.
    // If both boundary conditions are Dirichlet, it doesnt matter

    // Beta = 0.0 => Central fluxes
    if(Dirichlet == carrierProps.BCLeft.BCType)
    {
        if(Dirichlet == carrierProps.BCRight.BCType)
        {
            carrierProps.Beta = -1.0; //  changed from +1
            // \hat{q} = q^{-} 9/25/2014
            // beter result if \hat{q} following characteristic
        }
        else // RIGHT = NEUMANN OR ROBIN => \hat{q} = q^+
        {
            carrierProps.Beta = 1.0;
        }
    }
    else if(Dirichlet == carrierProps.BCRight.BCType)
    {
        if(Dirichlet == carrierProps.BCLeft.BCType)
        {
            carrierProps.Beta = 1.0; //1.0; //Doesnt matter
        }
        else // LEFT = NEUMANN OR ROBIN => \hat{q} = q^-
        {
            carrierProps.Beta = -1.0;
        }
    }
    else // IM NOT SURE, TAKE BETA = 1.0.. NOT SURE IF IT MAKES A DIFFERENCE
    {
        carrierProps.Beta = 1.0;
    }

    // set the penalty parameter for \tau * [u][v]
    carrierProps.Tau = 1.0; ///grid.DeltaxMin;

    ddpMakeDiffusiveFluxProperties(problem,
                                   grid,
                                   weightsSparse,
                                   Bijections,
                                   VandeMondeMatrices,
                                   DGFluxMatrices,
                                   carrierProps);

    return 0;
}


int
ddpCarrier_type::
setInitialConditions(ddpGrid_type & grid,
                     double (* InitialFunction) (const double & t) )
{
    // set this carriers initial condition function to
    // passsed function pointer InitialFunction
    InitialConditionsU = InitialFunction;

    ddpBijFlag_type BijFlag = DG;

    // Project this carrier's initial condition function onto the DG basis
    // and store result in this carriers.UDof
    ddpProjectFunction(InitialConditionsU,
                       grid,
                       Bijections.DGForward,
                       VandeMondeMatrices.globalVandeMondeDG,
                       Bijections.PTForward,
                       BijFlag,
                       weightsSparse,
                       carrierProps.MassU,
                       carrierState.uDof
                      );

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
// Functions To Set The Boundary Conditions
///////////////////////////////////////////////////////////////////////////////

//  The first thing that is done is to set the type of boundary condition
//  as to wheter its Dirichlet or Shottky-type. Testing functions are used for
//  solving mathematical problems and tests.  The others are set to do physical
//  problems. If a boundary condition function is not being employed we set it
//  to zero, so that we do not have to have an if statement within the
//  makeUdot function.


int
ddpCarrier_type::
setLeft_DirBC(double (* DopingProfile)(const double & x) )
{
    // get the mathematical type of boundary
    carrierProps.BCLeft.BCType = Dirichlet;

    // set the Density at infinity as the Ohmic value..
    // This is just a constant dirichlet value so its fine.
    carrierProps.BCLeft.OhmicValue =
        (*DopingProfile)(carrierProps.BCLeft.xLeftEndPoint);

    // set the allowance of Robin BC
    carrierProps.BCLeft.RobinValue = 0.0; // Automatically to zero

    // set the rest of the Boundary functions to zero function.
    carrierProps.BCLeft.TestDirichletValue = &ZeroFunction;
    carrierProps.BCLeft.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setRight_DirBC(double (* DopingProfile)(const double & x) )
{
    // get the mathematical type of boundary
    carrierProps.BCRight.BCType = Dirichlet;

    // set the Density at infinity as the Ohmic value..
    // This is just a constant dirichlet value so its fine.
    carrierProps.BCRight.OhmicValue = (*DopingProfile)(1.0);

    // set the allowance of Robin BC
    carrierProps.BCRight.RobinValue = 0.0; // Automatically to zero

    // set rest to zero function
    carrierProps.BCRight.TestDirichletValue = &ZeroFunction;
    carrierProps.BCRight.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setLeft_RobincBC()
{
    // get the mathematical type of boundary
    carrierProps.BCLeft.BCType = Robin;

    // set the allowance of Robin BC
    carrierProps.BCLeft.RobinValue = 0.0; // initially set to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCLeft.OhmicValue = 0.0;

    carrierProps.BCLeft.TestDirichletValue = &ZeroFunction;
    carrierProps.BCLeft.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setRight_RobinBC()
{
    // get the mathematical type of boundary
    carrierProps.BCRight.BCType = Robin;

    // set the allowance of Robin BC
    carrierProps.BCRight.RobinValue = 0.0; // initially to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCRight.OhmicValue = 0.0;

    carrierProps.BCRight.TestDirichletValue = &ZeroFunction;
    carrierProps.BCRight.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setLeft_RobinValue(const double & r)
{
    // Set the Left robin value
    carrierProps.BCLeft.RobinValue = r;

    return 0;
}

int
ddpCarrier_type::
setRight_RobinValue(const double & r)
{
    // Set the Right robin value
    carrierProps.BCRight.RobinValue = r;

    return 0;
}
// NOTE:  Now passing in function pointers to boundary conditions,
// 	  so function being passed in is a function of time.

int
ddpCarrier_type::
setLeftTestDirichletBC(double (* BCLeft_Input) (const double & t) )
{
    // get the mathematical type of boundary
    carrierProps.BCLeft.BCType = Dirichlet;

    // set the allowance of Robin BC
    carrierProps.BCLeft.RobinValue = 0.0; // Automatically to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCLeft.OhmicValue = 0.0;

    carrierProps.BCLeft.TestDirichletValue = BCLeft_Input;
    carrierProps.BCLeft.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setRightTestDirichletBC(double (* BCRight_Input) (const double & t) )
{
    // get the mathematical type of boundary
    carrierProps.BCRight.BCType = Dirichlet;

    // set the allowance of Robin BC
    carrierProps.BCRight.RobinValue = 0.0; // Automatically to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCRight.OhmicValue = 0.0;

    carrierProps.BCRight.TestDirichletValue = BCRight_Input;
    carrierProps.BCRight.TestRobinValue = &ZeroFunction;

    return 0;
}

int
ddpCarrier_type::
setLeftTestRobincBC(double (* BCLeft_Input) (const double & t) )
{
    // get the mathematical type of boundary
    carrierProps.BCLeft.BCType = Robin;

    // set the allowance of Robin BC
    carrierProps.BCLeft.RobinValue = 0.0; // Automatically to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCLeft.OhmicValue = 0.0;

    carrierProps.BCLeft.TestDirichletValue = &ZeroFunction;
    carrierProps.BCLeft.TestRobinValue = BCLeft_Input;

    return 0;
}

int
ddpCarrier_type::
setRightTestRobinBC(double (* BCRight_Input) (const double & t) )
{
    // get the mathematical type of boundary
    carrierProps.BCRight.BCType = Robin;

    // set the allowance of Robin BC
    carrierProps.BCRight.RobinValue = 0.0; // Automatically to zero

    // set boundary type to input function pointer and rest to zero
    carrierProps.BCRight.OhmicValue = 0.0;

    carrierProps.BCRight.TestDirichletValue = &ZeroFunction;
    carrierProps.BCRight.TestRobinValue = BCRight_Input;

    return 0;
}



int
ddpCarrier_type::
setEquilibriumDensity(double (* DopingProfile)(const double & x) )
{
    carrierProps.EquilibriumDensity = DopingProfile(1.0);
}

double
ddpCarrier_type::
getEquilibriumDensity()
{
    return carrierProps.EquilibriumDensity;
}

double
ddpCarrier_type::
getTransferRate()
{
    return carrierProps.TransferRate;
}


double
ddpCarrier_type::
getInterfaceDensity()
{
    // Get the carrier density point values
    ddpDenseVector_type DensityPTValues =
        VandeMondeMatrices.globalVandeMondeDG
        * carrierState.uDof;

    // Interface is on the right
    if(Semiconductor == carrierProps.Material)
    {
//		cout << "Semiconductor" << endl;
        return DensityPTValues(DensityPTValues.size() -1);
    }
    else // interface is on the right
    {
//		cout << "Electrolyte" << endl;
        return DensityPTValues(0);
    }
}

/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// STEADY STATE SOLVERS
///////////////////////////////////////////////////////////////////////////////

int
ddpCarrier_type::
Solve(ddpDenseVector_type const & ElecFieldDof,
      ddpProblemInfo_type const & problem,
      ddpGrid_type const & grid,
      double const & timeCurrent)
{
    /////////////////////////
    // MONOLITHIC SOLVER	//
/////////////////////////
    int Rows = carrierProps.TotalFluxFromBCRHS.rows();

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from "stiffness term" from electric field
    /////////////////////////////////////////////////////////////////////////////

    // Electric field point values

    // MAKE SURE TO ACCOUNT FOR THE SIGN OF THE CHARGE
    carrierState.ElecPTS =  carrierProps.Sign4Force *
                            VandeMondeMatrices.globalVandeMondeMX * ElecFieldDof;

    // Pointwise multiplication of quadrature weights times electric field
    // quadrature pointvalues
    carrierState.ElecPTSTimesWeights =
        carrierState.ElecPTS.cwiseProduct(weightsDense);

    // convert the dense vector into a sparse one
    carrierState.SparseElecPTSTimesWeights = carrierState.ElecPTSTimesWeights.sparseView();

    ddpBijFlag_type flagI = DG;
    ddpBijFlag_type flagJ = DG;

    ddpMakeGeneralMatrix(problem, grid, carrierState.SparseElecPTSTimesWeights,
                         flagI, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         flagJ, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         carrierState.PE);

    ///////////////////////////////////////////////////////////////////////////////
    // Update Boundary Conditions
    ///////////////////////////////////////////////////////////////////////////////

    // update Dirichlet Boundary condtions
    carrierState.BC_Dir_Input(0) = (carrierProps.BCLeft.OhmicValue +
                                    carrierProps.BCLeft.TestDirichletValue(timeCurrent) );

    carrierState.BC_Dir_Input(1) = (carrierProps.BCRight.OhmicValue +
                                    carrierProps.BCRight.TestDirichletValue(timeCurrent) );

    // Left Robin BC
    carrierState.BC_Robin_Input(0) = (carrierProps.BCLeft.TestRobinValue(timeCurrent)
                                      + carrierProps.BCLeft.RobinValue );

    // Right Robin BC
    carrierState.BC_Robin_Input(1) = (carrierProps.BCRight.TestRobinValue(timeCurrent)
                                      + carrierProps.BCRight.RobinValue );


    /////////////////////////////////////////////////////////////////////////////
    // Assemble Large LDG Matrix
    /////////////////////////////////////////////////////////////////////////////


    ddpSparseMatrix_type
    ATop,
    ATopTranspose,
    ABottom,
    ABottomTranspose,
    LDGABigTranspose;

    double d1 = 1.0;
    double d2 = 0.0;
    ddpConcatenateMatrices(carrierProps.MassQ,
                           carrierProps.B2 - carrierState.PE +  d2 * carrierProps.C,
                           ATop);

    ddpConcatenateMatrices(carrierProps.B1, d1 * carrierProps.C, ABottom);

    ATopTranspose = ATop.transpose();
    ABottomTranspose = ABottom.transpose();

    ddpConcatenateMatrices(ATopTranspose, ABottomTranspose, LDGABigTranspose);

    carrierState.LDG_ABig = LDGABigTranspose.transpose();



    carrierState.LDG_RHS = carrierProps.Dir_RHS * carrierState.BC_Dir_Input
                           + carrierProps.Robin_RHS * carrierState.BC_Robin_Input
                           + carrierState.Generation_RHS;


    // Solve the linear system
    carrierState.Solver_LDG.compute(carrierState.LDG_ABig);
    carrierState.BigSolution = carrierState.Solver_LDG.solve(carrierState.LDG_RHS);



    // get the u and q dof from the big solution
    carrierState.uDof = carrierState.BigSolution.tail(Rows);
    carrierState.qDof = carrierState.BigSolution.head(Rows);

    return 0;

}

int
ddpCarrier_type::
ShurSolve(ddpDenseVector_type const & ElecFieldDof,
          ddpProblemInfo_type const & problem,
          ddpGrid_type const & grid,
          double const & timeCurrent)
{

    int Rows = carrierProps.TotalFluxFromBCRHS.rows();

    /////////////////////////////////////////////////////////////////////////////
    // Compute contribution from "stiffness term" from electric field
    /////////////////////////////////////////////////////////////////////////////

    // Electric field point values

    // MAKE SURE TO ACCOUNT FOR THE SIGN OF THE CHARGE
    carrierState.ElecPTS =  carrierProps.Sign4Force *
                            VandeMondeMatrices.globalVandeMondeMX * ElecFieldDof;


    // Pointwise multiplication of quadrature weights times electric field
    // quadrature pointvalues
    carrierState.ElecPTSTimesWeights =
        carrierState.ElecPTS.cwiseProduct(weightsDense);

    // convert the dense vector into a sparse one
    carrierState.SparseElecPTSTimesWeights = carrierState.ElecPTSTimesWeights.sparseView();

    ddpBijFlag_type flagI = DG;
    ddpBijFlag_type flagJ = DG;

    ddpMakeGeneralMatrix(problem, grid, carrierState.SparseElecPTSTimesWeights,
                         flagI, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         flagJ, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         carrierState.PE);


    ///////////////////////////////////////////////////////////////////////////////
    // Update Boundary Conditions
    ///////////////////////////////////////////////////////////////////////////////

    // update Dirichlet Boundary condtions
    carrierState.BC_Dir_Input(0) = (carrierProps.BCLeft.OhmicValue +
                                    carrierProps.BCLeft.TestDirichletValue(timeCurrent) );

    carrierState.BC_Dir_Input(1) = (carrierProps.BCRight.OhmicValue +
                                    carrierProps.BCRight.TestDirichletValue(timeCurrent) );

    // Left Robin BC
    carrierState.BC_Robin_Input(0) = (carrierProps.BCLeft.TestRobinValue(timeCurrent)
                                      + carrierProps.BCLeft.RobinValue );

    // Right Robin BC
    carrierState.BC_Robin_Input(1) = (carrierProps.BCRight.TestRobinValue(timeCurrent)
                                      + carrierProps.BCRight.RobinValue );


    /////////////////////////////////////////////////////////////////////////////
    // Assemble Large LDG Matrix
    /////////////////////////////////////////////////////////////////////////////

    carrierState.LDG_ABig = ( carrierProps.C
                              - ( carrierProps.B1
                                  * carrierProps.InvMassQ
                                  * (carrierProps.B2 - carrierState.PE) ) );

//	cout << "LDG_ABig = \n" << carrierState.LDG_ABig << endl;

    carrierState.LDG_RHS = 	(carrierProps.TotalFluxFromBCRHS
                             * carrierState.BC_Robin_Input)
                            -
                            (carrierProps.B1
                             * carrierProps.InvMassQ
                             * carrierProps.TotalFluxFromBCRHS
                             * carrierState.BC_Dir_Input);

//	cout << "LDG_RHS = \n" << carrierState.LDG_RHS << endl;



    // solve the linear system
    carrierState.Solver_LDG.compute(carrierState.LDG_ABig);
    carrierState.uDof = carrierState.Solver_LDG.solve(carrierState.LDG_RHS);


//	cout << "uDof = " << carrierState.uDof << endl;
    return 0;

}

/////////////////////////////////////////////////////////////////////////////////
// TIME STEPPING FUNCTIONS.
/////////////////////////////////////////////////////////////////////////////////
int
ddpCarrier_type::
AssembleLDGSystem(ddpProblemInfo_type const & problem,
                  ddpGrid_type const & grid,
                  double const & interface_value,
                  double const & deltaT)
{
    /////////////////////////////////
    // MONOLITHIC BACKWARD EULER	//
/////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    // Assemble Large LDG Matrix
    /////////////////////////////////////////////////////////////////////////////
  //  grvy_timer_begin("Assemble LDG Matrix");

    double d1 = 1.0/grid.DeltaxMin;

    ddpSparseMatrix_type
    ATop,
    ATopTranspose,
    ABottom,
    ABottomTranspose,
    LDGABigTranspose,
    ImplicitInterfaceFlux;


    if(carrierProps.Material == Semiconductor)
        ImplicitInterfaceFlux = carrierProps.RightBoundary_LookInside;
    else if(carrierProps.Material == Electrolyte)
        ImplicitInterfaceFlux = carrierProps.LeftBoundary_LookInside;
    else
    {
        std::cerr << "Material is wrong" << std::endl;
        assert(false);
    }

    ImplicitInterfaceFlux *= interface_value;

    // Now construct the matrices
    // [ (1/scale) A ,						B^{T} + C  ]
    // [  			 	 B,	 					(1/DeltaT) M ]

    ddpConcatenateMatrices( (1.0/carrierProps.CurrentScale) * carrierProps.MassQ,
                            (carrierProps.B2 - carrierState.PE),
                            ATop);

    ddpConcatenateMatrices(carrierProps.B1,
                           (1/deltaT) * carrierProps.MassU + ImplicitInterfaceFlux,
                           ABottom);


    ATopTranspose = ATop.transpose();
    ABottomTranspose =  ABottom.transpose();

    ddpConcatenateMatrices(ATopTranspose, ABottomTranspose, LDGABigTranspose);

    carrierState.LDG_ABig = LDGABigTranspose.transpose();
//    grvy_timer_end("Assemble LDG Matrix");

//    grvy_timer_begin("Factorize LDG Matrix");
// Solve the linear system
    carrierState.Solver_LDG.compute(carrierState.LDG_ABig);
//    grvy_timer_end("Factorize LDG Matrix");


    return 0;
}

int
ddpCarrier_type::
AssembleLDGSystem(ddpProblemInfo_type const & problem,
                  ddpGrid_type const & grid,
                  double const & deltaT)
{
    /////////////////////////////////
    // MONOLITHIC BACKWARD EULER	//
/////////////////////////////////

    /////////////////////////////////////////////////////////////////////////////
    // Assemble Large LDG Matrix
    /////////////////////////////////////////////////////////////////////////////

 //   grvy_timer_begin("Assemble LDG Matrix");

    double d1 = 1.0/grid.DeltaxMin;

    ddpSparseMatrix_type
    ATop,
    ATopTranspose,
    ABottom,
    ABottomTranspose,
    LDGABigTranspose,
    ImplicitInterfaceFlux;


    // Now construct the matrices
    // [ (1/scale) A ,						B^{T} + C  ]
    // [  			 	 B,	 					(1/DeltaT) M ]
    double length = 1.0; //carrierProps.LengthScale;
    double scale = carrierProps.CurrentScale;// * length;

    ddpConcatenateMatrices((1.0/scale) * carrierProps.MassQ,
                           ( carrierProps.B2 - carrierState.PE),  //+ carrierProps.C,
                           ATop);

    ddpConcatenateMatrices(carrierProps.B1,
                           (1.0/deltaT) * carrierProps.MassU,// - d1 * carrierProps.C,
                           ABottom);

    ATopTranspose = ATop.transpose();
    ABottomTranspose =  ABottom.transpose();

    ddpConcatenateMatrices(ATopTranspose, ABottomTranspose, LDGABigTranspose);

    carrierState.LDG_ABig = LDGABigTranspose.transpose();
  //  grvy_timer_end("Assemble LDG Matrix");

 //   grvy_timer_begin("Factorize LDG Matrix");
    // Solve the linear system
    carrierState.Solver_LDG.compute(carrierState.LDG_ABig);
  //  grvy_timer_end("Factorize LDG Matrix");


    return 0;
}


int
ddpCarrier_type::
AssembleRHS(double const & timeCurrent,
            double const & deltaT)
{

  //  grvy_timer_begin("Assemble LDG RHS");

    // Update Boundary Conditions
    ///////////////////////////////////////////////////////////////////////////////

    // update Dirichlet Boundary condtions
    carrierState.BC_Dir_Input(0) = (carrierProps.BCLeft.OhmicValue +
                                    carrierProps.BCLeft.TestDirichletValue(timeCurrent) );

    carrierState.BC_Dir_Input(1) = (carrierProps.BCRight.OhmicValue +
                                    carrierProps.BCRight.TestDirichletValue(timeCurrent) );

    // Left Robin BC
    carrierState.BC_Robin_Input(0) = (carrierProps.BCLeft.TestRobinValue(timeCurrent)
                                      + carrierProps.BCLeft.RobinValue );

    // Right Robin BC
    carrierState.BC_Robin_Input(1) = (carrierProps.BCRight.TestRobinValue(timeCurrent)
                                      + carrierProps.BCRight.RobinValue );


    // Construct the right hand side that is
    // [ dirichlet boundary terms + PE * U^{old}															  ]
    // [ (1/delta t) M * U^{old} + generation + Recombination + interface terms ]

    carrierState.LDG_RHS = // carrierProps.CurrentScale *
        carrierProps.Dir_RHS *
        carrierState.BC_Dir_Input
        + carrierProps.Robin_RHS * carrierState.BC_Robin_Input;

    carrierState.LDG_RHS += (1.0/deltaT) * carrierProps.MassUToSystem * carrierState.uDof;

    carrierState.LDG_RHS += carrierState.Generation_RHS;
    carrierState.LDG_RHS += carrierProps.RecomboToSystem *
                            carrierState.RHSFromRecombination;

    carrierState.LDG_RHS += carrierState.DriftTerm;

 //   grvy_timer_end("Assemble LDG RHS");

    return 0;

}

int
ddpCarrier_type::
UpdateImplicitDriftTerm(ddpDenseVector_type & ElecFieldDof,
                        const ddpProblemInfo_type & problem,
                        const ddpGrid_type & grid)
{
    //////////////////////////////////////////////////////////////////////////////////
    // Compute the drift term from the electric field
    /////////////////////////////////////////////////////////////////////////////////
    //
    // Electric field point values

//    grvy_timer_begin("Assemble Drift Term");

    // MAKE SURE TO ACCOUNT FOR THE SIGN OF THE CHARGE
    carrierState.ElecPTS =  carrierProps.Sign4Force *
                            VandeMondeMatrices.globalVandeMondeMX * ElecFieldDof;
    // Pointwise multiplication of quadrature weights times electric field
    // quadrature pointvalues
    carrierState.ElecPTSTimesWeights =
        carrierState.ElecPTS.cwiseProduct(weightsDense);

    // convert the dense vector into a sparse one
    carrierState.SparseElecPTSTimesWeights = carrierState.ElecPTSTimesWeights.sparseView();

    ddpBijFlag_type flagI = DG;
    ddpBijFlag_type flagJ = DG;

    ddpMakeGeneralMatrix(problem, grid, carrierState.SparseElecPTSTimesWeights,
                         flagI, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         flagJ, Bijections.DGForward, VandeMondeMatrices.globalVandeMondeDG,
                         carrierState.PE);
    //grvy_timer_end("Assemble Drift Term");

    return 0;
}

int
ddpCarrier_type::
UpdateExplicitDriftTerm(ddpDenseVector_type & ElecFieldDof)
{
    //////////////////////////////////////////////////////////////////////////////////
    // Compute the drift term from the electric field
    /////////////////////////////////////////////////////////////////////////////////
    //
    // Electric field point values

   // grvy_timer_begin("Assemble Drift Term");
    // MAKE SURE TO ACCOUNT FOR THE SIGN OF THE CHARGE
    carrierState.ElecPTS =  carrierProps.Sign4Force *
                            VandeMondeMatrices.globalVandeMondeMX * ElecFieldDof;



    // Pointwise multiplication of quadrature weights times electric field
    // quadrature pointvalues
    carrierState.ElecPTSTimesWeights =
        carrierState.ElecPTS.cwiseProduct(weightsDense);

    carrierState.ElecTimesWeightsMatrix =
        carrierState.ElecPTSTimesWeights.asDiagonal();

    // do integration of \int p E u dx as matrix vector multiplication, this really
    // speed things up
    carrierState.DriftTerm =  carrierProps.DriftToSystem *
                              (VandeMondeMatrices.globalVandeMondeDGTransposed *
                               (carrierState.ElecTimesWeightsMatrix *
                                (VandeMondeMatrices.globalVandeMondeDG *
                                 carrierState.uDof) ) );

    //grvy_timer_end("Assemble Drift Term");

    return 0;
}

int
ddpCarrier_type::
BackwardEuler()
{
    //grvy_timer_begin("Solve LDG System");

    int Rows = carrierProps.TotalFluxFromBCRHS.rows();

    // set old solution equal to current solution and solve for new solution
    carrierState.qDof_prev = carrierState.BigSolution.head(Rows);

    // get the u and q dof from the big solution
    carrierState.BigSolution = carrierState.Solver_LDG.solve(carrierState.LDG_RHS);
    carrierState.uDof = carrierState.BigSolution.tail(Rows);
    carrierState.qDof = carrierState.BigSolution.head(Rows);
   // grvy_timer_end("Solve LDG System");

    return 0;
}

bool
ddpCarrier_type::
check_converged(const double & tol)
{
//	grvy_timer_begin("Check if semiconductor converged");
    ddpDenseVector_type q_values = VandeMondeMatrices.globalVandeMondeDG *
                                   carrierState.qDof;

    double left_value = q_values(0);
    for(int i = 1; i < q_values.size(); i++)
    {
        if( fabs(left_value - q_values(i)) > tol)
            return false;
    }

//`	grvy_timer_end("Check if semiconductor converged");
    return true;
}



