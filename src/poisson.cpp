#include "../include/poisson.hpp"

// May need to refer to documentation to see the defintion of
// some the matrices and vectors


int
ddpPoisson_type::
initialize(const ddpGrid_type & grid,
           const ddpProblemInfo_type & problem,
           const ddpCarrierConstants_type & carrierConstants )
{

    // Set the permativities of the materials
    PoissonProps.SemiconductorPerm = 1.0; // problem.semiCondRelativePerm;
    PoissonProps.ElectrolytePerm = 1.0; // problem.electrolyteRelativePerm;

    // For some reason need these to have two values 1/20/2013
    // scaled version of Semiconductors lambda (sqrt of debye length)
    PoissonProps.Lambda = (problem.vacuumPermittivity *
                           problem.semiCondRelativePerm *
                           problem.thermalVoltage) /
                          (problem.characteristicLength *
                           problem.characteristicLength *
                           problem.characteristicDensity *
                           problem.electronCharge);

    // make the weights and points
    ddpMakeWeightsAndPoints(grid,
                            problem,
                            weightsSparse,
                            PTSSparse,
                            PTSDense,
                            weightsDense);

    // Make the bijections for local DOF to global Dof
    ddpMakeAllBijections(grid, Bijections);

    // Make the global VandeMonde Matrice and Flux Matrices
    ddpMakeVandeMondeMatrices(grid,
                              problem,
                              Bijections,
                              VandeMondeMatrices,
                              DGFluxMatrices);

    // Make all the Poisson Matrices
    ddpMakePoissonProperties(problem,
                             grid,
                             weightsSparse,
                             Bijections,
                             VandeMondeMatrices,
                             DGFluxMatrices,
                             PoissonProps);


    // Initialize solver
    PoissonState.solverABig.compute(PoissonProps.ABig);

    // set the couplingStatus i.e. if the RHS has the electrons and holes
    // in it or whether its just he doping profile
    if(On == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 1.0;
    }
    else if(Off == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 0.0;
    }
    else
    {
        cout << "Error: Incorrect Coupling Status." <<endl;
        assert(false);
    }

    // make BCInput a 2x1 array
    PoissonState.BCInput = ddpDenseVector_type(2);

    // get MX DOF
    int RowsInTop = PoissonProps.VBig.rows();

    // get DG DOF
    int RowsInBottom = PoissonProps.C.rows();

    // get DG DOF  (C is square matrix and diagronal)
    int ColsInBottom = PoissonProps.C.cols();

    // Make a matrix with LambdaInverse values along its diagonal
    ddpDenseVector_type LambdaInverseVector(RowsInBottom);

    for(int i = 0; i < RowsInBottom; i++)
    {
        LambdaInverseVector[i] = 1.0/PoissonProps.Lambda;
    }

    ddpDiagonalMatrix_type LambdaInverse = LambdaInverseVector.asDiagonal();

//		cout << "LambdaInverse = \n " << LambdaInverseVector << endl;

    PoissonProps.C = LambdaInverse * PoissonProps.C;

    // make VRHS Matrix that is totalDOF x 2
    PoissonProps.VRHS = ddpSparseMatrix_type(RowsInTop+RowsInBottom,2);

    // make CRHS Matrix that is totalDOF x 2
    PoissonProps.CRHS =
        ddpSparseMatrix_type(RowsInTop + RowsInBottom, ColsInBottom);

    //copy sparse matrix over to larger one for V:
    for(int k = 0; k < PoissonProps.VBig.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.VBig,k);
                it; ++it)
        {
            PoissonProps.VRHS.insert(it.row(), k) = it.value();
        }
    }

//		cout << "VRHS = \n" << PoissonProps.VRHS << endl;

    //copy sparse matrix over to larger one for C:
    for(int k = 0; k < PoissonProps.C.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.C,k); it; ++it)
        {
            PoissonProps.CRHS.insert(RowsInTop + it.row(), k) = it.value();
        }
    }


    // make RHSFromBC vector that will be
    PoissonState.RHSFromBC = ddpDenseVector_type(RowsInTop + RowsInBottom);

    // make RHSFromCandU vector
    PoissonState.RHSFromCandU = ddpDenseVector_type(RowsInTop+RowsInBottom);


    PoissonState.elecDof = ddpDenseVector_type::Zero(RowsInTop);
    PoissonState.potDof = ddpDenseVector_type::Zero(RowsInBottom);

    /*
    	cout << "CRHS = \n" << PoissonProps.CRHS << endl;
    	cout << "A = " << PoissonProps.ABig << endl;

    	cout << "globalVandeMondeMatrixDG \n " << VandeMondeMatrices.globalVandeMondeDG
    			<< endl;
    	cout << "globalVandeMondeMatrixMX \n " << VandeMondeMatrices.globalVandeMondeMX
    			  << endl;
    	cout << "PTS = " << PTSDense << endl;
    */

    return 0;
}

int
ddpPoisson_type::
initializeWithInterface(const ddpGrid_type & left_grid,
                        const ddpGrid_type & right_grid,
                        const ddpProblemInfo_type & problem,
                        const ddpCarrierConstants_type & carrierConstants)
{
    /*
    	cout << "PTS DENSE " << PTSDense << endl;
    	cout << "Weights Dense " << weightsDense << endl;
    	cout << "C = " << PoissonProps.C << endl;
    	cout << "MX VandeMondeMatrix = "
    		   << VandeMondeMatrices.globalVandeMondeMX << endl;
    */
    // Set the permativities of the materials
    PoissonProps.SemiconductorPerm = problem.semiCondRelativePerm;
    PoissonProps.ElectrolytePerm = problem.electrolyteRelativePerm;

    // For some reason need these to have two values 1/20/2013
    // scaled version of Semiconductors lambda (sqrt of debye length)
    PoissonProps.Lambda = (problem.vacuumPermittivity *
                           problem.thermalVoltage) /
                          (problem.characteristicLength *
                           problem.characteristicLength *
                           problem.characteristicDensity *
                           problem.electronCharge);

//	cout << "1/Lambda = " << 1.0/PoissonProps.Lambda << endl;


    // set Poisson's grid to be double that of the transport grid
    ddpGlueGrids(problem, left_grid, right_grid, Poisson_grid);

    // make the weights and points
    ddpMakeWeightsAndPoints(Poisson_grid,
                            problem,
                            weightsSparse,
                            PTSSparse,
                            PTSDense,
                            weightsDense);

    //cout << "PTS DENSE " << PTSDense << endl;
//	cout << "Weights Dense " << weightsDense << endl;

    // Make the bijections for local DOF to global Dof
    ddpMakeAllBijections(Poisson_grid, Bijections);


    // Make the global VandeMonde Matrice and Flux Matrices
    ddpMakeVandeMondeMatrices(Poisson_grid,
                              problem,
                              Bijections,
                              VandeMondeMatrices,
                              DGFluxMatrices);

//	cout << "MX VandeMondeMatrix = "
    //<< VandeMondeMatrices.globalVandeMondeMX << endl;


    // Make all the Poisson Matrices
    ddpMakePoissonPropertiesWithInterface(problem,
                                          Poisson_grid,
                                          weightsSparse,
                                          Bijections,
                                          VandeMondeMatrices,
                                          DGFluxMatrices,
                                          PoissonProps);

    /*
    	cout << "Lambda = " << PoissonProps.Lambda << endl;
    	cout << "Lambda_{S} = "
    			 << PoissonProps.Lambda * PoissonProps.SemiconductorPerm << endl;
    	cout << "Lambda_{E} = "
    			 << PoissonProps.Lambda * PoissonProps.ElectrolytePerm << endl;
    */
//	cout << "C = " << PoissonProps.C << endl;
    /*
    	cout << "globalVandeMondeMatrixDG \n " << VandeMondeMatrices.globalVandeMondeDG
    			<< endl;
    	cout << "globalVandeMondeMatrixMX \n " << VandeMondeMatrices.globalVandeMondeMX
     << endl;
    	assert(false);
    */

    // Initialize solver
    PoissonState.solverABig.compute(PoissonProps.ABig);
//	PoissonState.solverABig.analyzePattern(A);
//	PoissonState.solverABig.factorize(A);

    // set the couplingStatus i.e. if the RHS has the electrons and holes
    // in it or whether its just he doping profile
    if(On == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 1.0;
    }
    else if(Off == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 0.0;
    }
    else
    {
        cout << "Error: Incorrect Coupling Status." <<endl;
        assert(false);
    }

    // make BCInput a 2x1 array
    PoissonState.BCInput = ddpDenseVector_type(2);

    // get MX DOF
    int RowsInTop = PoissonProps.VBig.rows();

    // get DG DOF
    int RowsInBottom = PoissonProps.C.rows();

    // get MX DOF  (C is square matrix and diagronal)
    int ColsInBottom = PoissonProps.C.cols();

    // Make a matrix with LambdaInverse values along its diagonal
    ddpDenseVector_type LambdaInverseVector(RowsInBottom);

    for(int i = 0; i < RowsInBottom; i++)
    {
        LambdaInverseVector[i] = 1.0/PoissonProps.Lambda;
    }

    ddpDiagonalMatrix_type LambdaInverse = LambdaInverseVector.asDiagonal();

//		cout << "LambdaInverse = \n " << LambdaInverseVector << endl;


    PoissonProps.C = LambdaInverse * PoissonProps.C;

    // make VRHS Matrix that is totalDOF x 2
    PoissonProps.VRHS = ddpSparseMatrix_type(RowsInTop+RowsInBottom,2);

    // make CRHS Matrix that is totalDOF x 2
    PoissonProps.CRHS =
        ddpSparseMatrix_type(RowsInTop + RowsInBottom, ColsInBottom);

    //copy sparse matrix over to larger one for V:
    for(int k = 0; k < PoissonProps.VBig.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.VBig,k); it; ++it)
        {
            PoissonProps.VRHS.insert(it.row(), k) = it.value();
        }
    }

//		cout << "VRHS = \n" << PoissonProps.VRHS << endl;

    //copy sparse matrix over to larger one for C:
    for(int k = 0; k < PoissonProps.C.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.C,k); it; ++it)
        {
            PoissonProps.CRHS.insert(RowsInTop + it.row(), k) = it.value();
        }
    }

//		cout << "CRHS = \n" << PoissonProps.CRHS << endl;

    // make RHSFromBC vector that will be
    PoissonState.RHSFromBC = ddpDenseVector_type(RowsInTop + RowsInBottom);

    // make RHSFromCandU vector
    PoissonState.RHSFromCandU = ddpDenseVector_type(RowsInTop+RowsInBottom);

    // Make the matrix to go from semiconductor dof to 2x dof
    int CarrierDof = RowsInBottom / 2;
    int OneSideElecFieldDof = CarrierDof+1;

    // build the matrices which doubles the degrees of freedom into a vector
    ddpSparseMatrix_type temp1(RowsInBottom, CarrierDof);
    ddpSparseMatrix_type temp2(RowsInBottom, CarrierDof);
    ddpSparseMatrix_type temp3(OneSideElecFieldDof, RowsInTop);
    ddpSparseMatrix_type temp4(OneSideElecFieldDof, RowsInTop);

    for(int k = 0; k < CarrierDof; k++)
    {
        temp1.insert(k,k) = 1.0;
        temp2.insert(CarrierDof + k, k ) = 1.0;
    }

    SemicondcutorToTotal = temp1;
    ElectrolyteToTotal = temp2;

    // transpose of these matirces which will return just the
    // semiconductor potential  dofs and the electrolyte
    // potential dofs from the total potetnial dofs
    PotTotalToSemiconductor = temp1.transpose();
    PotTotalToElectrolyte = temp2.transpose();

    // build the matirces which will return just the semiconductor
    // Elec field dofs and the electrolyte field dofs from the
    // total elec field dofs
    for(int k = 0; k < OneSideElecFieldDof; k++)
    {
        temp3.insert(k,k) = 1.0;
        temp4.insert(k,OneSideElecFieldDof-1 + k) = 1.0;
    }	// minus one is to double count the interface pt

    EFTotalToSemiconductor = temp3;
    EFTotalToElectrolyte = temp4;

    PoissonState.elecDof = ddpDenseVector_type::Zero(RowsInTop);
    PoissonState.potDof = ddpDenseVector_type::Zero(RowsInBottom);
    return 0;
}

int
ddpPoisson_type::
remakeProperties(const ddpProblemInfo_type & problem,
                 const ddpGrid_type & grid)
{

    // This remakes all the matrices so that the scaling parameter are 1
    PoissonProps.SemiconductorPerm = 1.0;
    PoissonProps.ElectrolytePerm = 1.0;
    PoissonProps.Lambda = 1.0;

    // Make all the Poisson Matrices
    ddpMakePoissonProperties(problem,
                             grid,
                             weightsSparse,
                             Bijections,
                             VandeMondeMatrices,
                             DGFluxMatrices,
                             PoissonProps);

    // Initialize solver
    PoissonState.solverABig.compute(PoissonProps.ABig);

//		cout << "ABig = " << PoissonProps.ABig << endl;


    // set the couplingStatus i.e. if the RHS has the electrons and holes
    // in it or whether its just he doping profile
    if(On == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 1.0;
    }
    else if(Off == problem.ElecFieldCouplingStatus)
    {
        this->PoissonProps.coupledOrNot = 0.0;
    }
    else
    {
        cout << "Error: Incorrect Coupling Status." <<endl;
        assert(false);
    }

    // make BCInput a 2x1 array
    PoissonState.BCInput = ddpDenseVector_type(2);

    // get MX DOF
    int RowsInTop = PoissonProps.VBig.rows();

    // get DG DOF
    int RowsInBottom = PoissonProps.C.rows();

    // get MX DOF  (C is square matrix and diagronal)
    int ColsInBottom = PoissonProps.C.cols();

    // make VRHS Matrix that is totalDOF x 2
    PoissonProps.VRHS = ddpSparseMatrix_type(RowsInTop+RowsInBottom,2);

    // make CRHS Matrix that is totalDOF x 2
    PoissonProps.CRHS =
        ddpSparseMatrix_type(RowsInTop + RowsInBottom, ColsInBottom);

    //copy sparse matrix over to larger one for V:
    for(int k = 0; k < PoissonProps.VBig.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.VBig,k);
                it; ++it)
        {
            PoissonProps.VRHS.insert(it.row(), k) = it.value();
        }
    }


    //copy sparse matrix over to larger one for C:
    for(int k = 0; k < PoissonProps.C.outerSize(); ++k)
    {
        for(ddpSparseMatrix_type::InnerIterator it(PoissonProps.C,k); it; ++it)
        {
            PoissonProps.CRHS.insert(RowsInTop + it.row(), k) = it.value();
        }
    }

//		cout << "MX_VRHS = " << PoissonProps.VRHS << endl;
//		cout << "MX_CRHS = " << PoissonProps.CRHS << endl;

    // make RHSFromBC vector that will be
    PoissonState.RHSFromBC = ddpDenseVector_type(RowsInTop + RowsInBottom);

    // make RHSFromCandU vector
    PoissonState.elecDof = ddpDenseVector_type::Zero(RowsInTop);
    PoissonState.RHSFromCandU = ddpDenseVector_type(RowsInTop+RowsInBottom);

    return 0;
}


int
ddpPoisson_type::
setDirichletBoundaryConditions(double (* testBCLeft) (const double & t),
                               double (* testBCRight) (const double & t) )
{
    /// This sets the boundary values to the passed in functions pointers
    // which are functions of time.  This is more for testing.

    leftVoltage = 0.0;
    rightVoltage = 0.0;
    testBCDirLeft = testBCLeft;
    testBCDirRight = testBCRight;

    return 0;
}

int
ddpPoisson_type::
setBias(const ddpProblemInfo_type & problem)
{
    // This sets the right bias to the applied bias.  The left side is
    // set to zero or grounded.

    leftVoltage = (problem.BuiltInBias - problem.appliedBias) / problem.thermalVoltage;
    rightVoltage = 0.0;
    testBCDirLeft = &ZeroFunction;
    testBCDirRight = &ZeroFunction;

    return 0;
}


int
ddpPoisson_type::
getElecFieldVals(ddpDenseVector_type & EPTS) const
{
    // This returns an array of the electric fields value at the
    // points in the domain.  It is used to compute the value of
    // DeltaT to make the method stable

    EPTS =  VandeMondeMatrices.globalVandeMondeMX
            * PoissonState.elecDof;

    return 0;
}

int
ddpPoisson_type::
getElectrolyteElecFieldVals(ddpDenseVector_type & EPTS) const
{
    // This returns an array of the electric fields value at the
    // points in the domain.  It is used to compute the value of
    // DeltaT to make the method stable


    // it sets the semiconductor fields to be zero
    EPTS =  VandeMondeMatrices.globalVandeMondeMX
            * PoissonState.elecDof;

//	cout << "Size EPTS = " << EPTS.size() << endl;

    for(int i = 0; i < (EPTS.size()/2); i++)
        EPTS(i) = 0.0;

    return 0;
}

int
ddpPoisson_type::
getSemiconductorElecFieldDOFS(ddpDenseVector_type & EPTS) const
{
    EPTS =	(1.0 / PoissonProps.SemiconductorPerm ) *
            EFTotalToSemiconductor * PoissonState.elecDof;

    return 0;
}

int
ddpPoisson_type::
getElectrolyteElecFieldDOFS(ddpDenseVector_type & EPTS) const
{
    EPTS =  (1.0 / PoissonProps.ElectrolytePerm ) *
            EFTotalToElectrolyte * PoissonState.elecDof;
    return 0;
}

int
ddpPoisson_type::
getSemiconductorPotentialDOFS(ddpDenseVector_type & EPTS) const
{
    EPTS = PotTotalToSemiconductor * PoissonState.potDof;

    return 0;
}

int
ddpPoisson_type::
getElectrolytePotentialDOFS(ddpDenseVector_type & EPTS) const
{
    EPTS = PotTotalToElectrolyte * PoissonState.potDof;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// BIPOLAR MODEL
///////////////////////////////////////////////////////////////////////////////

int
ddpPoisson_type::
solveSystem(ddpCarrier_type const & carrier1,
            ddpCarrier_type const & carrier2,
            ddpDopingProfile_type const & dopingProfile,
            const double & time,
            ddpProblemInfo_type const & problem)
{
    // This will make it so that the RHS of Poisson's equation is
    // C(x) - (n -p)
    // by using Sign4Poisson to change the constant in front of the carrier's
    // dof.  Note, the order of the carrier in the solveSystem does NOT
    // matter.

    PoissonState.carrierDof =
        carrier1.carrierProps.Sign4Poisson *
        carrier1.carrierState.uDof
        +
        carrier2.carrierProps.Sign4Poisson *
        carrier2.carrierState.uDof;

    // Solve thes system for the new DOF values of the electric field
    // and potential
    solveSystem(PoissonState.carrierDof,dopingProfile, time, problem);

    return 0;

}

///////////////////////////////////////////////////////////////////////////////
// UNIPOLAR MODEL
///////////////////////////////////////////////////////////////////////////////

int
ddpPoisson_type::
solveSystem(ddpCarrier_type const & carrier1,
            ddpDopingProfile_type const & dopingProfile,
            const double & time,
            ddpProblemInfo_type const & problem )
{

    // This will make it so that the RHS of Poisson's equation is
    // C(x) - n     OR    C(x) - (-p) = C(x) + p\
    // by using Sign4Poisson to change the constant infront of the carriers
    // dof

    PoissonState.carrierDof = carrier1.carrierProps.Sign4Poisson
                              * carrier1.carrierState.uDof;

    // Solve thes system for the new DOF values of the electric field
    // and potential
    solveSystem(PoissonState.carrierDof,dopingProfile, time, problem);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// BOTH MODELS (The one each calls)
///////////////////////////////////////////////////////////////////////////////

int
ddpPoisson_type::
solveSystem(ddpDenseVector_type const & carrier1Dof,
            ddpDopingProfile_type const & dopingProfile,
            const double & time,
            ddpProblemInfo_type const & problem )
{

// Both the unipolar and bipolar models will call this function
// which will solve for the current time values for the DOF
// for the electric field and potential

//    grvy_timer_begin("Build Poisson RHS");

    int numDofPotential = PoissonProps.C.rows();
    int numDofElec = PoissonProps.ABig.rows() - numDofPotential;

    // Update the boundary conditions (if they are time dependent)
    // and store the minto an array to be multiplied by a matrix
    // to make the RHS

    // Note in testing or setting the running, one of the terms in the
    // sum on each line will be zero
    PoissonState.BCInput(0) = (leftVoltage + testBCDirLeft(time) );
    PoissonState.BCInput(1) = (rightVoltage + testBCDirRight(time) );

    //Update the RHSFromBC vector to include new values of boundary conditons
    //if they are time dependent.
    PoissonState.RHSFromBC =  -PoissonProps.VRHS * PoissonState.BCInput;


    //Update the RHSFRomCandU vector taking into accound new values for the
    //carriers if the system is coupled.  The value couledOrNot takes into
    //account if the system is coupled
    PoissonState.RHSFromCandU = PoissonProps.CRHS *
                                (dopingProfile.Dof - PoissonProps.coupledOrNot* carrier1Dof);

    // Get the WHOLE RHS side from both the carriers and the boundary conditions
    PoissonState.RHSTotal = PoissonState.RHSFromBC - PoissonState.RHSFromCandU;

//    grvy_timer_end("Build Poisson RHS");

//    grvy_timer_begin("Solve Poisson System");

    // This line solvers Abig * PoissonState.Soln = PoissonStateRHSTotal
    // for PoissonState.Soln

    PoissonState.Soln = PoissonState.solverABig.solve(PoissonState.RHSTotal);

    // test to make sure solver worked.
    if(this->PoissonState.solverABig.info() != Eigen::Success)
    {
        cout << "failed in the solver step" << endl;
        assert(false);
    }


    // Extract the electric field dof from top of solution vector
    PoissonState.elecDof = PoissonState.Soln.head(numDofElec);

    // Extract the potential dof from the bottom of the solution vector
    PoissonState.potDof  = PoissonState.Soln.tail(numDofPotential);

  //  grvy_timer_end("Solve Poisson System");

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// SEMICONDUCTOR ELECTROLYTE MODEL
///////////////////////////////////////////////////////////////////////////////

int
ddpPoisson_type::
solveSystem(ddpCarrier_type const & electrons,
            ddpCarrier_type const & holes,
            ddpCarrier_type const & reductants,
            ddpCarrier_type const & oxidants,
            ddpDopingProfile_type const & dopingProfile,
            const double & time,
            ddpProblemInfo_type const & problem)
{
   // grvy_timer_begin("Assemble Poisson RHS");

    // This will make it so that the RHS of Poisson's equation is
    // [ C(x) - (\rho_{n} - \rho_{p} ) ]
    // [ 0 - (\rho_{r} - \rho_{o} ) ]
    //
    // by using Sign4Poisson to change the constant in front of the carrier's
    // dof.  Note, the order of the carrier in the solveSystem does NOT
    // matter.
    ddpDenseVector_type SemiconductorDof = SemicondcutorToTotal * (
            electrons.carrierProps.Sign4Poisson *
            electrons.carrierState.uDof
            +
            holes.carrierProps.Sign4Poisson *
            holes.carrierState.uDof );

    ddpDenseVector_type ElectrolyteDof = ElectrolyteToTotal * (
            reductants.carrierProps.Sign4Poisson *
            reductants.carrierState.uDof
            +
            oxidants.carrierProps.Sign4Poisson *
            oxidants.carrierState.uDof );
    /*
    	cout << "SDOF \n" <<  SemiconductorDof << endl;
    	cout << "EDOF \n" << ElectrolyteDof << endl;

    	assert(false);
    */
    PoissonState.carrierDof = SemiconductorDof + ElectrolyteDof;

// Both the unipolar and bipolar models will call this function
// which will solve for the current time values for the DOF
// for the electric field and potential


    int numDofPotential = PoissonProps.C.rows();
    int numDofElec = PoissonProps.ABig.rows() - numDofPotential;

    // Update the boundary conditions (if they are time dependent)
    // and store the minto an array to be multiplied by a matrix
    // to make the RHS

    // Note in testing or setting the running, one of the terms in the
    // sum on each line will be zero
    PoissonState.BCInput(0) = (leftVoltage + testBCDirLeft(time) );
    PoissonState.BCInput(1) = (rightVoltage + testBCDirRight(time) );

    //Update the RHSFromBC vector to include new values of boundary conditons
    //if they are time dependent.
    PoissonState.RHSFromBC =  -PoissonProps.VRHS * PoissonState.BCInput;

    //Update the RHSFRomCandU vector taking into accound new values for the
    //carriers if the system is coupled.  The value couledOrNot takes into
    //account if the system is coupled
    PoissonState.RHSFromCandU = PoissonProps.CRHS *
                                ( (SemicondcutorToTotal * dopingProfile.Dof)
                                  - (PoissonProps.coupledOrNot * PoissonState.carrierDof) );
    // NOTE:  carrierDof is already in the form of having twice DOF from
    // the beginning of this function, so no need to multiply it by the
    // ___ToTotal matrix

    //cout << "RHS" << PoissonState.RHSFromCandU << endl;

    // Get the WHOLE RHS side from both the carriers and the boundary conditions
    PoissonState.RHSTotal = PoissonState.RHSFromBC - PoissonState.RHSFromCandU;

   // grvy_timer_end("Assemble Poisson RHS");

   // grvy_timer_begin("Solve Poisson System");

    // This line solvers Abig * PoissonState.Soln = PoissonStateRHSTotal
    // for PoissonState.Soln

    PoissonState.Soln = PoissonState.solverABig.solve(PoissonState.RHSTotal);

    // test to make sure solver worked.
    if(this->PoissonState.solverABig.info() != Eigen::Success)
    {
        cout << "failed in the solver step" << endl;
        assert(false);
    }

    // Extract the electric field dof from top of solution vector
    PoissonState.elecDof = PoissonState.Soln.head(numDofElec);

    // Extract the potential dof from the bottom of the solution vector
    PoissonState.potDof  = PoissonState.Soln.tail(numDofPotential);

    //grvy_timer_end("Solve Poisson System");

    return 0;

}


