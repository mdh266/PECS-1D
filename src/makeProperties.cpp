#include "../include/makeProperties.hpp"

int
ddpMakeCarrierProperies(ddpProblemInfo_type const & problem,
                        ddpGrid_type const & grid,
                        ddpSparseVector_type const & weightsSparse,
                        ddpBijection_type const & Bijections,
                        ddpVandeMondeMatrices_type const & VandeMondeMatrices,
                        ddpDGFluxMatrices_type const & DGFluxMatrices,
                        ddpCarrierProperties_type & CarrierProps)
{


    // The "I" variables denote the test functions
    // While the "J" variable are the trial function.
    ddpBijFlag_type flagI;
    ddpBijFlag_type flagJ;
    ddpSparseMatrix_type vandI;
    ddpSparseMatrix_type vandJ;
    ddpBijForw_type BijI;
    ddpBijForw_type BijJ;

    // DGMassU
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDG;
    BijI = (Bijections.DGForward);

    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ = (Bijections.DGForward);

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI, vandI,
                         flagJ, BijJ, vandJ,
                         CarrierProps.MassU);

    //	cout << "MassU" << CarrierProps.MassU << endl;

    // DGMassQ
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDG;
    BijI = (Bijections.DGForward);

    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ = (Bijections.DGForward);

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI, vandI,
                         flagJ, BijJ, vandJ,
                         CarrierProps.MassQ);

    // InvMassU
    int n = CarrierProps.MassU.rows();
    ddpDenseVector_type x = ddpDenseVector_type::Ones(n);
    ddpDenseVector_type invMassU = CarrierProps.MassU * x;
    for(int i = 0; i < n; i++)
        invMassU(i) = 1.0/invMassU(i);

    CarrierProps.InvMassU = invMassU.asDiagonal();

    n = CarrierProps.MassQ.rows();
    x = ddpDenseVector_type::Ones(n);
    invMassU = CarrierProps.MassQ * x;
    for(int i = 0; i < n; i++)
        invMassU(i) = 1.0/invMassU(i);

    CarrierProps.InvMassQ = invMassU.asDiagonal();


    // StiffUFromQ
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDGPrime;
    BijI = (Bijections.DGForward);

    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ  = (Bijections.DGForward);

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI, vandI,
                         flagJ, BijJ, vandJ,
                         CarrierProps.StiffUFromQ);

    // DGStiffQFromU
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDGPrime;
    BijI  = Bijections.DGForward;


    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ  = Bijections.DGForward;

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI, vandI,
                         flagJ, BijJ, vandJ,
                         CarrierProps.StiffQFromU);

    return 0;
}

int
ddpMakeDiffusiveFluxProperties(
    ddpProblemInfo_type const & problem,
    ddpGrid_type const & grid,
    ddpSparseVector_type const & weightsSparse,
    ddpBijection_type const & Bijections,
    ddpVandeMondeMatrices_type const & VandeMondeMatrices,
    ddpDGFluxMatrices_type const & DGFluxMatrices,
    ddpCarrierProperties_type & CarrierProps)
{
//////////////////////////////////////////////////////////////////////////////
// CREATE FLUX MATRICES
//////////////////////////////////////////////////////////////////////////////

    ddpDirection_type Flag1;

    std::vector<ddpDirection_type> Flag2Vector(grid.NumElementsNoGhost + 1);

    // These will be the matrices of flux values and will multiply the DOFS
    ddpSparseMatrix_type
    LeftEnd_FluxPlus,
    LeftEnd_FluxMinus,
    RightEnd_FluxPlus,
    RightEnd_FluxMinus,
    LeftEnd_Jump,
    LeftEnd_Average,
    RightEnd_Jump,
    RightEnd_Average,
    LeftBoundary_LookInside,
    LeftBoundary_LookOutside,
    RightBoundary_LookInside,
    RightBoundary_LookOutside;

    // SET THE LDG flux Direction
    double Beta = CarrierProps.Beta;
    double Tau = CarrierProps.Tau;

////////////////////////////////////////////////////////////////////
// CREATE FLUXES ON INTERIOR CELLS
///////////////////////////////////////////////////////////////////

    // i is such that it is an inertior element face

    // \lim_{x -> x_{i-1/2}^{+}}
    Flag1 = Minus;
    for(int i = 0; i < grid.NumElementsNoGhost+1; ++i)
    {
        Flag2Vector[i] = Plus;
    }

    ddpMakeInteriorFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   LeftEnd_FluxPlus);

//	 cout << "LeftEnd_FluxPlus = " << LeftEnd_FluxPlus << endl;

    // \lim_{x -> x_{i-1/2}^{-}}
    Flag1 = Minus;
    for(int i = 0; i < grid.NumElementsNoGhost+1; ++i)
    {
        Flag2Vector[i] = Minus;
    }

    ddpMakeInteriorFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   LeftEnd_FluxMinus);


    // \lim_{x -> x_{i+1/2}^{+}}
    Flag1 = Plus;
    for(int i = 0; i < grid.NumElementsNoGhost+1; i++)
    {
        Flag2Vector[i] = Plus;
    }

    //NOTE:  This might not do anything because not picking up
    //boundary.  Was used before.
    Flag2Vector[grid.NumElementsNoGhost] = Minus;

    ddpMakeInteriorFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   RightEnd_FluxPlus);

    // \lim_{x -> x_{i+1/2}^{-}}
    Flag1 = Plus;
    for(int i = 0; i< grid.NumElementsNoGhost+1; ++i)
    {
        Flag2Vector[i] = Minus;
    }


    //NOTE:  This might not do anything because not picking up
    //boundary.  Was used before.
    Flag2Vector[ grid.NumElementsNoGhost ] = Plus;

    ddpMakeInteriorFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   RightEnd_FluxMinus);

/////////////////////////////////////////////////////////////////////////
// CREATE JUMPS AND AVERAGE MATRICES
////////////////////////////////////////////////////////////////////////

    LeftEnd_Jump = (1.0)*LeftEnd_FluxMinus + (-1.0)*LeftEnd_FluxPlus;
    LeftEnd_Average = 0.5*(LeftEnd_FluxMinus + LeftEnd_FluxPlus);
    RightEnd_Jump = (1.0)*RightEnd_FluxMinus + (-1.0)*RightEnd_FluxPlus;
    RightEnd_Average = 0.5*(RightEnd_FluxMinus + RightEnd_FluxPlus);


/////////////////////////////////////////////////////////////////////////////
// Alternating fluxes Choice:
//  				\hat{q} = { q } + 0.5 Beta [[q]]
//  				\hat{u} = { u } - 0.5 Beta [[u]]
//
// Note: Beta chosen in carrier.cpp under create LDG flux matrices
/////////////////////////////////////////////////////////////////////////////

    CarrierProps.FluxRightQFromU = RightEnd_Average + 0.5*Beta * RightEnd_Jump;
    CarrierProps.FluxRightUFromQ = RightEnd_Average - 0.5*Beta * RightEnd_Jump;
    CarrierProps.FluxLeftQFromU  = LeftEnd_Average + 0.5*Beta * LeftEnd_Jump;
    CarrierProps.FluxLeftUFromQ  = LeftEnd_Average - 0.5*Beta * LeftEnd_Jump;


/////////////////////////////////////////////////////////////////////////////
//	Penalization matrix for \hat{q} = \hat{q} - \tau [u] [v]
//
/////////////////////////////////////////////////////////////////////////////

    CarrierProps.C = -Tau * (RightEnd_Jump - LeftEnd_Jump);



// Print the matrices out to check them
    /*
    	cout << "LDG Matrices on interiour cells only:" << endl;
    	cout << "RightUFromQ = \n" << CarrierProps.FluxRightUFromQ << endl;
    	cout << "LeftUFromQ = \n" << CarrierProps.FluxLeftUFromQ << endl;
    	cout << "RightQFromU = \n" << CarrierProps.FluxRightQFromU << endl;
    	cout << "LeftQFromU = \n" << CarrierProps.FluxLeftQFromU << endl;
    	cout << "PenaltyOnU = \n" << CarrierProps.C << endl;
    	// assert(false);
    */

////////////////////////////////////////////////////////////////////////////////
// CREATE THE BOUNDARY FLUX MATRICES
////////////////////////////////////////////////////////////////////////////////
// NOTE:  Matrices which take into account fluxes only on the boundary,
// 				enteries are zero on the interior.
//
// Even though filling in rest of array only last and first matter.
// since calling boundary function.


    // \lim_{x -> x_{0-1/2}^+}
    Flag1 = Minus;
    for(int i = 0; i < grid.NumElementsNoGhost+1 ; ++i)
    {
        Flag2Vector[i] = Plus;
    }

    ddpMakeBoundaryFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   LeftBoundary_LookInside);

    // \lim_{x -> x_{0-1/2}^-}
    Flag1 = Minus;
    for(int i = 0; i < grid.NumElementsNoGhost+1 ; ++i)
    {
        Flag2Vector[i] = Minus;
    }

    //Note: this is zero matrix
    ddpMakeBoundaryFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   LeftBoundary_LookOutside);

    // \lim_{x -> x_{0+1/2}^-}
    Flag1 = Plus;
    for(int i = 0; i < grid.NumElementsNoGhost+1 ; ++i)
    {
        Flag2Vector[i] = Minus;
    }

    ddpMakeBoundaryFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   RightBoundary_LookInside);

    // \lim_{x -> x_{0+1/2}^+}
    // Note: This is zero matrix.
    Flag1 = Plus;
    for(int i = 0; i < grid.NumElementsNoGhost+1 ; ++i)
    {
        Flag2Vector[i] = Plus;
    }

    ddpMakeBoundaryFacesFluxMatrix(grid,
                                   Bijections.DGForward,
                                   DGFluxMatrices,
                                   Flag1,
                                   Flag2Vector,
                                   RightBoundary_LookOutside);




    ///////////////////////////////////////////////////////////////////////
    // TOTAL FLUX MATRICES =
    // 			ADD INTERIOUR FLUX MATRICES AND BOUNDARY FLUX MATRICES
    ///////////////////////////////////////////////////////////////////////


    // LEFT BOUNDARY CONDITION
    if(Dirichlet == CarrierProps.BCLeft.BCType)
    {
        // \hat{q}
        CarrierProps.FluxLeftUFromQ += LeftBoundary_LookInside;
        // \hat{u}
        CarrierProps.FluxLeftQFromU += LeftBoundary_LookOutside;
    }
    else // Neumann or Robin type
    {
        // \hat{q}
        CarrierProps.FluxLeftUFromQ += LeftBoundary_LookOutside;
        // \hat{u}
        CarrierProps.FluxLeftQFromU += LeftBoundary_LookInside;
    }

    // RIGHT BOUNDARY CONDITION
    if(Dirichlet == CarrierProps.BCRight.BCType)
    {
        // \hat{q}
        CarrierProps.FluxRightUFromQ += RightBoundary_LookInside;
        // \hat{u}
        CarrierProps.FluxRightQFromU += RightBoundary_LookOutside;
    }
    else // Neumann or Robin type
    {
        // \hat{q}
        CarrierProps.FluxRightUFromQ += RightBoundary_LookOutside;
        // \hat{u}
        CarrierProps.FluxRightQFromU += RightBoundary_LookInside;
    }


    // Matrices for full implicit robin condition
    CarrierProps.LeftBoundary_LookInside = LeftBoundary_LookInside;
    CarrierProps.RightBoundary_LookInside = RightBoundary_LookInside;




    ////////////////////////////////////////////////////////////////////////
    // BOUNDARY FLUX VECTORS
    ///////////////////////////////////////////////////////////////////////


    // Flux of basis functions on boundaries.  used to multiply the vector
    // that contains the boundary values to give the boundary projection
    // of the boundary functions onto the basis

    ddpSparseMatrix_type
    FluxRightEndPointFromBC,
    FluxLeftEndPointFromBC;

    int elementNumber;

    // DGFluxLeftFromBC
    // This multiplies left Dirichlet or Robin BC
    Flag1 = Minus;  // left side of element
    elementNumber = 0;
    ddpMakeSpecialFluxVector(grid,
                             elementNumber,
                             Bijections.DGForward,
                             DGFluxMatrices,
                             Flag1,
                             FluxLeftEndPointFromBC);

    // DGFluxRightFromBC
    // This multiplies right Dirichlet or Robin BC
    Flag1 = Plus;  //right side of element
    elementNumber = grid.NumElementsNoGhost - 1;
    ddpMakeSpecialFluxVector(grid,
                             elementNumber,
                             Bijections.DGForward,
                             DGFluxMatrices,
                             Flag1,
                             FluxRightEndPointFromBC);

    //concatenate the bondary vectors

    // The minus sign on FluxRigthtEndPoint matrix is due moving the flux vector
    // (which originally had plus sign) to the RHS of the equation.
    // The plus sign is if from moving the flux vector from the
    // LHS to the RHS as well.

    ddpConcatenateMatrices(FluxLeftEndPointFromBC,
                           -FluxRightEndPointFromBC,
                           CarrierProps.TotalFluxFromBCRHS);
//	cout << "TotalFLuxFromBCRHS =\n" << CarrierProps.TotalFluxFromBCRHS << endl;

    ddpConcatenateMatrices(FluxLeftEndPointFromBC,
                           FluxRightEndPointFromBC,
                           CarrierProps.TotalFluxFromBCPenalty);


    CarrierProps.TotalFluxFromBCPenalty *= Tau;
    // LDG Matrices with boundary fluxes

    /*
    	cout << "LDG Matrices with boundary fluxes:" << endl;
    	cout << "RightUFromQ = \n" << CarrierProps.FluxRightUFromQ << endl;
    	cout << "LeftUFromQ = \n" << CarrierProps.FluxLeftUFromQ << endl;
    	cout << "RightQFromU = \n" << CarrierProps.FluxRightQFromU << endl;
    	cout << "LeftQFromU = \n" << CarrierProps.FluxLeftQFromU << endl;
    	cout << "PenaltyOnU = \n" << CarrierProps.C << endl;
    	assert(false);
    */
    return 0;
}



int
ddpMakePoissonProperties(ddpProblemInfo_type const & problem,
                         ddpGrid_type const & grid,
                         ddpSparseVector_type const & weightsSparse,
                         ddpBijection_type const & Bijections,
                         ddpVandeMondeMatrices_type const & VandeMondeMatrices,
                         ddpDGFluxMatrices_type const & DGFluxMatrices,
                         ddpPoissonProperties_type & PoissonProperties)
{

    /// Look at documentation for matrix defintions.


    ddpBijFlag_type flagI;
    ddpBijFlag_type flagJ;
    ddpSparseMatrix_type  vandI;
    ddpSparseMatrix_type  vandJ;
    ddpBijForw_type  BijI;
    ddpBijForw_type  BijJ;



    // MX A00 Matrix
    flagI = Mixed;
    vandI = VandeMondeMatrices.globalVandeMondeMX;
    BijI  = Bijections.MXForward;

    flagJ = Mixed;
    vandJ = VandeMondeMatrices.globalVandeMondeMX;
    BijJ  = Bijections.MXForward;
    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI,vandI,
                         flagJ, BijJ,vandJ,
                         PoissonProperties.A00);

// 	cout << "A00 = \n" << PoissonProperties.A00 << endl;

    // Note:  Was called A01 before
    // MX A10
    flagI = Mixed;
    vandI = VandeMondeMatrices.globalVandeMondeMXPrime;
    BijI  = Bijections.MXForward;

    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ  = Bijections.DGForward;

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI,vandI,
                         flagJ, BijJ,vandJ,
                         PoissonProperties.A10);


//	cout << "A10 =\n " << PoissonProperties.A10 << endl;

    // NOTE: Was called A10 before

    // MX A01
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDG;
    BijI  = (Bijections.DGForward);

    flagJ = Mixed;
    vandJ = VandeMondeMatrices.globalVandeMondeMXPrime;
    BijJ  = Bijections.MXForward;

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI,vandI,
                         flagJ, BijJ,vandJ,
                         PoissonProperties.A01);

//	cout << "A01 =\n " << PoissonProperties.A01 << endl;


    // MXVLeft
    ddpSparseMatrix_type MXVLeftFromBC;

    // This depends on a specific choice of the upsilon
    // basis functions.  This part needs to be changed if a different choice
    // of basis functions is made.

    ddpSparseMatrix_type globalVandeMondeFluxMXTranspose =
        VandeMondeMatrices.globalVandeMondeFluxMX.transpose();

    //  VBCLeft
    MXVLeftFromBC = globalVandeMondeFluxMXTranspose.col(0);

    // MXVRight
    ddpSparseMatrix_type MXVRightFromBC;

    MXVRightFromBC =
        globalVandeMondeFluxMXTranspose.col(grid.NumElementsWithGhost-1);

    // MXVBig
    ddpSparseMatrix_type MXVBig_temp;

    // Concatenate these matrices.

//	 cout << "Entering Concatenate BC" << endl;
    ddpConcatenateMatrices(-MXVLeftFromBC, MXVRightFromBC, MXVBig_temp);
//	 cout << "\tExiting Concatenate BC" << endl;
    PoissonProperties.VBig = MXVBig_temp;
//	 cout << "VBig = \n " << PoissonProperties.VBig	<< endl;

    // MX C Matrix
    flagI = DG;
    vandI = VandeMondeMatrices.globalVandeMondeDG;
    BijI  = (Bijections.DGForward);

    flagJ = DG;
    vandJ = VandeMondeMatrices.globalVandeMondeDG;
    BijJ  = (Bijections.DGForward);

    ddpMakeGeneralMatrix(problem, grid, weightsSparse,
                         flagI, BijI,vandI,
                         flagJ, BijJ,vandJ,
                         PoissonProperties.C);

//		cout << "C = " << PoissonProperties.C << endl;
//		assert(false);

    //MXABig
    ddpSparseMatrix_type MXABig;
    // We make this a block so the extra matrices go out of scope
    // and the memory is freed.
    {
        ddpSparseMatrix_type
        ATop,
        ATopTranspose,
        ABottom,
        ABottomTranspose,
        MXABigTranspose;

        // makes a zero matrix of size A10.rows() x A01.cols() matrix
        ddpSparseMatrix_type zeros(PoissonProperties.A01.rows(),
                                   PoissonProperties.A10.cols() );

        ddpConcatenateMatrices(PoissonProperties.A00,
                               -PoissonProperties.A10, ATop);

        ddpConcatenateMatrices(-PoissonProperties.A01, zeros, ABottom);


        // We have to transpose them and then untranspose them since
        // we can only concatenate matrices column wise
        ATopTranspose = ATop.transpose();
        ABottomTranspose = ABottom.transpose();

        ddpConcatenateMatrices(ATopTranspose, ABottomTranspose, MXABigTranspose);

        PoissonProperties.ABig = MXABigTranspose.transpose();

    }

//	cout << "ABig = " << PoissonProperties.ABig << endl;
    return 0;
}

int
ddpMakePoissonPropertiesWithInterface(ddpProblemInfo_type const & problem,
                                      ddpGrid_type const & grid,
                                      ddpSparseVector_type & weightsSparse,
                                      ddpBijection_type const & Bijections,
                                      ddpVandeMondeMatrices_type const & VandeMondeMatrices,
                                      ddpDGFluxMatrices_type const & DGFluxMatrices,
                                      ddpPoissonProperties_type & PoissonProperties)
{
    // Do the matrices which dont have to change with the interface and
    // then redo A00
    ddpMakePoissonProperties(problem,
                             grid,
                             weightsSparse,
                             Bijections,
                             VandeMondeMatrices,
                             DGFluxMatrices,
                             PoissonProperties);

    ddpBijFlag_type flagI;
    ddpBijFlag_type flagJ;
    ddpSparseMatrix_type  vandI;
    ddpSparseMatrix_type  vandJ;
    ddpBijForw_type  BijI;
    ddpBijForw_type  BijJ;


    // Set up so that
    // A00(I,J) = \int{Support(J)} \epslion^{-}(x) \Psi{I} \Psi_{J} dx
    // by assigning a new sparse vector the multiplication of epsilon(x)
    // times the sparse weights
    int size = weightsSparse.size();
    ddpSparseVector_type EpsilonInvTimesWeights(size);

    for(int i = 0; i < size/2; i++)
    {
        EpsilonInvTimesWeights.insert(i) =  (weightsSparse.coeffRef(i)
                                             / PoissonProperties.SemiconductorPerm );
    }
    for(int i = size/2; i < size; i++)
    {
        EpsilonInvTimesWeights.insert(i) =  (weightsSparse.coeffRef(i)
                                             / PoissonProperties.ElectrolytePerm );
    }

    // MX A00 Matrix
    flagI = Mixed;
    vandI = VandeMondeMatrices.globalVandeMondeMX;
    BijI  = Bijections.MXForward;

    flagJ = Mixed;
    vandJ = VandeMondeMatrices.globalVandeMondeMX;
    BijJ  = Bijections.MXForward;
    ddpMakeGeneralMatrix(problem, grid, EpsilonInvTimesWeights,
                         flagI, BijI,vandI,
                         flagJ, BijJ,vandJ,
                         PoissonProperties.A00);

    //cout << "A00 = \n" << PoissonProperties.A00 << endl;

    // REMAKE THE BIG MATRIX
    // | A00 -A01 |
    // | -A10 0   |

    ddpSparseMatrix_type MXABig;
    // We make this a block so the extra matrices go out of scope
    // and the memory is freed.
    {
        ddpSparseMatrix_type
        ATop,
        ATopTranspose,
        ABottom,
        ABottomTranspose,
        MXABigTranspose;

        // makes a zero matrix of size A10.rows() x A01.cols() matrix
        ddpSparseMatrix_type zeros(PoissonProperties.A01.rows(),
                                   PoissonProperties.A10.cols() );

        ddpConcatenateMatrices(PoissonProperties.A00,
                               -PoissonProperties.A10, ATop);

        ddpConcatenateMatrices(-PoissonProperties.A01, zeros, ABottom);


        // We have to transpose them and then untranspose them since
        // we can only concatenate matrices column wise
        ATopTranspose = ATop.transpose();
        ABottomTranspose = ABottom.transpose();

        ddpConcatenateMatrices(ATopTranspose, ABottomTranspose, MXABigTranspose);

        PoissonProperties.ABig = MXABigTranspose.transpose();
    }
//	 cout << "ABig = " << PoissonProperties.ABig << endl;

    return 0;
}
