#include "../include/matrixSetUp.hpp"

int
ddpMakeGeneralMatrix(ddpProblemInfo_type const & problem,
                     ddpGrid_type const & grid,
                     ddpSparseVector_type const & sparseWeights,
                     ddpBijFlag_type const & flagI,
                     ddpBijForw_type const & ForwardI,
                     ddpSparseMatrix_type const & globalVandeMondeI,
                     ddpBijFlag_type const & flagJ,
                     ddpBijForw_type const & ForwardJ,
                     ddpSparseMatrix_type const & globalVandeMondeJ,
                     ddpSparseMatrix_type & OutputMatrix
                    )
{
// This function returns computes
// M(I,J) = \int_{Support(J)} Psi_{I} Psi_{J} dx

// for GLOBAL DOF I and J.
// We loop over every element to figure out its maximum polynomial
// order and then loop over the neighboring elements polynomials
// and preform quadrature.  The result will be a block diagonal matrix

// NOTE:  This does this for both the DG basis functions and the
// CG basis functions.

// NOTE: Since we have assumed Legrendre polynomials we will have
// diagonal matrix.


    ddpSparseVector_type localWeights = sparseWeights;

    int numElements = grid.NumElementsWithGhost;
    int numTestFunctions = ForwardI.size();  //Test Functions
    int numDof = ForwardJ.size();            //Trial Functions


    ddpSparseMatrix_type MatrixTemp( numTestFunctions, numDof);
    //cout << "will attempt to make MatrixTemp to be of size ("
    //  << numTestFunctions << "x" <<  numDof << ")"  << endl;

    typedef Eigen::Triplet < double > T;
    std::vector<T> triplets;

    // We need to estimate the number of nonzero entries in the matrix.
    // Without thinking too hard, we decide to just use the maximum of
    // DGandMXOrder

    // TODO: Think of a better way to do this.  (if memory becomes an issue).
    int maxOfDGandMXOrder = std::max(problem.MaxOrderDG, problem.MaxOrderMX);

    int estimateOfSize = numTestFunctions * maxOfDGandMXOrder *
                         maxOfDGandMXOrder;

    triplets.reserve(estimateOfSize);

    //We look over all the elements and local orders.
    //Note that we are looping over the ghost elements even for the DG
    //matrices. This is ok since since the MaxOrderDG for the ghost
    // element is set  to something like -1.
    for(int elem1 = 0; elem1 < numElements; ++ elem1)
    {
        // find out if we are doing this for DG or CG basis functions
        // and get the appropriate order
        int localOrderMaxJ;
        switch (flagJ)
        {
        case DG:
            localOrderMaxJ = grid.ElementList[elem1].OrderDG;
            break;
        case Mixed:
            localOrderMaxJ = grid.ElementList[elem1].OrderMX;
            break;
        default:
            cout << "this shouldn't happen on line " << __LINE__ << endl;
            return 1;
        }

        //Now loop over the neighbors of elem1.
        //Note this only does nearest neighbors, if a bigger support is needed
        //will have to change the bounds of the for loop for elem2


        // The plus 2 is required here.
        for(  int elem2 = std::max(0,elem1-1);
                elem2 < std::min(numElements, elem1+2);
                ++elem2)
        {
            // find out if we are doing this for DG or CG basis functions
            // and get the appropriate order
            int localOrderMaxI;
            switch (flagI)
            {
            case DG:
                localOrderMaxI = grid.ElementList[elem2].OrderDG;
                break;
            case Mixed:
                localOrderMaxI = grid.ElementList[elem2].OrderMX;
                break;
            default:
                cout << "this shouldn't happen on line " << __LINE__ << endl;
                return 1;
            }

            // now that we have elem1 and elem2 we want to integrate over
            // each polynomial order of both elem1 and elem2
            for( int order1 = 0; order1 <= localOrderMaxJ; ++order1)
            {
                for( int order2 = 0; order2 <= localOrderMaxI; ++order2)
                {
                    // Get the local DOFS of basis function (elem1, order1)
                    // and (elem2, order2)
                    //
                    std::pair<int, int> elemAndOrderI;
                    std::pair<int, int> elemAndOrderJ;
                    elemAndOrderI.first = elem2;
                    elemAndOrderI.second = order2;
                    elemAndOrderJ.first = elem1;
                    elemAndOrderJ.second = order1;

                    // Find their Gloabl DOF
                    int globalI = ForwardI.find(elemAndOrderI)->second;
                    int globalJ = ForwardJ.find(elemAndOrderJ)->second;

                    // So now we know the global where to find the
                    // values for the functions
                    ddpSparseVector_type temp1 = globalVandeMondeI.col(globalI);
                    ddpSparseVector_type temp2 = globalVandeMondeJ.col(globalJ);

                    //We loop over the nonzero elements of temp1
                    //to preform the quadratuer
                    double partial_sum = 0.0;
                    for(ddpSparseVector_type::InnerIterator it(temp1); it; ++it)
                    {
                        int ind = it.index();
                        double t1val = temp1.coeffRef(ind);
                        double t2val = temp2.coeffRef(ind);
                        double weightval = localWeights.coeffRef(ind);

                        partial_sum += t1val * t2val * weightval;
                    }

                    double valueIJ = partial_sum;


                    // if the result is less than a certain tolerance set it
                    // to zero otherwise store in the list of triplets
                    double absValueIJ = fabs(valueIJ);

                    double Threeeps = 12.0 *
                                      std::numeric_limits<double>::epsilon();

                    if( absValueIJ > Threeeps )
                    {
                        triplets.push_back( T(globalI, globalJ, valueIJ) );
                    }

                } // end for loop over order 2
            } // end for loop over order 1
        } // end for loop over elem2
    } // end for loop over elem1

    // conver the list of triplets into matrix form.
    MatrixTemp.setFromTriplets(triplets.begin(), triplets.end() );
    MatrixTemp.makeCompressed();

    OutputMatrix = MatrixTemp;

    return 0;
}



int
ddpMakeGeneralFluxMatrix(ddpGrid_type const & grid,
                         ddpDenseVector_type const & Vpos,
                         ddpDenseVector_type const & Vneg,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type  const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         ddpSparseMatrix_type & output)
{

// This function will return a matrix of the fluxes of the element
// across the cell interface that will be multiply the DOFs.

// This calculating:
//
// lim_{x -> x_{i \pm 1/2}}  \widehat{ V(x) \Psi_{I}(x)}
//
//     = V(x_{i \pm 1/2}) lim_{x -> x_{i \pm 1/2}^{-}} Psi(I)(x)
//                                        if V(x_{i\pm1/2}) >  0
//
//     = V(x_{i \pm 1/2}) lim_{x -> x_{i \pm 1/2}^{+}} Psi(I)(x)
//                                        if V(x_{i\pm1/2}) <  0
//
//
//
// Flag1=Plus on {i+1/2} fluxes across right boundary
// Flag2=Minus on {i-1/2} fluxes across left boundary
//
// Across all the element nodal points.
//
// Vpos and Vneg are the positive and negative parts of the function V(x)
// which is the point value of the characteristic function (Electric field)
// at ALL the element Nodes.
//
//
// The "upwinding" at each point depends on where the passed in vectors
// Vpos and Vneg are non-zero at those points.  Vpos is the vector only
// the positive values of V and 0 at points where V is negative.

// This helps us by giving us upwinding direction for information to flow
// is to the right.  Vneg is when the upwinding direction for
// information to flow is from the left and contains only values of V
// which are negative.
//
// NOTE:  This takes in a vector of V(x) at ALL the cell nodes.

// NOTE: This is used to make both the diffusive and convective  fluxes.


    int numPoints = grid.NumElementsWithGhost;
    int numElementsNoGhost = grid.NumElementsNoGhost;

    ddpDenseVector_type VposShort, VnegShort;


    const ddpSparseMatrix_type * Mneg;
    const ddpSparseMatrix_type * Mpos;
    const ddpSparseMatrix_type * Mtest;


    // This controls whether or not we take the last or first cell
    // boundary depending on whether or not we are taking flux
    // across left or right side of elements.
    int offset;

    long long int NumTestFunctions = FuncForward.size();


    //  cout << "will attempt to make MatrixTemp to be of size ("
    // << numTestFunctions << "x" <<  numDof << ")"  << endl;


    if(Plus == Flag1)  // we're on the x_{i+1/2}
    {
        // trial functions
        Mneg = &(DGFluxMatrices.plusminus);
        Mpos = &(DGFluxMatrices.plusplus);

        // test functions
        Mtest = &(DGFluxMatrices.plusminus);

        offset = 1;  // we're going to take the
        // last numPoints-1 entries of vector
    }
    else if( Minus == Flag1 ) // were on x_{i-1/2}
    {
        // trial functions
        Mneg = &(DGFluxMatrices.minusminus);
        Mpos = &(DGFluxMatrices.minusplus);

        // test functions
        Mtest = &(DGFluxMatrices.minusplus);

        offset = 0; // we're going to take the
        //  first numPoints-1 entries of vector
    }
    else
    {
        cout << "can't happen on line " << __LINE__ << endl;
    }

    // Subtract out the cell node you dont need
    // and convert to matrix format so easier to use.
    VposShort = Vpos.segment(offset, numPoints - 1);
    VnegShort = Vneg.segment(offset, numPoints - 1);

    ddpDiagonalMatrix_type VposShortDiagMatrix(VposShort.size() );
    ddpDiagonalMatrix_type VnegShortDiagMatrix(VnegShort.size() );

    ddpSparseMatrix_type TEMP1, TEMP2;

    VposShortDiagMatrix = VposShort.asDiagonal();
    VnegShortDiagMatrix = VnegShort.asDiagonal();


    ddpSparseMatrix_type
    temp1,
    temp2;


    // \Psi_I^{\mp} ({x \pm 1/2} ) ( POS(V) * Psi_{J}^{-}(x \pm 1/2)
    temp1 = Mtest->transpose() * VposShortDiagMatrix * (*Mneg);

    // \Psi_I^{\mp} ({x \pm 1/2} ) ( NEG(V) * Psi_{J}^{+}(x \pm 1/2)
    temp2 = Mtest->transpose() * VnegShortDiagMatrix * (*Mpos);

    output = (temp1 + temp2);


    return 0;
}

int
ddpMakeSpecialFluxMatrix(ddpGrid_type const & grid,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         std::vector<ddpDirection_type> const& Flag2Vector,
                         ddpSparseMatrix_type & output )
{
    // This is used in making the diffusive fluxes.  Since there is no general
    // vector for the characteristics we defined a vector of flags which tell
    // at given point what where the information should be flowing form.


    // FLag2Vector[i] = is the direction at which we information should becoming
    // from at node i.  If its Minus we information is flowing to the right and give
    // Vpos[i] = 1 and Vneg[i] = 0.
    // If its plus then information is flowing from the left and fill in
    // Vpos[i] = 0 and Vneg[i] = 1.

    int numElementsWithGhost = grid.NumElementsWithGhost;

    ddpDenseVector_type dummyVPos = ddpDenseVector_type::Zero(numElementsWithGhost);
    ddpDenseVector_type dummyVNeg = ddpDenseVector_type::Zero(numElementsWithGhost);

    // We loop through each element boundary and find out which way the flux is going
    // and fill in the vectors Vpos and Vneg to be passed to ddpMakeGeneralFluxMatrix
    // to do the actual building of the matrix.

    for(int i = 0; i < numElementsWithGhost; ++i)
    {
        if(Plus == Flag2Vector[i])
        {
            dummyVPos(i) = 0.0;
            dummyVNeg(i) = 1.0;
        }
        else if(Minus == Flag2Vector[i])
        {
            dummyVPos(i) = 1.0;
            dummyVNeg(i) = 0.0;
        }
        else
        {
            cout << "error on line " << __LINE__ << endl;
            assert(false);
        }
    }


    // send fake characteristic vectors to make matrix.
    ddpMakeGeneralFluxMatrix( grid,
                              dummyVPos,
                              dummyVNeg,
                              FuncForward,
                              DGFluxMatrices,
                              Flag1,
                              output);
    return 0;
}


// TODO:  Think about merging following to with the functionality of
//        of next two with general flux by setting V = 0 at appropriate
//        places and then send it to ddpMakeGeneralFluxMatrix.
//        NOTE: Might not be worth time.

int
ddpMakeInteriorFacesFluxMatrix(ddpGrid_type const & grid,
                               ddpBijForw_type const & FuncForward,
                               ddpDGFluxMatrices_type  const & DGFluxMatrices,
                               ddpDirection_type const & Flag1,
                               std::vector<ddpDirection_type> const& Flag2Vector,
                               ddpSparseMatrix_type & output)
{

    // This does the same thing as ddpMakeSpecialFluxes but Zeros out V on
    // the boundaries to only make a matrix with interior fluxes and zero for
    // the fluxes on the boundary.
    //

    int numPoints = grid.NumElementsWithGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;

    ddpDenseVector_type Vpos = ddpDenseVector_type::Zero(numElementsWithGhost);
    ddpDenseVector_type Vneg = ddpDenseVector_type::Zero(numElementsWithGhost);

    for(int i = 0; i < numElementsWithGhost; ++i)
    {
        if(Plus == Flag2Vector[i])
        {
            Vpos(i) = 0.0;
            Vneg(i) = 1.0;
        }
        else if(Minus == Flag2Vector[i])
        {
            Vpos(i) = 1.0;
            Vneg(i) = 0.0;
        }
        else
        {
            cout << "error on line " << __LINE__ << endl;
            assert(false);
        }
    }



    ddpDenseVector_type VposShort, VnegShort;

    const ddpSparseMatrix_type * Mneg;
    const ddpSparseMatrix_type * Mpos;
    const ddpSparseMatrix_type * Mtest;

    int offset;

    long long int NumTestFunctions = FuncForward.size();

    if(Plus == Flag1)  // we're on the x_{i+1/2}
    {
        // trial functions
        Mneg = &(DGFluxMatrices.plusminus);
        Mpos = &(DGFluxMatrices.plusplus);

        // test functions
        Mtest = &(DGFluxMatrices.plusminus);

        offset = 1;  // we're going to take the
        // last numPoints-1 entries of vector
    }
    else if( Minus == Flag1 ) // were on x_{i-1/2}
    {
        // trial functions
        Mneg = &(DGFluxMatrices.minusminus);
        Mpos = &(DGFluxMatrices.minusplus);

        // test functions
        Mtest = &(DGFluxMatrices.minusplus);

        offset = 0; // we're going to take the
        //  first numPoints-1 entries of vector
    }
    else
    {
        cout << "can't happen on line " << __LINE__ << endl;
    }

    VposShort = Vpos.segment(offset, numPoints - 1);
    VnegShort = Vneg.segment(offset, numPoints - 1);


    // Set V = 0 on the left boundary point if where taking flux across
    // left side of cells.  Set V = 0 on the right boundary point if were
    // taking flux form right.
    if(Minus == Flag1)
    {
        VposShort(0) = 0.0;
        VnegShort(0) = 0.0;
    }
    else // Plus == Flag1
    {
        VposShort( VposShort.size() - 1) = 0.0;
        VnegShort( VnegShort.size() - 1) = 0.0;
    }


    ddpDiagonalMatrix_type VposShortDiagMatrix(VposShort.size() );
    ddpDiagonalMatrix_type VnegShortDiagMatrix(VnegShort.size() );
    ddpSparseMatrix_type TEMP1, TEMP2;

    VposShortDiagMatrix = VposShort.asDiagonal();
    VnegShortDiagMatrix = VnegShort.asDiagonal();


    ddpSparseMatrix_type
    temp1,
    temp2;


// \Psi_I^{\mp} ({x \pm 1/2} ) ( POS(V) * Psi_{J}^{-}(x \pm 1/2)
    temp1 = Mtest->transpose() * VposShortDiagMatrix * (*Mneg);

// \Psi_I^{\mp} ({x \pm 1/2} ) ( NEG(V) * Psi_{J}^{+}(x \pm 1/2)
    temp2 = Mtest->transpose() * VnegShortDiagMatrix * (*Mpos);

    output = (temp1 + temp2);

    return 0;
}

int
ddpMakeBoundaryFacesFluxMatrix(ddpGrid_type const & grid,
                               ddpBijForw_type const & FuncForward,
                               ddpDGFluxMatrices_type  const & DGFluxMatrices,
                               ddpDirection_type const & Flag1,
                               std::vector<ddpDirection_type> const& Flag2Vector,
                               ddpSparseMatrix_type & output)
{
    // Now this should only compute the fluxes on the boundary elements
    // Flux across boundary cell interface
    int numPoints = grid.NumElementsWithGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;

    ddpDenseVector_type Vpos = ddpDenseVector_type::Zero(numElementsWithGhost);
    ddpDenseVector_type Vneg = ddpDenseVector_type::Zero(numElementsWithGhost);

    for(int i = 0; i < numElementsWithGhost; ++i)
    {
        if(Plus == Flag2Vector[i])
        {
            Vpos(i) = 0.0;
            Vneg(i) = 1.0;
        }
        else if(Minus == Flag2Vector[i])
        {
            Vpos(i) = 1.0;
            Vneg(i) = 0.0;
        }
        else
        {
            cout << "error on line " << __LINE__ << endl;
            assert(false);
        }
    }

    ddpDenseVector_type VposShort, VnegShort;

    const ddpSparseMatrix_type * Mneg;
    const ddpSparseMatrix_type * Mpos;
    const ddpSparseMatrix_type * Mtest;

    int offset;

    long long int NumTestFunctions = FuncForward.size();

    int start, end;


    if(Plus == Flag1)  // we're on the last element
    {
        // trial functions
        Mneg = &(DGFluxMatrices.plusminus);
        Mpos = &(DGFluxMatrices.plusplus);

        // test functions
        Mtest = &(DGFluxMatrices.plusminus);

        offset = 1;  // we're going to take the
        // last numPoints-1 entries of vector

    }
    else if( Minus == Flag1 ) // were on the first element
    {
        // trial functions
        Mneg = &(DGFluxMatrices.minusminus);
        Mpos = &(DGFluxMatrices.minusplus);

        // test functions
        Mtest = &(DGFluxMatrices.minusplus);


        offset = 0; // we're going to take the
        //  first numPoints-1 entries of vector
    }
    else
    {
        cout << "can't happen on line " << __LINE__ << endl;
    }

    VposShort = Vpos.segment(offset, numPoints - 1);
    VnegShort = Vneg.segment(offset, numPoints - 1);



    // set V = 0 on the non boundary Points
    start = 1-offset;
    end = numPoints - (1 + offset);//+offset);

    for(int i = start; i < end; i++)
    {
        VposShort(i) = 0.0;
        VnegShort(i) = 0.0;
    }

    ddpDiagonalMatrix_type VposShortDiagMatrix(VposShort.size() );
    ddpDiagonalMatrix_type VnegShortDiagMatrix(VnegShort.size() );
    ddpSparseMatrix_type TEMP1, TEMP2;

    VposShortDiagMatrix = VposShort.asDiagonal();
    VnegShortDiagMatrix = VnegShort.asDiagonal();


    ddpSparseMatrix_type temp1,
                         temp2;

    // \Psi_I^{\mp} ({x \pm 1/2} ) ( POS(V) * Psi_{J}^{-}(x \pm 1/2)
    temp1 = Mtest->transpose() * VposShortDiagMatrix * (*Mneg);

    // \Psi_I^{\mp} ({x \pm 1/2} ) ( NEG(V) * Psi_{J}^{+}(x \pm 1/2)
    temp2 = Mtest->transpose() * VnegShortDiagMatrix * (*Mpos);

    output = (temp1 + temp2);

    return 0;
}
