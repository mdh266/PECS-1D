#include "../include/vectorSetUp.hpp"


// NOTE: Not really a point in having a general flux vector
// and special flux vector maker since the general one is only called
// by the special one at one point in the code.
//
// These are used to pick up the boundary values through q


int
ddpMakeGeneralFluxVector(ddpGrid_type const & grid,
                         ddpDenseVector_type const & VPosOrNeg,
                         int const & ElementNumber,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type  const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         ddpSparseMatrix_type & output)
{

    // This function is to make a vector of the flux of the DG basis funtions
    // on the boundary elements.  This is so we can multiply them by the Dirichlet
    // boundary values so q can pick them up.  It takes  in vector VposOrNeg like
    // to mulitply the vector, but really its just a vector of ones so it doesnt
    // matter


    int numElementsNoGhost = grid.NumElementsNoGhost;

    const ddpSparseMatrix_type * Mneg;
    const ddpSparseMatrix_type * Mpos;
    const ddpSparseMatrix_type * Mtest;

    int offset;  // helps us get the point value of the left
    // or right side of the elements

    long long int NumTestFunctions = FuncForward.size();

    // we build it as a sparse matrix of just one column.
    // so we set up it up as a list of triples (row, col, value)
    // and then convert it to a sparse matrix after we are done.

    int one = 1;

    ddpSparseMatrix_type VectorTemp( NumTestFunctions, one);

    typedef
    Eigen::Triplet < double > T;

    std::vector<T> triplets;
    long long int estimateofSize = 20;
    triplets.reserve(estimateofSize);


    // If flag1 == Plus we are taking the flux across the
    // right side of the boundary Element

    // If flag1 == Minus, we are taking the flux across the
    // left side of the boundary element
    if(Plus == Flag1)
    {
        Mneg = &(DGFluxMatrices.plusminus);
        Mpos = &(DGFluxMatrices.plusplus);
        Mtest = &(DGFluxMatrices.plusminus);
        offset = 1;
    }
    else if( Minus == Flag1 )
    {
        Mneg = &(DGFluxMatrices.minusminus);
        Mpos = &(DGFluxMatrices.minusplus);
        Mtest = &(DGFluxMatrices.minusplus);
        offset = 0;
    }
    else
    {
        cout << "can't happen on line " << __LINE__ << endl;
        assert(false);
    }

    // making a column vector so only nonzero values will be in column 0
    int GlobalJ = 0;

    // only fill in for the element that straddles the boundary of
    // the boundary condition we are trying to pick up
    for(int elem1 = ElementNumber; elem1 < ElementNumber+1; ++elem1)
    {
        // loop over all the polynomial orders of the basis function
        // of the element you are in.
        int maxLocalOrderDG1 = (grid.ElementList[elem1]).OrderDG;
        for(int order1 = 0; order1 <= maxLocalOrderDG1; ++ order1)
        {
            //Compute GlobalI
            std::pair<int, int> elemAndOrderI;

            elemAndOrderI.first = elem1;
            elemAndOrderI.second = order1;

            // find the global DOF index of (elemNumber, Order)
            int GlobalI = FuncForward.find(elemAndOrderI)->second;

            double ValueI;
            double termPosOrNeg;

            double termTest;

            // this gets point value of the left or right side of the
            // element you are in.
            int elem1PlusOffset = elem1 + offset;

            // Get the passed in function is value at that end point.
            // ( vector that was made and passed in) It should be 1... I think...
            termPosOrNeg = VPosOrNeg.coeff(elem1PlusOffset);

            // value of the test term at the element for DOF
            termTest = Mtest->coeff(elem1, GlobalI);

            // value of the function at the end point times the flux of the
            // basis functions across that cell end point.
            ValueI = termPosOrNeg * termTest;

            double absValueI = fabs(ValueI);

            // only store if bigger than machine epsilon
            double Threeeps  =
                12.0 * std::numeric_limits<double>::epsilon();

            if( absValueI > Threeeps )
            {
                triplets.push_back( T(GlobalI,GlobalJ, ValueI) );
            }

        } // end for order1
    } // end for elem

    VectorTemp.setFromTriplets(triplets.begin(), triplets.end() );
    VectorTemp.makeCompressed();
    output = VectorTemp;

    return 0;
}

int
ddpMakeSpecialFluxVector(ddpGrid_type const & grid,
                         int const & ElementNumber,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         ddpSparseMatrix_type & output)
{
    // create a fake electric field like we did for ddpMakeSpecialFluxMatrix
    // that is we pass to ddpMakeGeneralFluxVector to make the actual flux
    // vector

    int numElementsWithGhost = grid.NumElementsWithGhost;
    ddpDenseVector_type dummyV(numElementsWithGhost);

    for(int i = 0; i < numElementsWithGhost; ++i)
    {
        dummyV(i) = 1.0;
    }

    ddpMakeGeneralFluxVector(grid,
                             dummyV,
                             ElementNumber,
                             FuncForward,
                             DGFluxMatrices,
                             Flag1,
                             output);

    return 0;
}


