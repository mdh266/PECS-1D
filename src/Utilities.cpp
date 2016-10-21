#include "../include/includes.hpp"
#include "../include/main.hpp"


int
ddpConcatenateMatrices(ddpSparseMatrix_type const & A,
                       ddpSparseMatrix_type const & B,
                       ddpSparseMatrix_type & C)
{
    // Concatnates to columnwise matrices A, B. to a new matrix C = [A B]

    assert(A.rows() == B.rows());
    C.resize(A.rows(), A.cols() + B.cols());
    C.middleCols(0,A.cols()) = A;
    C.middleCols(A.cols(),B.cols()) = B;
    return 0;
}



int
ddpPositivePart(ddpDenseVector_type const & input, ddpDenseVector_type & output)
{

    // takes input vector input and returns a new vector ouput
    // which has replaced inputs negative values with 0
    ddpDenseVector_type tempVector(input.size());

    for(int i = 0; i < input.size(); ++i)
    {
        tempVector(i) = std::max(0.0, input.coeffRef(i) );
    }

    output = tempVector;
    return 0;

}

int
ddpNegativePart(ddpDenseVector_type const & input, ddpDenseVector_type & output)
{
    // takes input vector input and returns a new vector ouput
    // which has replaced inputs positive values with 0

    ddpDenseVector_type tempVector(input.size());

    for(int i = 0; i < input.size(); ++i)
    {
        tempVector(i) = std::min(0.0, input.coeffRef(i) );
    }

    output = tempVector;
    return 0;
}

int
ddpProjectFunction(
    double (ddpGenerationFunction_type::*function)(double const &),
    ddpGenerationFunction_type & GenFunObject,
    ddpGrid_type const & grid,
    ddpBijForw_type const & FuncForward,
    ddpSparseMatrix_type const & globalVandeMondeJ,
    ddpBijForw_type const & PTForward,
    ddpBijFlag_type const & BijFlag,
    ddpSparseVector_type const & sparseWeights,
    ddpSparseMatrix_type const & MassMatrix,
    ddpDenseVector_type & output)
{

// Projects a member function of the ddpGenerationFunction_type
// which will be based as  function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed
// using the BijFlag.  Notice we need the MassMatrix for this,
// because we have not assumed that the basis functions are orthogonal

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);
    ddpDenseVector_type tempOutput(numTestFunctions);


    ddpDenseVector_type FunctionPointVals( PTForward.size() );

    // Build a vector that contains that the function to be projected,
    // with its value at each sample point.


    // loop over the number of elements
    for(int elem = 0; elem < lastElement; ++elem)
    {
        int localNumGaussLegendrePoints =
            grid.ElementList[elem].GaussLegendreNumPoints;

        // loop over the number of points in each element and sample function there
        for(int pointNum = 0; pointNum < localNumGaussLegendrePoints; ++pointNum)
        {
            double xi = grid.ElementList[elem].GaussLegendrePoints[pointNum];

            std::pair<int, int> elemAndPointJ;
            elemAndPointJ.first = elem;
            elemAndPointJ.second = pointNum;

            // get global Dof corresponding to point and use it to fill in
            int globalJ = PTForward.find(elemAndPointJ)->second;
            FunctionPointVals(globalJ) = (GenFunObject.*function)(xi);

        }
    }

    ddpProjectFunction(FunctionPointVals,
                       grid,
                       FuncForward,
                       globalVandeMondeJ,
                       PTForward,
                       BijFlag,
                       sparseWeights,
                       MassMatrix,
                       output);



    return 0;
}

int
ddpProjectFunction( double (*function)(double const &),
                    ddpGrid_type const & grid,
                    ddpBijForw_type const & FuncForward,
                    ddpSparseMatrix_type const & globalVandeMondeJ,
                    ddpBijForw_type const & PTForward,
                    ddpBijFlag_type const & BijFlag,
                    ddpSparseVector_type const & sparseWeights,
                    ddpSparseMatrix_type const & MassMatrix,
                    ddpDenseVector_type & output)
{
// Projects a function name function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed
// using the BijFlag.  Notice we need the MassMatrix for this,
// because we have not assumed that the basis functions are orthogonal

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);
    ddpDenseVector_type tempOutput(numTestFunctions);


    ddpDenseVector_type FunctionPointVals( PTForward.size() );

    // Build a vector that contains that the function to be projected,
    // with its value at each sample point.


    // loop over the number of elements
    for(int elem = 0; elem < lastElement; ++elem)
    {
        int localNumGaussLegendrePoints =
            grid.ElementList[elem].GaussLegendreNumPoints;

        // loop over the number of points in each element and sample function there
        for(int pointNum = 0; pointNum < localNumGaussLegendrePoints; ++pointNum)
        {
            double xi = grid.ElementList[elem].GaussLegendrePoints[pointNum];

            std::pair<int, int> elemAndPointJ;
            elemAndPointJ.first = elem;
            elemAndPointJ.second = pointNum;

            // get global Dof corresponding to point and use it to fill in
            int globalJ = PTForward.find(elemAndPointJ)->second;
            FunctionPointVals(globalJ) = function(xi);

        }
    }


    ddpProjectFunction(FunctionPointVals,
                       grid,
                       FuncForward,
                       globalVandeMondeJ,
                       PTForward,
                       BijFlag,
                       sparseWeights,
                       MassMatrix,
                       output);

    return 0;
}

int
ddpProjectFunction( ddpDenseVector_type & FunctionPointVals,
                    ddpGrid_type const & grid,
                    ddpBijForw_type & FuncForward,
                    ddpSparseMatrix_type & globalVandeMondeJ,
                    ddpBijForw_type & PTForward,
                    ddpBijFlag_type & BijFlag,
                    ddpSparseVector_type & sparseWeights,
                    ddpSparseMatrix_type & MassMatrix,
                    ddpDenseVector_type & output)
{
// Projects a function name function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed
// using the BijFlag.  Notice we need the MassMatrix for this,
// because we have not assumed that the basis functions are orthogonal

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);
    ddpDenseVector_type tempOutput(numTestFunctions);



    //Build the vector of the functions when projected onto the
    // Chosen basis (DG or Mixed CG) from the functions point values

    int maxOrder;

    // loop over each element and find out what kind of basis function
    for(int elem = 0; elem < lastElement; ++elem)
    {
        //cout << "elem = " << elem << endl;

        if(DG == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderDG;
        }
        else if(Mixed == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderMX;
        }
        else
        {
            assert(false);
            return 1;
        }

        // loop over the order of the polynomials of the basis function's
        // element
        //
        for(int order = 0; order <= maxOrder; ++order)
        {
            // get the global DOF from the local DOF
            //cout << "\torder = " << order << endl;
            std::pair<int, int> elemAndOrderJ;
            elemAndOrderJ.first = elem;
            elemAndOrderJ.second = order;
            int globalJ =  FuncForward.find(elemAndOrderJ)->second;

            ddpSparseVector_type temp = globalVandeMondeJ.col(globalJ);

            //We loop over the nonzero elements of temp1 and project function
            //onto the shape function/basis function by integrating that function
            //against the basis function
            double partial_sum = 0.0;
            for(ddpSparseVector_type::InnerIterator it(temp); it; ++it)
            {
                int ind = it.index();
                double tempval = temp.coeffRef(ind);
                double weightval = localWeights.coeffRef(ind);
                double FunctionVal = FunctionPointVals.coeffRef(ind);
                partial_sum += tempval * FunctionVal * weightval;
            }
            double valueJ = partial_sum;

            double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

            // only store that value in the vector if bigger than a certain
            // tolerance
            if(fabs(valueJ) > Threeeps)
            {
                RHSVector(globalJ) = valueJ;
            }
            else
            {
                RHSVector(globalJ) = 0.0;
            }
        }
    }


    //perform a linear solve to find the function's degrees of freedom
    // This is because the basis functions are not necessarily ORTHONORMAL
    ddpSolver_type solver0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver1;


    solver0.compute(MassMatrix);
    if(solver0.info() != Eigen::Success)
    {
        cout << "solver0 compute step didn't work" << endl;
    }

    tempOutput = solver0.solve(RHSVector);

    if(solver0.info() == Eigen::Success)
    {
        //cout << "Solver0 did work" << endl;
        output = tempOutput;
        return 0;
    }
    else
    {
        // output the matrix and RHS if neither worked
        tempOutput = solver1.compute(MassMatrix).solve(RHSVector);
        if(solver1.info() == Eigen::Success)
        {
            //cout << "Solver1 did work" << endl;
            output = tempOutput;
            return 0;
        }
        else
        {
            cout << "MassMatrix = \n" << MassMatrix  << endl;
            cout << "RHSVector = \n" << RHSVector << endl;
            cout << "Neither solver worked" <<endl;
            assert(false);
            return 1;
        }
    }

    return 0;
}

int
ddpProjectFunction( ddpDenseVector_type const & FunctionPointVals,
                    ddpGrid_type const & grid,
                    ddpBijForw_type const & FuncForward,
                    ddpSparseMatrix_type const & globalVandeMondeJ,
                    ddpBijForw_type const & PTForward,
                    ddpBijFlag_type const & BijFlag,
                    ddpSparseVector_type const & sparseWeights,
                    ddpSparseMatrix_type const & MassMatrix,
                    ddpDenseVector_type & output)
{
// Projects a function name function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed
// using the BijFlag.  Notice we need the MassMatrix for this,
// because we have not assumed that the basis functions are orthogonal

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);
    ddpDenseVector_type tempOutput(numTestFunctions);



    //Build the vector of the functions when projected onto the
    // Chosen basis (DG or Mixed CG) from the functions point values

    int maxOrder;

    // loop over each element and find out what kind of basis function
    for(int elem = 0; elem < lastElement; ++elem)
    {
        //cout << "elem = " << elem << endl;

        if(DG == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderDG;
        }
        else if(Mixed == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderMX;
        }
        else
        {
            assert(false);
            return 1;
        }

        // loop over the order of the polynomials of the basis function's
        // element
        //
        for(int order = 0; order <= maxOrder; ++order)
        {
            // get the global DOF from the local DOF
            //cout << "\torder = " << order << endl;
            std::pair<int, int> elemAndOrderJ;
            elemAndOrderJ.first = elem;
            elemAndOrderJ.second = order;
            int globalJ =  FuncForward.find(elemAndOrderJ)->second;

            ddpSparseVector_type temp = globalVandeMondeJ.col(globalJ);

            //We loop over the nonzero elements of temp1 and project function
            //onto the shape function/basis function by integrating that function
            //against the basis function
            double partial_sum = 0.0;
            for(ddpSparseVector_type::InnerIterator it(temp); it; ++it)
            {
                int ind = it.index();
                double tempval = temp.coeffRef(ind);
                double weightval = localWeights.coeffRef(ind);
                double FunctionVal = FunctionPointVals.coeffRef(ind);
                partial_sum += tempval * FunctionVal * weightval;
            }
            double valueJ = partial_sum;

            double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

            // only store that value in the vector if bigger than a certain
            // tolerance
            if(fabs(valueJ) > Threeeps)
            {
                RHSVector(globalJ) = valueJ;
            }
            else
            {
                RHSVector(globalJ) = 0.0;
            }
        }
    }


    //perform a linear solve to find the function's degrees of freedom
    // This is because the basis functions are not necessarily ORTHONORMAL
    ddpSolver_type solver0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver1;


    solver0.compute(MassMatrix);
    if(solver0.info() != Eigen::Success)
    {
        cout << "solver0 compute step didn't work" << endl;
    }

    tempOutput = solver0.solve(RHSVector);

    if(solver0.info() == Eigen::Success)
    {
        //cout << "Solver0 did work" << endl;
        output = tempOutput;
        return 0;
    }
    else
    {
        // output the matrix and RHS if neither worked
        tempOutput = solver1.compute(MassMatrix).solve(RHSVector);
        if(solver1.info() == Eigen::Success)
        {
            //cout << "Solver1 did work" << endl;
            output = tempOutput;
            return 0;
        }
        else
        {
            cout << "MassMatrix = \n" << MassMatrix  << endl;
            cout << "RHSVector = \n" << RHSVector << endl;
            cout << "Neither solver worked" <<endl;
            assert(false);
            return 1;
        }
    }

    return 0;
}

// Needs to be Utilities.cpp so that it can be used in carrier.cpp
double ZeroFunction(const double & t)
{
    return 0.0;
}


int
ddpProjectRHSFunction( double (*function)(double const &),
                       ddpGrid_type const & grid,
                       ddpBijForw_type const & FuncForward,
                       ddpSparseMatrix_type const & globalVandeMondeJ,
                       ddpBijForw_type const & PTForward,
                       ddpBijFlag_type const & BijFlag,
                       ddpSparseVector_type const & sparseWeights,
                       ddpDenseVector_type & output)
{
// Projects a function name function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed
// using the BijFlag.

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);
    ddpDenseVector_type tempOutput(numTestFunctions);

    ddpDenseVector_type FunctionPointVals( PTForward.size() );

    // Build a vector that contains that the function to be projected,
    // with its value at each sample point.


    // loop over the number of elements
    for(int elem = 0; elem < lastElement; ++elem)
    {
        int localNumGaussLegendrePoints =
            grid.ElementList[elem].GaussLegendreNumPoints;

        // loop over the number of points in each element and sample function there
        for(int pointNum = 0; pointNum < localNumGaussLegendrePoints; ++pointNum)
        {
            double xi = grid.ElementList[elem].GaussLegendrePoints[pointNum];

            std::pair<int, int> elemAndPointJ;
            elemAndPointJ.first = elem;
            elemAndPointJ.second = pointNum;

            // get global Dof corresponding to point and use it to fill in
            int globalJ = PTForward.find(elemAndPointJ)->second;
            FunctionPointVals(globalJ) = function(xi);

        }
    }


    ddpProjectRHSFunction(FunctionPointVals,
                          grid,
                          FuncForward,
                          globalVandeMondeJ,
                          PTForward,
                          BijFlag,
                          sparseWeights,
                          output);

} // end projectRHS


int
ddpProjectRHSFunction( ddpDenseVector_type const & FunctionPointVals,
                       ddpGrid_type const & grid,
                       ddpBijForw_type const & FuncForward,
                       ddpSparseMatrix_type const & globalVandeMondeJ,
                       ddpBijForw_type const & PTForward,
                       ddpBijFlag_type const & BijFlag,
                       ddpSparseVector_type const & sparseWeights,
                       ddpDenseVector_type & output)
{
// Projects a function name function on to a given basis
// using quadrature.  The basis id determined to be either DG or Mixed

    int lastElement;

    if(DG == BijFlag)
    {
        lastElement = grid.NumElementsNoGhost;
    }
    else if(Mixed == BijFlag)
    {
        lastElement = grid.NumElementsWithGhost;
    }
    else
    {
        assert(false);
        return 1;
    }

    ddpSparseVector_type localWeights = sparseWeights;
    int numTestFunctions = FuncForward.size();
    ddpDenseVector_type RHSVector(numTestFunctions);


    //Build the vector of the functions when projected onto the
    // Chosen basis (DG or Mixed CG) from the functions point values

    int maxOrder;

    // loop over each element and find out what kind of basis function
    for(int elem = 0; elem < lastElement; ++elem)
    {
        //cout << "elem = " << elem << endl;

        if(DG == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderDG;
        }
        else if(Mixed == BijFlag)
        {
            maxOrder = (grid.ElementList[elem]).OrderMX;
        }
        else
        {
            assert(false);
            return 1;
        }

        // loop over the order of the polynomials of the basis function's
        // element
        //
        for(int order = 0; order <= maxOrder; ++order)
        {
            // get the global DOF from the local DOF
            //cout << "\torder = " << order << endl;
            std::pair<int, int> elemAndOrderJ;
            elemAndOrderJ.first = elem;
            elemAndOrderJ.second = order;
            int globalJ =  FuncForward.find(elemAndOrderJ)->second;

            ddpSparseVector_type temp = globalVandeMondeJ.col(globalJ);

            //We loop over the nonzero elements of temp1 and project function
            //onto the shape function/basis function by integrating that function
            //against the basis function
            double partial_sum = 0.0;
            for(ddpSparseVector_type::InnerIterator it(temp); it; ++it)
            {
                int ind = it.index();
                double tempval = temp.coeffRef(ind);
                double weightval = localWeights.coeffRef(ind);
                double FunctionVal = FunctionPointVals.coeffRef(ind);
                partial_sum += tempval * FunctionVal * weightval;
            }
            double valueJ = partial_sum;

            double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

            // only store that value in the vector if bigger than a certain
            // tolerance
            if(fabs(valueJ) > Threeeps)
            {
                RHSVector(globalJ) = valueJ;
            }
            else
            {
                RHSVector(globalJ) = 0.0;
            }
        } // order
    } // elem
    output = RHSVector;
}

