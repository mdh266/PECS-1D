#include "../include/includes.hpp"
#include "../include/basis.hpp"
#include "../include/main.hpp"
#include <iomanip>


// TODO:  The integration points are attached here in the make uniform grid,
//        This will cause problems if we ever want to use non-uniform grids.
//        The proper fix is to do something else to attach the grid points

int
ddpMakeUniformGrid(ddpProblemInfo_type const & problem,
                   ddpDomain_type const & domain,
                   ddpGrid_type & grid)
{
    // Make uniform grid. Although we make a uniform grid, what follows does
    // not require the grid to be uniform

    // Prepare quadrature point and weights
    int gaussLegendreNumPoints = problem.GaussLegendreNumPoints;

    // Get quadrature points and weights for Gauss and Legrendre
    // qudrature
    double * gaussLegendrePoints = new double[gaussLegendreNumPoints];
    double * gaussLegendreWeights = new double [gaussLegendreNumPoints];

    // contained in quadrule.cpp
    legendre_set (gaussLegendreNumPoints, gaussLegendrePoints,
                  gaussLegendreWeights);


//  lobatto_set (gaussLegendreNumPoints, gaussLegendrePoints,
//								gaussLegendreWeights);


    // Get the number of elements and right and left end point of the domain
    int maxOrderDG = problem.MaxOrderDG;
    int maxOrderMX = problem.MaxOrderMX;
    int numElementsNoGhost = grid.NumElementsNoGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;

    // Create TOTAL element ARRAY LIST
    ddpElement_type * tempElementList =
        new ddpElement_type[numElementsWithGhost];

    // Get the left and right end point of first subdomain and create a uniform
    // grid on it with numBoundaryElements number of elements
    double LocalLeftEndPoint = domain.LeftEndPoint;
    double LocalRightEndPoint = domain.RightEndPoint;
    int LocalNumElements = numElementsNoGhost;
    int NumElementsBefore = 0;

    // need to set the DeltaxMin first, can be arbitrary
    grid.DeltaxMin = 1.0;

    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);


    // We still need to add the upsilon element at the extreme right hand side
    // We do that here.  I.E. this is almost a ghost cell.

    ddpElement_type tempElement;

    // Add the extreme right hand side element
    // tempElement.Delta = NAN;
    // This is an ugly hack, need to fix this later.
    tempElement.Left =  domain.RightEndPoint;
    tempElement.Right = domain.RightEndPoint;
    tempElement.OrderMX = 0;
    tempElement.OrderDG = -1;
    tempElement.GaussLegendreNumPoints = -1;
    tempElement.GaussLegendrePoints = NULL;
    tempElement.GaussLegendreWeights = NULL;

    // Set the last element to the temp element
    tempElementList[numElementsWithGhost - 1] = tempElement;
    grid.ElementList = tempElementList;

    delete [] gaussLegendrePoints;
    delete [] gaussLegendreWeights;

    return 0;


}

int
ddpMakeRightRefinedBoundaryGrid(ddpProblemInfo_type const & problem,
                                ddpDomain_type const & domain,
                                ddpGrid_type & grid)
{

    // Prepare quadrature point and weights
    int gaussLegendreNumPoints = problem.GaussLegendreNumPoints;

    // Get quadrature points and weights for Gauss and Legrendre
    // qudrature
    double * gaussLegendrePoints = new double[gaussLegendreNumPoints];
    double * gaussLegendreWeights = new double [gaussLegendreNumPoints];

    // contained in quadrule.cpp
    legendre_set (gaussLegendreNumPoints, gaussLegendrePoints,
                  gaussLegendreWeights);

// lobatto_set (gaussLegendreNumPoints, gaussLegendrePoints,
    //							gaussLegendreWeights);


    // Get the number of elements and right and left end point of the domain
    int maxOrderDG = problem.MaxOrderDG;
    int maxOrderMX = problem.MaxOrderMX;
    int numElementsNoGhost = grid.NumElementsNoGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numBoundaryElements = grid.NumBoundaryElements;
    int numInteriorElements = numElementsNoGhost -numBoundaryElements;
    double boundaryLayerWidth = grid.BoundaryLayerWidth;


    // need to set the DeltaxMin first, can be arbitrary
    grid.DeltaxMin = 1.0;

// Create TOTAL element ARRAY LIST
    ddpElement_type * tempElementList =
        new ddpElement_type[numElementsWithGhost];

    // Get the left end point and the location of start of
    // the right boundary layer and create a uniform
    // grid on it with numBoundaryElements number of elements
    double LocalLeftEndPoint = domain.LeftEndPoint;
    double LocalRightEndPoint = domain.RightEndPoint - boundaryLayerWidth;
    int LocalNumElements = numInteriorElements;
    int NumElementsBefore = 0;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // Get the left and right end point of right boundary layer and create a uniform
    // grid on it with numBoundaryElements number of elements
    LocalLeftEndPoint = domain.RightEndPoint- boundaryLayerWidth;
    LocalRightEndPoint = domain.RightEndPoint;
    LocalNumElements = numBoundaryElements;
    NumElementsBefore = numInteriorElements;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // We still need to add the upsilon element at the extreme right hand side
    // We do that here.  I.E. this is almost a ghost cell.

    ddpElement_type tempElement;

    // Add the extreme right hand side element
    // tempElement.Delta = NAN;
    // This is an ugly hack, need to fix this later.
    tempElement.Left =  domain.RightEndPoint;
    tempElement.Right = domain.RightEndPoint;
    tempElement.OrderMX = 0;
    tempElement.OrderDG = -1;
    tempElement.GaussLegendreNumPoints = -1;
    tempElement.GaussLegendrePoints = NULL;
    tempElement.GaussLegendreWeights = NULL;

    // Set the last element to the temp element
    tempElementList[numElementsWithGhost - 1] = tempElement;
    grid.ElementList = tempElementList;

    delete [] gaussLegendrePoints;
    delete [] gaussLegendreWeights;

    return 0;
}

int
ddpMakeLeftRefinedBoundaryGrid(ddpProblemInfo_type const & problem,
                               ddpDomain_type const & domain,
                               ddpGrid_type & grid)
{

    // Prepare quadrature point and weights
    int gaussLegendreNumPoints = problem.GaussLegendreNumPoints;

    // Get quadrature points and weights for Gauss and Legrendre
    // qudrature
    double * gaussLegendrePoints = new double[gaussLegendreNumPoints];
    double * gaussLegendreWeights = new double [gaussLegendreNumPoints];

    // contained in quadrule.cpp
    legendre_set (gaussLegendreNumPoints, gaussLegendrePoints,
                  gaussLegendreWeights);

//	 lobatto_set (gaussLegendreNumPoints, gaussLegendrePoints,
//								 gaussLegendreWeights);


    // Get the number of elements and right and left end point of the domain
    int maxOrderDG = problem.MaxOrderDG;
    int maxOrderMX = problem.MaxOrderMX;
    int numElementsNoGhost = grid.NumElementsNoGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numBoundaryElements = grid.NumBoundaryElements;
    int numInteriorElements = numElementsNoGhost - numBoundaryElements;
    double boundaryLayerWidth = grid.BoundaryLayerWidth;


    // need to set the DeltaxMin first, can be arbitrary
    grid.DeltaxMin = 1.0;

// Create TOTAL element ARRAY LIST
    ddpElement_type * tempElementList =
        new ddpElement_type[numElementsWithGhost];

    // Get the left and right end point of left boundary and create a uniform
    // grid on it with numBoundaryElements number of elements
    double LocalLeftEndPoint = domain.LeftEndPoint;
    double LocalRightEndPoint = domain.LeftEndPoint + boundaryLayerWidth;
    int LocalNumElements = numBoundaryElements;
    int NumElementsBefore = 0;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // Get the left and right end point of right mesh and create a uniform
    // grid on it with numBoundaryElements number of elements
    LocalLeftEndPoint = domain.LeftEndPoint + boundaryLayerWidth;
    LocalRightEndPoint = domain.RightEndPoint;
    LocalNumElements = numInteriorElements;
    NumElementsBefore = numBoundaryElements;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // We still need to add the upsilon element at the extreme right hand side
    // We do that here.  I.E. this is almost a ghost cell.

    ddpElement_type tempElement;

    // Add the extreme right hand side element
    // tempElement.Delta = NAN;
    // This is an ugly hack, need to fix this later.
    tempElement.Left =  domain.RightEndPoint;
    tempElement.Right = domain.RightEndPoint;
    tempElement.OrderMX = 0;
    tempElement.OrderDG = -1;
    tempElement.GaussLegendreNumPoints = -1;
    tempElement.GaussLegendrePoints = NULL;
    tempElement.GaussLegendreWeights = NULL;

    // Set the last element to the temp element
    tempElementList[numElementsWithGhost - 1] = tempElement;
    grid.ElementList = tempElementList;

    delete [] gaussLegendrePoints;
    delete [] gaussLegendreWeights;

    return 0;
}
int
ddpMakeBothEndsRefinedGrid(ddpProblemInfo_type const & problem,
                           ddpDomain_type const & domain,
                           ddpGrid_type & grid)
{

    // Prepare quadrature point and weights
    int gaussLegendreNumPoints = problem.GaussLegendreNumPoints;

    // Get quadrature points and weights for Gauss and Legrendre
    // qudrature
    double * gaussLegendrePoints = new double[gaussLegendreNumPoints];
    double * gaussLegendreWeights = new double [gaussLegendreNumPoints];

    // contained in quadrule.cpp
    legendre_set (gaussLegendreNumPoints, gaussLegendrePoints,
                  gaussLegendreWeights);


    // Get the number of elements and right and left end point of the domain
    int maxOrderDG = problem.MaxOrderDG;
    int maxOrderMX = problem.MaxOrderMX;
    int numElementsNoGhost = grid.NumElementsNoGhost;
    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numBoundaryElements = grid.NumBoundaryElements;
    int numInteriorElements = numElementsNoGhost - 2 * numBoundaryElements;
    double boundaryLayerWidth = grid.BoundaryLayerWidth;


    // need to set the DeltaxMin first, can be arbitrary
    grid.DeltaxMin = 1.0;

// Create TOTAL element ARRAY LIST
    ddpElement_type * tempElementList =
        new ddpElement_type[numElementsWithGhost];

    // Get the left and right end point of left boundary and create a uniform
    // grid on it with numBoundaryElements number of elements
    double LocalLeftEndPoint = domain.LeftEndPoint;
    double LocalRightEndPoint = domain.LeftEndPoint + boundaryLayerWidth;
    int LocalNumElements = numBoundaryElements;
    int NumElementsBefore = 0;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // Get the left and right end point of interior mesh and create a uniform
    // grid on it with numBoundaryElements number of elements
    LocalLeftEndPoint = domain.LeftEndPoint + boundaryLayerWidth;
    LocalRightEndPoint = domain.RightEndPoint - boundaryLayerWidth;
    LocalNumElements = numInteriorElements;
    NumElementsBefore = numBoundaryElements;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);

    // Get the left and right end point of right boundary layer and create a uniform
    // grid on it with numBoundaryElements number of elements
    LocalLeftEndPoint = domain.RightEndPoint - boundaryLayerWidth;
    LocalRightEndPoint = domain.RightEndPoint;
    LocalNumElements = numBoundaryElements;
    NumElementsBefore = numInteriorElements + numBoundaryElements;


    ddpMakeUniformSubgrid(grid,
                          tempElementList,
                          maxOrderDG,
                          maxOrderMX,
                          LocalLeftEndPoint,
                          LocalRightEndPoint,
                          LocalNumElements,
                          NumElementsBefore,
                          gaussLegendrePoints,
                          gaussLegendreWeights,
                          gaussLegendreNumPoints);


    // We still need to add the upsilon element at the extreme right hand side
    // We do that here.  I.E. this is almost a ghost cell.

    ddpElement_type tempElement;

    // Add the extreme right hand side element
    // tempElement.Delta = NAN;
    // This is an ugly hack, need to fix this later.
    tempElement.Left =  domain.RightEndPoint;
    tempElement.Right = domain.RightEndPoint;
    tempElement.OrderMX = 0;
    tempElement.OrderDG = -1;
    tempElement.GaussLegendreNumPoints = -1;
    tempElement.GaussLegendrePoints = NULL;
    tempElement.GaussLegendreWeights = NULL;

    // Set the last element to the temp element
    tempElementList[numElementsWithGhost - 1] = tempElement;
    grid.ElementList = tempElementList;

    delete [] gaussLegendrePoints;
    delete [] gaussLegendreWeights;

    return 0;
}

int
ddpMakeUniformSubgrid(ddpGrid_type & grid,
                      ddpElement_type * tempElementList,
                      const int & maxOrderDG,
                      const int & maxOrderMX,
                      const double & LocalLeftEndPoint,
                      const double & LocalRightEndPoint,
                      const int & LocalNumElements,
                      const int & NumElementsBefore,
                      double * gaussLegendrePoints,
                      double * gaussLegendreWeights,
                      const int & gaussLegendreNumPoints)
{

    // Make uniform local mesh in each subdomain [LocalLeftEndPoint, LocalRightEndPoint]

    double LocalLength = LocalRightEndPoint - LocalLeftEndPoint;
    double LocaldeltaX = LocalLength / (LocalNumElements);

    // check to make sure that this subdomain doenst have the highest order
    if(grid.OrderDGMax < maxOrderDG)
        grid.OrderDGMax = maxOrderDG;

    // check to make sure that this subdomain doenst have the smallest dx
    if(grid.DeltaxMin > LocaldeltaX)
        grid.DeltaxMin = LocaldeltaX;

    // loop over the number of elements, and fill in each elements
    // information into a temporary element and then set the
    // elemen to that temporary element

    ddpElement_type tempElement;

    for(int i = 0; i < LocalNumElements; ++i)
    {
        // get and store each elements metadata
        tempElement.Left = LocalLeftEndPoint + (i ) * LocaldeltaX;
        tempElement.Right = LocalLeftEndPoint + (i+1) * LocaldeltaX;
        tempElement.Delta = tempElement.Right - tempElement.Left;
        tempElement.OrderMX = maxOrderMX;
        tempElement.OrderDG = maxOrderDG;
        tempElement.GaussLegendreNumPoints = gaussLegendreNumPoints;

        // Fill in quadrature point and weight values for Gauss-Lobatto
        // Legerendre quadrature.  This is scaled so that it is in each
        // elements end points

        tempElement.GaussLegendrePoints = new double[gaussLegendreNumPoints];
        tempElement.GaussLegendreWeights = new double[gaussLegendreNumPoints];

        for(int j = 0; j < gaussLegendreNumPoints; ++j)
        {
            (tempElement.GaussLegendrePoints)[j] = tempElement.Left + LocaldeltaX
                                                   / 2 * (gaussLegendrePoints[j] + 1);

            (tempElement.GaussLegendreWeights)[j] = 0.5 * (tempElement.Right -
                                                    tempElement.Left) * gaussLegendreWeights[j];

        }

        // Now place in the temp element list so that it is the correct spot
        tempElementList[NumElementsBefore + i] = tempElement;
    }

    return 0;
}

int
ddpGlueGrids(ddpProblemInfo_type const & problem,
             ddpGrid_type const & left_grid,
             ddpGrid_type const & right_grid,
             ddpGrid_type & GluedGrid)
{
    // copy the grid so that is has (-left_grid, right_grid)

    GluedGrid.DeltaxMin = std::min(left_grid.DeltaxMin, right_grid.DeltaxMin);
    GluedGrid.OrderDGMax = std::max(left_grid.OrderDGMax, right_grid.OrderDGMax);


    GluedGrid.Domain.LeftEndPoint =  left_grid.Domain.LeftEndPoint;
    GluedGrid.Domain.RightEndPoint = right_grid.Domain.RightEndPoint;
    GluedGrid.NumElementsNoGhost = left_grid.NumElementsNoGhost + right_grid.NumElementsNoGhost;
    GluedGrid.NumElementsWithGhost = GluedGrid.NumElementsNoGhost + 1;

    // needs to be the same
    GluedGrid.NumBoundaryElements = left_grid.NumBoundaryElements;
    GluedGrid.BoundaryLayerWidth = left_grid.BoundaryLayerWidth;

    // create new array of elements that is twice as long
    ddpElement_type * tempElementList =
        new ddpElement_type[GluedGrid.NumElementsWithGhost];

//	cout << "GluedGrid NumElementsWithGhost = " << GluedGrid.NumElementsWithGhost << endl;

    int copierIndex = 0;
    ddpElement_type tempElement;

//	cout << "First copy over" << endl;
    for(int i = 0; i < left_grid.NumElementsNoGhost; i++)
    {
        // copy over eveything here from the left grid should be the same
        tempElement.Delta = left_grid.ElementList[i].Delta;
        tempElement.OrderMX = left_grid.ElementList[i].OrderMX;
        tempElement.OrderDG = left_grid.ElementList[i].OrderDG;
        tempElement.GaussLegendreNumPoints =
            left_grid.ElementList[i].GaussLegendreNumPoints;

        // Left end point of element should = -1 x right end of grid element
        // right end point of element should = -1 x left end of grid element
        tempElement.Left = left_grid.ElementList[i].Left;
        tempElement.Right = left_grid.ElementList[i].Right;

        int N = left_grid.ElementList[i].GaussLegendreNumPoints;
//		cout << "LEGENDRE: N  = " << N << endl;
        tempElement.GaussLegendrePoints = new double[N];
        tempElement.GaussLegendreWeights = new double[N];

        // copy over the guass-legendre points
        for(int j = 0; j < N; j++)
        {
            tempElement.GaussLegendrePoints[j] =
                left_grid.ElementList[i].GaussLegendrePoints[j];
        }
        // copy over the guass-legendre weights
        for(int j = 0; j < N; j++)
        {
            tempElement.GaussLegendreWeights[j] =
                left_grid.ElementList[i].GaussLegendreWeights[j];
        }

//		cout << "Copier index = " << copierIndex << endl;
        tempElementList[copierIndex] = tempElement;

        // update copier index
        copierIndex++;
    }


//	cout << "2nd copy over" << endl;
    // copy over everything to form right half of grid
    for(int i = 0; i < right_grid.NumElementsNoGhost; i++)
    {
//		cout << "i = " << i << endl;

        // copy over eveything here from the left grid should be the same
        tempElement.Delta = right_grid.ElementList[i].Delta;
        tempElement.OrderMX = right_grid.ElementList[i].OrderMX;
        tempElement.OrderDG = right_grid.ElementList[i].OrderDG;
        tempElement.GaussLegendreNumPoints =
            right_grid.ElementList[i].GaussLegendreNumPoints;

        // Left end point of element should = left end of grid element
        // right end point of element should = right end of grid element
        tempElement.Left = right_grid.ElementList[i].Left;
        tempElement.Right = right_grid.ElementList[i].Right;

        int	N = right_grid.ElementList[i].GaussLegendreNumPoints;
//		cout << "LEGENDRE: N  = " << N << endl;
        tempElement.GaussLegendrePoints = new double[N];
        tempElement.GaussLegendreWeights = new double[N];

        // copy over the guass-legendre points
        for(int j = 0; j < N; j++)
        {
            tempElement.GaussLegendrePoints[j] =
                right_grid.ElementList[i].GaussLegendrePoints[j];
        }
        // copy over the guass-legendre weights
        for(int j = 0; j < N; j++)
        {
            tempElement.GaussLegendreWeights[j] =
                right_grid.ElementList[i].GaussLegendreWeights[j];
        }

//		cout << "Copier Index = " << copierIndex << endl;
        // copy over this element to the list
        tempElementList[copierIndex] = tempElement;

        // update the copier index
        copierIndex++;
    }

    // tempElement.Delta = NAN;
    // This is an ugly hack, need to fix this later.
    tempElementList[copierIndex].Left =  GluedGrid.Domain.RightEndPoint;
    tempElementList[copierIndex].Right = GluedGrid.Domain.RightEndPoint;
    tempElementList[copierIndex].OrderMX = 0;
    tempElementList[copierIndex].OrderDG = -1;

    // Add the extreme right hand side element
    tempElementList[copierIndex].GaussLegendreNumPoints = -1;
    tempElementList[copierIndex].GaussLegendrePoints = NULL;
    tempElementList[copierIndex].GaussLegendreWeights = NULL;

    // Now copy over the element list
    GluedGrid.ElementList = tempElementList;
//	cout << "Here " << endl;

    return 0;
}


int
ddpMakeTimeStamps(ddpProblemInfo_type const & problem,
                  std::vector< double > & timeStamps)
{

    // to collect the right times that we should be outputtin data

    double tInitial;
    if(Original == problem.IVP_Type)
    {
        tInitial = problem.TimeInitial;
    }
    else
    {
        int N = problem.NumFrames -1;
        char fileNumber[3];
        sprintf(fileNumber, "%4.4d", N);

        std::stringstream ss;
        ss << "State" << fileNumber << ".dat";
        std::string str = ss.str();
        std::ifstream read;
        read.open(str.c_str());
        double time;
        int counter = 0;
        while(!read.eof() && counter < 8)
        {
            read >> time;
            counter++;
        }

        tInitial = time;
        read.close();
    }
    double tFinal = problem.TimeFinal;
    cout << "tInitial: " << tInitial << endl;
    cout << "tFinal: " << tFinal << endl;



    int numFrames = problem.NumFrames;

//  cout << "numFrames = " << numFrames << endl;
    timeStamps.reserve(numFrames);

    double FrameDeltaT = (tFinal - tInitial) / (numFrames - 1);
    // Kents Note:  think about the right way to handle the case where we only want
    //        one frame, the initial condition.
    //        The "minus one" is for the fact that we will output the
    //        initial condition.
    //  Solution:  Don't think in terms of frames, think in terms of time steps
    //             (fencepost problem)

    //double * timeStampsTemp = new double[numFrames];
    for(int j = 0; j < numFrames; ++j)
    {
        timeStamps.push_back( tInitial + FrameDeltaT * j);
    }

    return 0;
}


int
ddpMakeBijs(ddpGrid_type const & grid,
            ddpBijFlag_type const & flag,
            ddpBijForw_type& Forward,
            ddpBijBack_type& Backward)
{
// This makes the bijections from local DOFs
// To global dofs FOR The SPECIFIC FEM METHOD

    int numElements = grid.NumElementsWithGhost;
    ddpElement_type currentElement;
    int counter = 0;
    int LocalMax = 0;

    std::pair < int, int > tempPair;
    std::pair < std::pair < int, int >, int > tempDataForward;
    std::pair < int, std::pair <int, int > > tempDataBackward;

    // Here, we loop over every element in the grid,
    // and depending on how the flag is set, we perform an inner loop.

    // Note that we loop over all the elements, including the ghost elements.
    // This doens't cause a problem because of the way the element data is set
    // for the ghost element.
    for(int i = 0; i < numElements; ++i)
    {
        /// get the type of FEM method and elements local maximum poly order
        currentElement = grid.ElementList[i];
        switch(flag)
        {
        case DG:
            LocalMax = currentElement.OrderDG+1; // plus one is on purpose
            break;
        case Mixed:
            LocalMax = currentElement.OrderMX+1;  // plus one is on purpose
            break;
        case Points:
            LocalMax = (currentElement.GaussLegendreNumPoints);
            // lack of plus one is on purpose
            break;
        default:
            cout << "Can't happen on line number "  << __LINE__ << endl;
        }

        // get the local and global degree of freedom and put into the bijective
        // MAP
        for(int j = 0; j < LocalMax; ++j)
        {
            tempPair.first  = i;
            tempPair.second = j;

            tempDataForward.first   = tempPair;
            tempDataForward.second  = counter;

            tempDataBackward.first  = counter;
            tempDataBackward.second = tempPair;

            Forward.insert (tempDataForward );
            Backward.insert (tempDataBackward);

            ++counter;
        }
    }

    return 0;
}

int
ddpMakeWeightsAndPoints(ddpGrid_type const & grid,
                        ddpProblemInfo_type const & problem,
                        ddpSparseVector_type & weightsSparse,
                        ddpSparseVector_type & PTSSparse,
                        ddpDenseVector_type & PTSDense,
                        ddpDenseVector_type & weightsDense)
{

// Makes weights and points vectors for quadrature

    int weightVectorLength = (grid.NumElementsNoGhost)
                             * problem.GaussLegendreNumPoints;

// cout << "Entering MakeWeightVectorSparse" << endl;
    ddpSparseVector_type
    temp_weightsSparse(weightVectorLength),
                       temp_PTSSparse(weightVectorLength);

    ddpDenseVector_type
    temp_PTSDense(weightVectorLength),
                  temp_weightsDense(weightVectorLength);

    ddpMakeWeightVectorSparse(grid, temp_weightsSparse);

    temp_weightsDense = temp_weightsSparse;

//	    cout << "\tExiting MakeWeightVectorSparse" << endl;

//	  cout << "Entering MakePTVectorSparse"  << endl;
    ddpMakePTVectorSparse(grid, temp_PTSSparse);

//	  cout << "\tExiting MakePTVectorSparse" << endl;

//	  cout << "Entering MakePTVectorDense" << endl;

    ddpMakePTVectorDense(grid, temp_PTSDense);
    //cout << "PTSDense = " << PTSDense << endl;

//	  cout << "\tExiting MakePTVectorSparse" << endl;

    weightsSparse = temp_weightsSparse;
    PTSSparse = temp_PTSSparse;
    PTSDense = temp_PTSDense;
    weightsDense = temp_weightsDense;

    return 0;
}

int
ddpMakeAllBijections(ddpGrid_type const & grid,
                     ddpBijection_type & Bijections)
{

// Make bijection depending on flag type.  Sets flag type
// before calling subroutine.

    ddpBijFlag_type flag;

//	  cout << "Entering MakeBijs" << endl;
    flag = DG;
    ddpMakeBijs(grid, flag, Bijections.DGForward, Bijections.DGBackwrd);

    flag = Mixed;
    ddpMakeBijs(grid, flag, Bijections.MXForward, Bijections.MXBackwrd);

    flag = Points;
    ddpMakeBijs(grid, flag, Bijections.PTForward, Bijections.PTBackwrd);
//	  cout << "\tExiting MakeBijs" << endl;

    return 0;
}


int
ddpMakeWeightVectorDense(ddpGrid_type const & grid,
                         ddpDenseVector_type &weights)
{
    /// Make a dense vector of all the domains quadrature weights
    ddpDenseVector_type temp1;
    Eigen::VectorXd temp;

    int numElementsNoGhost = grid.NumElementsNoGhost;
    std::vector<double> weightsVector;
    int pointCount = 0;

    for(int i = 0; i < numElementsNoGhost; ++i)
    {
        int numLocalPoints = (grid.ElementList[i]).GaussLegendreNumPoints;
        for(int j = 0; j < numLocalPoints; ++j)
        {
            double localWeight = (grid.ElementList[i]).GaussLegendreWeights[j];
            weights(pointCount) = localWeight;
            pointCount++;
        }
    }

    return 0;
}

int
ddpMakePTVectorDense(ddpGrid_type const & grid,
                     ddpDenseVector_type & points)
{

    // Make a dense vector of the all the domains quadrature pointts
    int numElements = grid.NumElementsNoGhost;

    int pointCount = 0;

    for(int i = 0; i < numElements; ++i)
    {
        int numLocalPoints = (grid.ElementList[i]).GaussLegendreNumPoints;
        for(int j = 0; j < numLocalPoints; ++j)
        {
            double localPT = (grid.ElementList[i]).GaussLegendrePoints[j];
            points(pointCount)  = localPT;
            pointCount++;
        }
    }
    return 0;
}



int
ddpMakeWeightVectorSparse(ddpGrid_type const & grid,
                          ddpSparseVector_type & weights)
{

    /// Make a sparse vector of all the domains quadrature weights
    int numElements = grid.NumElementsNoGhost;
    int pointCount = 0;

    for(int i = 0; i < numElements; ++i)
    {
        int numLocalPoints = (grid.ElementList[i]).GaussLegendreNumPoints;
        for(int j = 0; j < numLocalPoints; ++j)
        {
            double localWeight = (grid.ElementList[i]).GaussLegendreWeights[j];
            weights.insert(pointCount)  = localWeight;
            pointCount++;
        }
    }

    return 0;
}

int
ddpMakePTVectorSparse(ddpGrid_type const & grid,
                      ddpSparseVector_type & points)
{
    // Make a sparse vector of the all the domains quadrature pointts

    int numElements = grid.NumElementsNoGhost;
    int pointCount = 0;

    for(int i = 0; i < numElements; ++i)
    {
        int numLocalPoints = (grid.ElementList[i]).GaussLegendreNumPoints;
        for(int j = 0; j < numLocalPoints; ++j)
        {
            double localPT = (grid.ElementList[i]).GaussLegendrePoints[j];
            points.insert(pointCount)  = localPT;
            pointCount++;
        }
    }

    return 0;
}



