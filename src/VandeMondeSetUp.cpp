#include "../include/VandeMondeSetUp.hpp"


int
ddpMakeGlobalVandeMonde(ddpGrid_type const & grid,
                        ddpProblemInfo_type const & problem,
                        ddpBijFlag_type const & flag,
                        ddpBijForw_type const & FuncForward,
                        ddpBijForw_type const & PTForward,
                        ddpNumDeriv_type const & numDeriv,
                        ddpSparseMatrix_type & vandeMonde)
{
    // Makes general VandeMondeMatrix for each basis function
    // and its deriviatves.  Takes into account whether
    // DG or Mixed CG.
    //
    //  When you multiply the DOF of a solution by this matirx
    //  it will give you the solutions point values.

    // M(I,J) = \Psi_{J}(x_I)
    // J = global degree of freedom depending on element and the LOCAL
    // maximum polynomial order of that element
    // I = mesh point
    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numElementsNoGhost = grid.NumElementsNoGhost;


    int localOrderMax;
    int localNumPoints;
    long int numDofs = FuncForward.size();
    long int numPtsPerDof = problem.GaussLegendreNumPoints;
    long int numPts = PTForward.size();


    std::pair<int, int> elemAndOrderJ;
    std::pair<int, int> elemAndPointI;

    int GlobalI;  // global matrix row
    int GlobalJ; // global matrix column
    double ValueIJ; // value at matrix(GlobalI,GlobalJ)
    double pt;

    // Triplets is just used to build the sparse matrix
    // it stores the (row, col, value) as a list
    // and then we will use it to fill the sparse matrix later

    typedef Eigen::Triplet < double > T;
    std::vector<T> triplets;

    ddpElement_type localElement;

    long long int estimateofsize;
    estimateofsize = numDofs * numPtsPerDof;

    triplets.reserve(estimateofsize);

    Eigen::SparseMatrix<double> vandeMondeAux(numPts,numDofs);


    // Loop over just the nearest neighbors
    // NOTE:  This supposes support is just nearest neighbors
    for(int elem1 = 0; elem1 < numElementsWithGhost; ++elem1)
    {

        for( int elem2 = std::max(0,elem1 - 1);
                elem2 < std::min(numElementsNoGhost, elem1 + 1); ++elem2)
        {
            // Figure out if your doing Mixed CG or DG
            switch (flag)
            {
            case DG:
                localOrderMax = grid.ElementList[elem2].OrderDG;
                break;
            case Mixed:
                localOrderMax = grid.ElementList[elem2].OrderMX;
                break;
            default:
                cout << "this cant happen on line" << __LINE__ << endl;
                return 1;
            }

            localElement = grid.ElementList[elem2];

            // loop over the orders of your basis functions
            for(int order = 0; order <= localOrderMax; ++order)
            {
                localNumPoints = localElement.GaussLegendreNumPoints;

                // loop over the number of points in each element
                for(int pointNumber = 0; pointNumber < localNumPoints; ++pointNumber)
                {
                    elemAndOrderJ.first = elem1;
                    elemAndOrderJ.second = order;
                    elemAndPointI.first = elem2;
                    elemAndPointI.second = pointNumber;

                    pt = grid.ElementList[elem2].GaussLegendrePoints[pointNumber];
//			 			 cout << "pt = " << pt << endl;
                    // These are the rows and columns correspond to point
                    // value in the mesh while the columns the global DOFs
                    GlobalJ = FuncForward.find(elemAndOrderJ)->second;
                    GlobalI = PTForward.find(elemAndPointI)->second;

                    // Heres the work where it decides the basis function and
                    // number of derivativs.
                    switch(flag)
                    {
                    case DG:
                        ValueIJ = ddpPsi(grid, elem1, order, numDeriv, pt);
                        break;
                    case Mixed:
                        ValueIJ = ddpUpsilon(grid, elem1, order, numDeriv,pt);
                        break;
                    default:
                        cout << "this shouldn't happen on line"  << __LINE__ << endl;
                        return 1;
                        break;
                    }

                    // only store the value if its above the machine epsilon
                    double absValueIJ = fabs(ValueIJ);
                    double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

//						cout << "pt = " << pt << endl;

//						cout << "(" << GlobalI << ", " << GlobalJ << ", " << ValueIJ << ")" << endl;
                    if( absValueIJ > Threeeps )
                    {
                        triplets.push_back( T(GlobalI,GlobalJ,ValueIJ) );
                    }

                } // loop over pointNumber
            } // loop over order
        } // loop over elem2
    } //loop over elem1

//	cout << "here " << endl;
    // convert from triplets to sparse matrix
    vandeMondeAux.setFromTriplets(triplets.begin(), triplets.end() );
//	cout << "here2" << endl;
    vandeMondeAux.makeCompressed();
    vandeMonde = vandeMondeAux;

    return 0;
}



int
ddpMakeGlobalVandeMondeFluxDG(ddpGrid_type const & grid,
                              ddpProblemInfo_type const & problem,
                              ddpBijForw_type const & FuncForward,
                              ddpDirection_type const & Flag1,
                              ddpDirection_type const & Flag2,
                              ddpSparseMatrix_type & vandeMonde)
{

    // Makes a VandeMondeMatrix for the Fluxe values of the
    // DG basis functions and its deriviatve.
    // When you multiply the functions DOF by this matrix you will get
    // the function flux values at the elements interface nodes.
    //

    // M(I,J) = \Psi_{J}(x_I)
    // J = global degree of freedom depending on element and the LOCAL
    // maximum polynomial order of that element
    // I = end point of the elements

    // A COLUMN IS THE DIFFERENT BASIS FUNCTION
    // A ROW IS THE CHOSEN END POINT OF THE ELEMENT
    // OVER WHICH THE BASIS HAS SUPPORT

    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numElementsNoGhost = grid.NumElementsNoGhost;

    int numFunctions = FuncForward.size();
    int numEndpoints = numElementsNoGhost + 1;

    ddpSparseMatrix_type vandeMondeAux(numEndpoints - 1,numFunctions);

    std::pair<int, int> elemAndOrderI;
    std::pair<int, int> elemAndOrderJ;

    int GlobalI;
    int GlobalJ;
    double ValueIJ;

    // Triplets is just used to build the sparse matrix
    // it stores the (row, col, value) as a list
    // and then we will use it to fill the sparse matrix later
    typedef
    Eigen::Triplet < double > T;

    std::vector<T> triplets;

    long long int estimateofsize;
    estimateofsize = numEndpoints;
    triplets.reserve(3 * estimateofsize);

    int endpointNumber;

    // loop  over the elements
    for(int elem1 = 0; elem1 < numElementsNoGhost; ++elem1)
    {

        // get the maximum polynomial order for the basis functions on
        // this element
        int MaxDGOrder = grid.ElementList[elem1].OrderDG;

        // loop over the nearest neightbors to find out the basis function
        // (elem1, order)'s flux across its end points

        // NOTE: Assumes that Legrendre basis function with support over
        // only nearest neighbors, if the basis functions
        // have larger support will need to change.
        for(int elem2 = std::max(0, elem1 - 1 );
                elem2 <= std::min(numElementsNoGhost - 1, elem1 + 1);
                ++elem2)
        {

            // loop over the polynomial order of the basis function
            // on elem1
            for(int order = 0; order <= MaxDGOrder; ++order)
            {

                // This is used to get left pt value or right pt
                // value of elem2's boundarys
                int offset;


                // Flag1 =  Plus -> We are taking the flux across
                // 		  the right sides of the elements
                // Flag1 =  Minus -> We are taking the flux across
                // 		  the left sides of the elements
                if(Plus == Flag1)
                {
                    offset = 1;

                }
                else if(Minus == Flag1)
                {
                    offset = 0;
                }
                else
                {
                    cout << "this can't happen on line "  << __LINE__ << endl;
                    assert(false);
                }

                // Get the global DOF from the basis function on local dof
                // (elem1, order)
                elemAndOrderJ.first = elem1;
                elemAndOrderJ.second = order;
                GlobalJ = FuncForward.find(elemAndOrderJ)->second;

                // Get the point value of chosen side the elem2's boundary
                GlobalI = elem2;
                endpointNumber = elem2 + offset;

                // get the flux of the polynomial of (elem, order) across the
                // elem's chosen boundary from the direction Flag2
                // Flag2 = Minus From Left
                // Flag2 = Plus from Right
                ValueIJ = ddpLimPsi(grid, elem1, order, Flag2, endpointNumber);

                // only store the value if its above the machine epsilon
                double absValueIJ = fabs(ValueIJ);
                double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

                if( absValueIJ > Threeeps )
                {
                    triplets.push_back( T(GlobalI,GlobalJ,ValueIJ) );
                }

            } // loop over order
        } // loop over elem2
    } // loop over elem1

    /// covert from triplets list to sparse matrix format.
    vandeMondeAux.setFromTriplets(triplets.begin(), triplets.end() );
    vandeMondeAux.makeCompressed();
    vandeMonde = vandeMondeAux;

    return 0;
}

int
ddpMakeGlobalVandeMondeFluxMX(ddpGrid_type const & grid,
                              ddpProblemInfo_type const & problem,
                              ddpBijForw_type const & FuncForward,
                              ddpSparseMatrix_type & vandeMonde)
{
    // Makes a VandeMondeMatrix for the values of the
    // Mixed CG basis functions at the element end points
    //
    // When you multiply the Electric field DOF by this
    // matrix you will get the Electric field point values
    //  values at the elements end points.
    //

    // M(I,J) = \Upsilon_{J}(x_I)
    // J = global degree of freedom depending on element and the LOCAL
    // maximum polynomial order of that element
    // I = end point of the elements

    // NOTE: There is no real flux for the electric field since its a continuous
    // function.  However, we this matrix is used to get the value of the
    // electric field at the elements end point and is USED in making the
    // drift/convective flux.  This is the reason for the
    // the term flux in the function title.

    // A COLUMN IS THE DIFFERENT BASIS FUNCTION
    // A ROW IS THE CHOSEN END POINT OF THE ELEMENT
    // OVER WHICH THE BASIS HAS SUPPORT

    int numElementsWithGhost = grid.NumElementsWithGhost;
    int numElementsNoGhost = grid.NumElementsNoGhost;

    int numFunctions = FuncForward.size();
    int numEndpoints = numElementsNoGhost + 1;

    ddpSparseMatrix_type vandeMondeAux(numEndpoints, numFunctions);

    // Triplets is just used to build the sparse matrix
    // it stores the (row, col, value) as a list
    // and then we will use it to fill the sparse matrix later
    typedef
    Eigen::Triplet < double > T;

    std::vector<T> triplets;

    long long int estimateofsize;
    estimateofsize = numEndpoints;

    triplets.reserve(estimateofsize);

    // loop over the number of elemens
    for(int elem = 0; elem < numElementsWithGhost - 1; ++elem)
        // minus one on NumElements is on purpose to deal with the last
        // element specifically.
    {

        int MaxMXOrder = grid.ElementList[elem].OrderMX;

        // loop over the polynomial order of each elements basis function
        for(int order = 0; order <= MaxMXOrder; ++order)
        {
            std::pair<int, int> elemAndOrderJ;

            int GlobalI;
            int GlobalJ;
            double ValueIJ;
            elemAndOrderJ.first = elem;
            elemAndOrderJ.second = order;

            double pt = grid.ElementList[elem].Left;
            // loop over all "possible" points
            // "possible" in here means possible for upsilon to be
            // nonzero.

            GlobalJ = FuncForward.find(elemAndOrderJ)->second;
            GlobalI = elem;


            ddpNumDeriv_type nD = zero;
            ValueIJ = ddpUpsilon(grid, elem, order, nD, pt);

            //This is an ugly hack to account for the overcounting
            // of the interior
            // element endpoints.  A less hacky way to do this
            // would be to
            // implement some sort of "limit" functional like
            // in Mathematica.
            // For now, this works.
            if(0 != elem)
            {
                ValueIJ *= 0.5;
            }

            // only store values that are bigger than machine epsilon
            double absValueIJ = fabs(ValueIJ);
            double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

            if( absValueIJ > Threeeps )
            {
                triplets.push_back( T(GlobalI,GlobalJ,ValueIJ) );
            }

        } // loop over order
    } // loop over elements

    //Need to add the last function, the last "half tent" function;
    int elem = numElementsWithGhost - 1; //The last element
    int order = 0;

    int GlobalI;
    int GlobalJ;
    double ValueIJ;
    std::pair<int, int> elemAndOrderJ;
    double pt = (grid.ElementList[numElementsNoGhost - 1 ]).Right;

    elemAndOrderJ.first = elem;
    elemAndOrderJ.second = order;

    // The upsilon we will consider here is the last function,
    // half of its domain is on the last element and
    // the other half of its domain is the "ghost element".
    GlobalI = numElementsNoGhost;
    GlobalJ = FuncForward.find(elemAndOrderJ)->second;
    ddpNumDeriv_type nD = zero;
    ValueIJ = ddpUpsilon(grid, elem, order, nD, pt);


    // only store values which are above machine epislon
    double absValueIJ = fabs(ValueIJ);
    double Threeeps = 12.0*std::numeric_limits<double>::epsilon();

    if( absValueIJ > Threeeps )
    {
        triplets.push_back( T(GlobalI,GlobalJ,ValueIJ) );
    }


    // convert from list of triplets to compressed matrix
    vandeMondeAux.setFromTriplets(triplets.begin(), triplets.end() );
    vandeMondeAux.makeCompressed();
    vandeMonde = vandeMondeAux;

    return 0;
}

int
ddpMakeVandeMondeMatrices(ddpGrid_type const & grid,
                          ddpProblemInfo_type const & problem,
                          ddpBijection_type const & Bijections,
                          ddpVandeMondeMatrices_type & VandeMondeMatrices,
                          ddpDGFluxMatrices_type & DGFluxMatrices)
{
    // This function will create all the  global vandemonde matrices
    // for mixed and dg spaces that are used in solving the
    // drift diffusion-Poisson equations. It calls
    // the functions defined above this.


    ddpSparseMatrix_type
    temp_globalVandeMondeDG,
    temp_globalVandeMondeDGPrime,
    temp_globalVandeMondeMX,
    temp_globalVandeMondeMXPrime,
    temp_globalVandeMondeFluxMX,
    workingMatrix;

    ddpBijFlag_type flag;
    ddpNumDeriv_type numDeriv;

    //////////////////////////////////////////////////////////////////
    // Make Regualar VandeMonde Matrices
    //////////////////////////////////////////////////////////////////

//	  cout << "entering MakeGlobalVandeMonde" << endl;

    // DG basis functions vandemonde matrix
    // M(I,J) = \psi_{J}(x_{I})
    // J = Global DOF for (elem, order)
    // I = point value in mesh
    flag = DG;
    numDeriv = zero;
    ddpMakeGlobalVandeMonde(grid, problem, flag,
                            Bijections.DGForward, Bijections.PTForward, numDeriv,
                            temp_globalVandeMondeDG);


    // Derivative of DG basis functions vandemonde matrix
    // M(I,J) = \frac{d \psi_{J}(x_{I})}{dx}
    //
    // J = Global DOF for (elem, order)
    // I = point value in mesh
    flag = DG;
    numDeriv = one;
    ddpMakeGlobalVandeMonde(grid, problem, flag,
                            Bijections.DGForward, Bijections.PTForward, numDeriv,
                            temp_globalVandeMondeDGPrime);


    // Mixed CG basis functions vandemonde matrix
    // M(I,J) = \Upsilon_{J}(x_{I})
    // J = Global DOF for (elem, order)
    // I = point value in mesh
    flag = Mixed;
    numDeriv = zero;
    ddpMakeGlobalVandeMonde(grid, problem, flag,
                            Bijections.MXForward, Bijections.PTForward, numDeriv,
                            temp_globalVandeMondeMX);


    // Dertivative of Mixed CG basis functions vandemonde matrix
    // M(I,J) = \frac{d\Upsilon_{J}(x_{I})}{dx}
    // J = Global DOF for (elem, order)
    // I = point value in mesh
    flag = Mixed;
    numDeriv = one;
    ddpMakeGlobalVandeMonde(grid, problem, flag,
                            Bijections.MXForward, Bijections.PTForward, numDeriv,
                            temp_globalVandeMondeMXPrime);

//	  cout << "\tExiting MakeGlobalVandeMonde" << endl;
//	  cout << "Entering MakeGlobalVandeMondeFluxMX" << endl;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////


    // Mixed CG basis functions vandemonde matrix at element end point
    // M(I,J) = \Upsilon_{J}(x_{I})
    // J = Global DOF for (elem, order)
    // I = point value elements end point
    ddpMakeGlobalVandeMondeFluxMX(grid,
                                  problem,
                                  Bijections.MXForward,
                                  temp_globalVandeMondeFluxMX);


//	 cout << "\tExiting MakeGlobalVandeMondeFluxMX" << endl;
//	  cout << "Entering MakeGlobalVandeMondeFluxDG" << endl;

    // Vandemonde matrices corresponding to the flux of baiss functions
    // across left or right element end point and from left or right
    // Flag1 = Plus -> Right Side of Element
    // Flag1 = Minus -> Left Side of Element
    // Flag2 = Plus -> limit/flux from right
    // Flag2 = Minus -> limit/flux from left


    ddpDGFluxMatrices_type temp_DGFluxMatrices;

    // Flux of basis function at right end point of elements from
    // the right.
    // M(i,j) = lim_{x-> x_{i+1/2} && x > x_{i+1/2}} \psi_{j}(x)
    // j = Global DOF for (elem, order)
    // i = point value elements end point
    ddpDirection_type Flag1 = Plus;
    ddpDirection_type Flag2 = Plus;

    ddpMakeGlobalVandeMondeFluxDG(grid,
                                  problem,
                                  Bijections.DGForward,
                                  Flag1,
                                  Flag2,
                                  workingMatrix);

    temp_DGFluxMatrices.plusplus = workingMatrix;
    temp_DGFluxMatrices.plusplusTransposed = workingMatrix.transpose();


    // Flux of basis function at right end point of elements from
    // the left.
    // M(i,j) = lim_{x-> x_{i+1/2} && x < x_{i+1/2}} \psi_{j}(x)
    // j = Global DOF for (elem, order)
    // i = point value elements end point
    Flag1 = Plus;
    Flag2 = Minus;

    ddpMakeGlobalVandeMondeFluxDG(grid,
                                  problem,
                                  Bijections.DGForward,
                                  Flag1,
                                  Flag2,
                                  workingMatrix);

    temp_DGFluxMatrices.plusminus = workingMatrix;
    temp_DGFluxMatrices.plusminusTransposed = workingMatrix.transpose();


    // Flux of the basis function at the left end point of elements from
    // the right.
    // M(i,j) = lim_{x-> x_{i-1/2} && x > x_{i-1/2}} \psi_{j}(x)
    // j = Global DOF for (elem, order)
    // i = point value elements end point
    Flag1 = Minus;
    Flag2 = Plus;

    ddpMakeGlobalVandeMondeFluxDG(grid,
                                  problem,
                                  Bijections.DGForward,
                                  Flag1,
                                  Flag2,
                                  workingMatrix);

    temp_DGFluxMatrices.minusplus = workingMatrix;
    temp_DGFluxMatrices.minusplusTransposed = workingMatrix.transpose();


    // Flux of the basis function at the left end point of elements from
    // the left.
    // M(i,j) = lim_{x-> x_{i-1/2} && x < x_{i-1/2}} \psi_{j}(X)
    // j = Global DOF for (elem, order)
    // i = point value elements end point
    Flag1 = Minus;
    Flag2 = Minus;

    ddpMakeGlobalVandeMondeFluxDG(grid,
                                  problem,
                                  Bijections.DGForward,
                                  Flag1,
                                  Flag2,
                                  workingMatrix);


    temp_DGFluxMatrices.minusminus = workingMatrix;
    temp_DGFluxMatrices.minusminusTransposed = workingMatrix.transpose();

//	  cout << "\tExiting MakeGlobalVandeMondeFluxDG" << endl;

    // regular version
    VandeMondeMatrices.globalVandeMondeDG = temp_globalVandeMondeDG;
    VandeMondeMatrices.globalVandeMondeDGPrime = temp_globalVandeMondeDGPrime;
    VandeMondeMatrices.globalVandeMondeMX = temp_globalVandeMondeMX;
    VandeMondeMatrices.globalVandeMondeMXPrime = temp_globalVandeMondeMXPrime;
    VandeMondeMatrices.globalVandeMondeFluxMX =  temp_globalVandeMondeFluxMX;

    //NOTE:  Storing the transposed versions at the beginnning
    //	   gives significant speed up instead of computing them
    //	   when used in calculating convective flux contribution at
    //	   every time step

    // transposed versions
    VandeMondeMatrices.globalVandeMondeDGTransposed =
        temp_globalVandeMondeDG.transpose();
    VandeMondeMatrices.globalVandeMondeDGPrimeTransposed =
        temp_globalVandeMondeDGPrime.transpose();

    DGFluxMatrices = temp_DGFluxMatrices;

    return 0;
}
