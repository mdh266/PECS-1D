#include "../include/basis.hpp"


double
ddpPsi( ddpGrid_type const & grid,
        int const & elem,
        int const & order,
        ddpNumDeriv_type const & numDerivs,
        double const & x)
{

    // This returns the value of the elem-th DG basis function of order order
    // with numDerive number of derivatives taken at the point x in the grid
    //
    // NOTE:  WE HAVE SPECIFICALLY ASSUMED LEGRENDRE BASIS FUNCTIONS

    double xLeft = (grid.ElementList[elem]).Left;
    double xRight = (grid.ElementList[elem]).Right;
    double Delta = (grid.ElementList[elem]).Delta;
    double xScaled = 2.0 * (x - xLeft) / Delta - 1;
    double temp1 = 0;

//	cout << "elem = " << elem << endl;
//	cout << "order = " << order << endl;
//	cout << "xLeft = " << xLeft << endl;
//	cout << "xRight = " << xRight << endl;

//  This is the jacobian to account for the change of variables from the
    //  reference element to the local element.
    double jacobian;
    jacobian = 2.0 / Delta;;

    // holds the values of the basis functions of order 0,1, .. ,
    // order at the point x
    double basis_values[order];
    // holds the values of the basis functions of order 0,1, .. ,
    // order at the point x
    double basis_deriv_values[order];
// 	std::cout << "xScaled = " << xScaled << std::endl;
    if(  (x > xRight) || (x < xLeft) )
    {
        //assert (fabs(xScaled) >= 1.0);
        return 0.0;
    }
    if(Delta == 0.0)
    {
        return 0;
    }
    else
    {
        // we need to cast the order into a double for the gsl to use it.
        double dorder = (double) order;
        gsl_sf_legendre_Pl_deriv_array(order, xScaled, basis_values, basis_deriv_values);

        switch (numDerivs)
        {
        case zero:
            temp1 = gsl_sf_legendre_Pl(order, xScaled); // gets value of basis function with order order
            break;
        case one:
            temp1 = basis_deriv_values[order]; // gets value of deriv of basis function
            temp1 *= jacobian;                 //with order order
            break;
        default:
            cout << "Can't happen on line number "  << __LINE__ << endl;
            break;
            return 1;
        }
        return (temp1);
    }

}



double
ddpLimPsi(ddpGrid_type const & grid,
          int const & elem,
          int const & order,
          ddpDirection_type const & Flag,
          int const & endPointNumber
         )
{
    // This returns the value of the elem-th DG basis function of order order
    // with numDerive number of derivatives taken as the limit going to point value
    // x from the left (flag1 = minus) OR right (flag1 = plus)


    // Note:  These fluxes are for the legendre polynomials.
    // If a different choice of basis function is chosen, this has to be modified
    // appropriately.

    int i = elem;
    int istar = endPointNumber;
    double temp;
    int numElementsNoGhost = grid.NumElementsNoGhost;
    if(Flag == Plus)
    {
        if( istar == i) // Then we are at the left endpoint of the i'th elemement
        {
            temp = pow( -1.0, (double) order);
            return temp;
        }
        else
            return 0.0;
    }
    else if(Flag == Minus)
    {
        if(i + 1 == istar )  // Then we are the right endpoint of the i'th element
        {
            return 1.0;
        }
        else
            return 0.0;
    }
    else
    {
        cout << "can't happen on line " << __LINE__ << endl;
        return NAN;
        assert(false);
    }
}



double
ddpUpsilon( ddpGrid_type const & grid,
            int const & elem,
            int const & order,
            ddpNumDeriv_type const & numDerivs,
            double const & x
          )
{
    // This returns the value of the elem-th CG Mixed FEM basis function of
    // order order with numDerive number of derivatives at the point x
    //

    // NOTE: This makes the basis hat/bubble functions from the DG basis
    // fucnctions.

    int NumElementsNoGhost = grid.NumElementsNoGhost;
    int J = (order+1) % 2;

    if( ( 0 == order ) && (0 == elem) )
    {
        return (0.5) *(ddpPsi(grid,0,0,numDerivs,x)
                       - ddpPsi(grid,0,1,numDerivs,x) );
    }
    else if(  (0 == order) && (elem == (NumElementsNoGhost)))
    {
        return (0.5) * (ddpPsi(grid,NumElementsNoGhost - 1, 0,numDerivs,x)
                        + ddpPsi(grid,NumElementsNoGhost - 1, 1,numDerivs,x) );
        //The minus two here is something we don't really like, but it is
        // on purpose.
    }
    else if (0 == order)
    {
        return (0.5) * ( ddpPsi(grid, elem - 1,0,numDerivs,x)
                         + ddpPsi(grid, elem - 1,1,numDerivs,x)
                         + ddpPsi(grid,elem, 0, numDerivs, x)
                         -ddpPsi(grid,elem,1,numDerivs,x) );
    }
    else
    {
        return ddpPsi(grid, elem, order+1, numDerivs, x)
               - ddpPsi(grid,elem, J,numDerivs, x);
    }

}
