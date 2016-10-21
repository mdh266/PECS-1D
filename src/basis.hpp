#ifndef __BASIS_H_
#define __BASIS_H_
#include "includes.hpp"

/** \brief  DG Basis function. */
/** DG basis function is taken to be Legrendre polynomial.
 * NOTE: This uses the gsl/gsl_legrendre_sf function to get
 * values at the quadrature points. */
double
ddpPsi( ddpGrid_type const & grid,
        int const & elem,
        int const & order,
        ddpNumDeriv_type const & numDerivs,
        double const & x);


/** \brief Function used in takeing limits of DG Basis functiom. */
double
ddpLimPsi(ddpGrid_type const & grid,
          int const & elem,
          int const & order,
          ddpDirection_type const & Flag,
          int const & endPointNumber);

/** \brief MX CG FEM Basis function. */
/** Note these are continuous basis functions, they are the normal
 * hat functions for FEM that live in H^{1}.  They are built up
 * from the DG basis functions. */
double
ddpUpsilon( ddpGrid_type const & grid,
            int const & elem,
            int const & order,
            ddpNumDeriv_type const & numDerivs,
            double const & x);
#endif
