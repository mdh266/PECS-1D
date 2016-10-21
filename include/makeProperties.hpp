#ifndef __MAKEPROPERTIES_H_
#define __MAKEPROPERTIES_H_

#include "includes.hpp"
#include "matrixSetUp.hpp"
#include "vectorSetUp.hpp"
#include "Utilities.hpp"


int
ddpMakeCarrierProperies(
    ddpProblemInfo_type const & problem,
    ddpGrid_type const & grid,
    ddpSparseVector_type const & weightsSparse,
    ddpBijection_type const & Bijections,
    ddpVandeMondeMatrices_type const & VandeMondeMatrices,
    ddpDGFluxMatrices_type const & DGFluxMatrices,
    ddpCarrierProperties_type & CarrierProps);
int
ddpMakeDiffusiveFluxProperties(
    ddpProblemInfo_type const & problem,
    ddpGrid_type const & grid,
    ddpSparseVector_type const & weightsSparse,
    ddpBijection_type const & Bijections,
    ddpVandeMondeMatrices_type const & VandeMondeMatrices,
    ddpDGFluxMatrices_type const & DGFluxMatrices,
    ddpCarrierProperties_type & CarrierProps);

int
ddpMakePoissonProperties( ddpProblemInfo_type const & problem,
                          ddpGrid_type const & grid,
                          ddpSparseVector_type const & weightsSparse,
                          ddpBijection_type const & Bijections,
                          ddpVandeMondeMatrices_type const & VandeMondeMatrices,
                          ddpDGFluxMatrices_type const & DGFluxMatrices,
                          ddpPoissonProperties_type & PoissonProperties);

int
ddpMakePoissonPropertiesWithInterface(ddpProblemInfo_type const & problem,
                                      ddpGrid_type const & grid,
                                      ddpSparseVector_type & weightsSparse,
                                      ddpBijection_type const & Bijections,
                                      ddpVandeMondeMatrices_type const & VandeMondeMatrices,
                                      ddpDGFluxMatrices_type const & DGFluxMatrices,
                                      ddpPoissonProperties_type & PoissonProperties);


#endif
