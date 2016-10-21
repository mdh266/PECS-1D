#ifndef __VANDEMONDESETUP_H__
#define __VANDEMONDESETUP_H__


#include "includes.hpp"
#include "basis.hpp"

int
ddpMakeGlobalVandeMonde(ddpGrid_type const & grid,
                        ddpProblemInfo_type const & problem,
                        ddpBijFlag_type const & flag,
                        ddpBijForw_type const & FuncForward,
                        ddpBijForw_type const & PTForward,
                        ddpNumDeriv_type const & numDeriv,
                        ddpSparseMatrix_type & vandeMonde);

int
ddpMakeGlobalVandeMondeFluxMX(ddpGrid_type const & grid,
                              ddpProblemInfo_type const & problem,
                              ddpBijForw_type const & FuncForward,
                              ddpSparseMatrix_type & vandeMonde);

int
ddpMakeGlobalVandeMondeFluxDG(ddpGrid_type const & grid,
                              ddpProblemInfo_type const & problem,
                              ddpBijForw_type const & FuncForward,
                              ddpDirection_type const & Flag1,
                              ddpDirection_type const & Flag2,
                              ddpSparseMatrix_type & vandeMonde);


int
ddpMakeVandeMondeMatrices(ddpGrid_type const & grid,
                          ddpProblemInfo_type const & problem,
                          ddpBijection_type const & Bijections,
                          ddpVandeMondeMatrices_type & VandeMondeMatrices,
                          ddpDGFluxMatrices_type & DGFluxMatrices);

#endif
