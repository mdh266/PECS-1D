#ifndef __VECTORSETUP_H__
#define __VECTORSETUP_H__

#include"includes.hpp"


int ddpMakeGeneralFluxVector(ddpGrid_type const & grid,
                             ddpDenseVector_type const & VposOrVneg,
                             int const & ElementNumber,
                             ddpBijForw_type const & FuncForward,
                             ddpDGFluxMatrices_type  const & DGFluxMatrices,
                             ddpDirection_type const & Flag1,
                             ddpSparseMatrix_type & output);


int ddpMakeSpecialFluxVector(ddpGrid_type const & grid,
                             int const & ElementNumber,
                             ddpBijForw_type const & FuncForward,
                             ddpDGFluxMatrices_type const & DGFluxMatrices,
                             ddpDirection_type const & Flag1,
                             ddpSparseMatrix_type & output);


#endif
