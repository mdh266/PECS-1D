#ifndef __MATRIXSETUP_H_
#define __MATRIXSETUP_H_

#include "includes.hpp"


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
                    );

int
ddpMakeGeneralFluxMatrix(ddpGrid_type const & grid,
                         ddpDenseVector_type const & Vpos,
                         ddpDenseVector_type const & Vneg,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         ddpSparseMatrix_type & output);

int
ddpMakeSpecialFluxMatrix(ddpGrid_type const & grid,
                         ddpBijForw_type const & FuncForward,
                         ddpDGFluxMatrices_type const & DGFluxMatrices,
                         ddpDirection_type const & Flag1,
                         std::vector<ddpDirection_type> const& Flag2Vector,
                         ddpSparseMatrix_type & output);

int
ddpMakeInteriorFacesFluxMatrix(ddpGrid_type const & grid,
                               ddpBijForw_type const & FuncForward,
                               ddpDGFluxMatrices_type const & DGFluxMatrices,
                               ddpDirection_type const & Flag1,
                               std::vector<ddpDirection_type> const& Flag2Vector,
                               ddpSparseMatrix_type & output);


int
ddpMakeBoundaryFacesFluxMatrix(ddpGrid_type const & grid,
                               ddpBijForw_type const & FuncForward,
                               ddpDGFluxMatrices_type  const & DGFluxMatrices,
                               ddpDirection_type const & Flag1,
                               std::vector<ddpDirection_type> const& Flag2Vector,
                               ddpSparseMatrix_type & output);

#endif
