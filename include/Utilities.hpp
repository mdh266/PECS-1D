#ifndef __UTILITIES_H_
#define __UTILITIES_H_

#include "includes.hpp"

int
ddpConcatenateMatrices(ddpSparseMatrix_type const & A,
                       ddpSparseMatrix_type const & B,
                       ddpSparseMatrix_type & C);

int
ddpPositivePart(ddpDenseVector_type const & input,
                ddpDenseVector_type & output);

int
ddpNegativePart(ddpDenseVector_type const & input,
                ddpDenseVector_type & output);


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
    ddpDenseVector_type & output);

int
ddpProjectFunction( double (*function)(double const &),
                    ddpGrid_type const & grid,
                    ddpBijForw_type const & FuncForward,
                    ddpSparseMatrix_type const & globalVandeMondeJ,
                    ddpBijForw_type const & PTForward,
                    ddpBijFlag_type const & BijFlag,
                    ddpSparseVector_type const & sparseWeights,
                    ddpSparseMatrix_type const & MassMatrix,
                    ddpDenseVector_type & output);
int
ddpProjectFunction( ddpDenseVector_type const & FunctionPointVals,
                    ddpGrid_type const & grid,
                    ddpBijForw_type const & FuncForward,
                    ddpSparseMatrix_type const & globalVandeMondeJ,
                    ddpBijForw_type const & PTForward,
                    ddpBijFlag_type const & BijFlag,
                    ddpSparseVector_type const & sparseWeights,
                    ddpSparseMatrix_type const & MassMatrix,
                    ddpDenseVector_type & output);

int
ddpProjectFunction( ddpDenseVector_type & FunctionPointVals,
                    ddpGrid_type const & grid,
                    ddpBijForw_type & FuncForward,
                    ddpSparseMatrix_type & globalVandeMondeJ,
                    ddpBijForw_type & PTForward,
                    ddpBijFlag_type & BijFlag,
                    ddpSparseVector_type & sparseWeights,
                    ddpSparseMatrix_type & MassMatrix,
                    ddpDenseVector_type & output);

double ZeroFunction(const double & t);

int
ddpProjectRHSFunction( double (*function)(double const &),
                       ddpGrid_type const & grid,
                       ddpBijForw_type const & FuncForward,
                       ddpSparseMatrix_type const & globalVandeMondeJ,
                       ddpBijForw_type const & PTForward,
                       ddpBijFlag_type const & BijFlag,
                       ddpSparseVector_type const & sparseWeights,
                       ddpDenseVector_type & output);

int
ddpProjectRHSFunction( ddpDenseVector_type const & FunctionPointVals,
                       ddpGrid_type const & grid,
                       ddpBijForw_type const & FuncForward,
                       ddpSparseMatrix_type const & globalVandeMondeJ,
                       ddpBijForw_type const & PTForward,
                       ddpBijFlag_type const & BijFlag,
                       ddpSparseVector_type const & sparseWeights,
                       ddpDenseVector_type & output);

#endif
