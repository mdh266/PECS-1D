#ifndef __PROBLEMSETUP_H_
#define __PROBLEMSETUP_H_

#include "includes.hpp"

///////////////////////////////////////////////////////////////////////////////
// In problemSetUp.cpp
///////////////////////////////////////////////////////////////////////////////

int
ddpMakeUniformGrid(ddpProblemInfo_type const & problem,
                   ddpDomain_type const & domain,
                   ddpGrid_type & grid);
int
ddpMakeLeftRefinedBoundaryGrid(ddpProblemInfo_type const & problem,
                               ddpDomain_type const & domain,
                               ddpGrid_type & grid);
int
ddpMakeRightRefinedBoundaryGrid(ddpProblemInfo_type const & problem,
                                ddpDomain_type const & domain,
                                ddpGrid_type & grid);


int
ddpMakeBothEndsRefinedGrid(ddpProblemInfo_type const & problem,
                           ddpDomain_type const & domain,
                           ddpGrid_type & grid);

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
                      const int & gaussLegendreNumPoints);

int
ddpGlueGrids(ddpProblemInfo_type const & problem,
             ddpGrid_type const & left_grid,
             ddpGrid_type const & right_grid,
             ddpGrid_type & GluedGrid);

int
ddpMakeTimeStamps(ddpProblemInfo_type const & problem,
                  std::vector< double > & timeStamps);

int
ddpMakeBijs(ddpGrid_type const & grid,
            ddpBijFlag_type const & flag,
            ddpBijForw_type& Forward,
            ddpBijBack_type& Backward);

int
ddpPrintBij_Backward(ddpBijBack_type const & bij);


int
ddpMakeWeightVectorDense(ddpGrid_type const & grid,
                         ddpDenseVector_type & weights);

int
ddpMakeWeightVectorSparse(ddpGrid_type const & grid,
                          ddpSparseVector_type & weights);
int
ddpMakePTVectorSparse(ddpGrid_type const & grid,
                      ddpSparseVector_type & points);
int
ddpMakePTVectorDense(ddpGrid_type const & grid,
                     ddpDenseVector_type & points);
int
ddpMakeWeightsAndPoints(ddpGrid_type const & grid,
                        ddpProblemInfo_type const & problem,
                        ddpSparseVector_type & weightsSparse,
                        ddpSparseVector_type & PTSSparse,
                        ddpDenseVector_type & PTSDense,
                        ddpDenseVector_type & weightsDense);

int
ddpMakeAllBijections(ddpGrid_type const & grid,
                     ddpBijection_type & Bijections);



#endif
