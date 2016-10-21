#ifndef _ddpComputeDeltaT_H_
#define _ddpComputeDeltaT_H_

#include"carrier.hpp"

// unipolar model
int
ddpComputeDeltaT(ddpGrid_type const & grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier,
                 double & DeltaT
                );

// bipolar model
int
ddpComputeDeltaT(ddpGrid_type const & grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier1,
                 ddpCarrier_type const & carrier2,
                 double & DeltaT
                );

// interface model
int
ddpComputeDeltaT(ddpGrid_type const & semiconductor_grid,
                 ddpGrid_type const & electrolyte_grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier1,
                 ddpCarrier_type const & carrier2,
                 ddpCarrier_type const & carrier3,
                 ddpCarrier_type const & carrier4,
                 double & DeltaT
                );
#endif
