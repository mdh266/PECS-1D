#ifndef _DDPRECOMBINATION_H__
#define  _DDPRECOMBINATION_H__
#include "includes.hpp"
#include "carrier.hpp"

typedef class ddpRecombination_type
{
private:
    int size;

    double IntrinsicDensity;
    double RecombinationTime1;
    double RecombinationTime2;

    ddpBijFlag_type BijFlag;

    ddpDenseVector_type RecombinationPTVals;
    ddpDenseVector_type	Carrier1PTVals;
    ddpDenseVector_type	Carrier2PTVals;

public:
    ddpRecombination_type(ddpCarrier_type & carrier1,
                          ddpCarrier_type & carrier2);

    ~ddpRecombination_type();


    int updateRecombination(ddpGrid_type const & grid,
                            ddpCarrier_type & carrier1,
                            ddpCarrier_type & carrier2);

} ddpRecombination_type;
#endif
