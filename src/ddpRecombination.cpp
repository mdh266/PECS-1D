#include "../include/ddpRecombination.hpp"

ddpRecombination_type::
ddpRecombination_type(ddpCarrier_type & carrier1,
                      ddpCarrier_type & carrier2)
{

    size = carrier1.weightsDense.size();

    IntrinsicDensity = carrier1.carrierProps.IntrisincDensity;
    RecombinationTime1 = carrier1.carrierProps.RecombinationTime;
    RecombinationTime2 = carrier2.carrierProps.RecombinationTime;

    BijFlag = DG;

    RecombinationPTVals = ddpDenseVector_type::Zero(size);
}

ddpRecombination_type::
~ddpRecombination_type()
{
}

int
ddpRecombination_type::
updateRecombination(ddpGrid_type const & grid,
                    ddpCarrier_type & carrier1,
                    ddpCarrier_type & carrier2)
{

    // Get the carrier density point values
    Carrier1PTVals = carrier1.VandeMondeMatrices.globalVandeMondeDG *
                     carrier1.carrierState.uDof;

    Carrier2PTVals = carrier2.VandeMondeMatrices.globalVandeMondeDG *
                     carrier2.carrierState.uDof;


    // make the recombination point values
    #pragma omp parallel for
    for(int i = 0; i < size; i++)
    {
        RecombinationPTVals(i) =
            ( (IntrinsicDensity*IntrinsicDensity) - (Carrier1PTVals(i)*Carrier2PTVals(i)) )
            /
            ( (RecombinationTime1*(Carrier1PTVals(i) + IntrinsicDensity) ) +
              (RecombinationTime2*(Carrier2PTVals(i) + IntrinsicDensity) )
            );
    }

    std::pair<int, int> elemAndOrderJ;
    int globalJ;
    ddpSparseVector_type column;
    double Threepsilon = 12.0*std::numeric_limits<double>::epsilon();

    int elem;
    int order;
    double partial_sum;
    int numElem = grid.NumElementsNoGhost;
    #pragma omp parallel for private(order, partial_sum, globalJ, elemAndOrderJ, column)
    for(elem=0; elem<numElem; elem++)
    {
        for(order=0; order<=grid.ElementList[elem].OrderDG; order++)
        {
            elemAndOrderJ.first = elem;
            elemAndOrderJ.second = order;
            globalJ =
                carrier1.Bijections.DGForward.find(elemAndOrderJ)->second;
            column  =
                carrier1.VandeMondeMatrices.globalVandeMondeDG.col(globalJ);
            partial_sum = 0.0;
            for(ddpSparseVector_type::InnerIterator it(column); it; ++it)
            {
                partial_sum += column.coeffRef(it.index()) *
                               carrier1.weightsSparse.coeffRef(it.index()) *
                               RecombinationPTVals.coeffRef(it.index());
            }
            if(fabs(partial_sum) > Threepsilon)
                carrier1.carrierState.RHSFromRecombination(globalJ) = partial_sum;
            else
                carrier1.carrierState.RHSFromRecombination(globalJ) = 0.0;
        }
    }

    /*
    	  // project the recombination pt values onto the basis fucntions
    	  // \int_{I} v R(rho_n, rho_p) dx
    		ddpProjectRHSFunction(RecombinationPTVals,
    											 grid,
    											 carrier1.Bijections.DGForward,
    											 carrier1.VandeMondeMatrices.globalVandeMondeDG,
    											 carrier1.Bijections.PTForward,
    											 BijFlag,
    											 carrier1.weightsSparse,
    											 carrier1.carrierState.RHSFromRecombination);
     */
    // copy this vector into RHS vector for carrier2
    carrier2.carrierState.RHSFromRecombination =
        carrier1.carrierState.RHSFromRecombination;

    return 0;
}
