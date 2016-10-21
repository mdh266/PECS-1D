#include "../include/ddpComputeDeltaT.hpp"

///////////////////////////////////////////////////////////////////////////////
// Compute DeltaT for Unipolar Model
///////////////////////////////////////////////////////////////////////////////

// Computes DeltaT's for each carrier such that LDG
// scheme will be stable with forward Euler

int
ddpComputeDeltaT(ddpGrid_type const & grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier,
                 double & DeltaT
                )
{
    double maxV = electricFieldPTS.lpNorm<Eigen::Infinity>();

    double
    theta = (1.0 + grid.OrderDGMax) * (1.0 + grid.OrderDGMax) /
            (0.0 + grid.DeltaxMin);

    /*
     	cout << "Theta = " << theta << endl;
      cout << "DeltaXMin = " << grid.DeltaxMin << endl;
      cout << "OrderDG = " << grid.OrderDGMax<< endl;
    */
    // Note:  This uses the max and min of the diffusivity and mobility
    // for electron or reductant
    double numerator = 2.0 * carrier.carrierProps.MaxDiffusivity;

    double sqrtOfDenominator =  carrier.carrierProps.MaxMobility *maxV + 2 *
                                carrier.carrierProps.MaxDiffusivity * theta;

    double denominator = sqrtOfDenominator * sqrtOfDenominator;
    DeltaT = numerator / denominator;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Compute DeltaT for Bipolar Model
///////////////////////////////////////////////////////////////////////////////

// Takes the minimum of the DeltaT's for each carrier such that LDG
// scheme will be stable with forward Euler

int
ddpComputeDeltaT(ddpGrid_type const & grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier1,
                 ddpCarrier_type const & carrier2,
                 double & DeltaT
                )
{
    double maxV = electricFieldPTS.lpNorm<Eigen::Infinity>();

    double
    theta = (1.0 + grid.OrderDGMax) * (1.0 + grid.OrderDGMax) /
            (0.0 + grid.DeltaxMin);

    // Note:  This uses the max and min of the diffusivity and mobility
    // for electron or reductant and hole or oxidant

    double numerator1 = 2.0 * carrier1.carrierProps.MaxDiffusivity;
    double numerator2 = 2.0 * carrier2.carrierProps.MaxDiffusivity;

    double sqrtOfDenominator1 =  carrier1.carrierProps.MaxMobility *maxV + 2 *
                                 carrier1.carrierProps.MaxDiffusivity * theta;

    double sqrtOfDenominator2 =  carrier2.carrierProps.MaxMobility *maxV + 2 *
                                 carrier2.carrierProps.MaxDiffusivity * theta;

    double denominator1 = sqrtOfDenominator1 * sqrtOfDenominator1;
    double denominator2 = sqrtOfDenominator2 * sqrtOfDenominator2;

    double DeltaT1 = numerator1 / denominator1;
    double DeltaT2 = numerator2 / denominator2;

    DeltaT =  std::min(DeltaT1, DeltaT2);

    return 0;
}

int
ddpComputeDeltaT(ddpGrid_type const & semiconductor_grid,
                 ddpGrid_type const & electrolyte_grid,
                 ddpDenseVector_type const & electricFieldPTS,
                 ddpCarrier_type const & carrier1,
                 ddpCarrier_type const & carrier2,
                 ddpCarrier_type const & carrier3,
                 ddpCarrier_type const & carrier4,
                 double & DeltaT
                )
{
    // Computes DeltaT for the whole system as the minimum value for the delteT's
    // of each system.
    double DeltaT1, DeltaT2;

    ddpComputeDeltaT(semiconductor_grid,
                     electricFieldPTS,
                     carrier1,
                     carrier2,
                     DeltaT1);


    ddpComputeDeltaT(electrolyte_grid,
                     electricFieldPTS,
                     carrier3,
                     carrier4,
                     DeltaT2);

//	cout << "Semiconductor dt = " << DeltaT1 << endl;
//	cout << "Electrolyte dt = " << DeltaT2 << endl;

    DeltaT = std::min(DeltaT1, DeltaT2);

    return 0;
}


