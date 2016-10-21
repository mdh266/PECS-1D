#include "../include/SimulationConditions.hpp"

/////////////////////////////////////////////////////////////////////////
// Devices!  These are their doping profiles and
/////////////////////////////////////////////////////////////////////////

// These are interms of the characteristic denisty

/////////////////////////////////////////////////////////////////////
// Testing Ones
/////////////////////////////////////////////////////////////////////

double
StepDoping(double const & x)
{
    if(x < 0.5)
        return 1.0; // - (x)*(x);
    else
        return 0.0;// - (x-1)*(x-1) ;
}

double
Normal(double const & x)
{
    double mean = 0.5;
    double sigma = 0.1;

    double result = (1.0/(sigma * sqrt(2.0*M_PI)));
    return 1.0 + result * exp(-(x-mean)*(x-mean)/(2.0*sigma*sigma));
}



/////////////////////////////////////////////////////////////////////
// n+ - n - n+ device
/////////////////////////////////////////////////////////////////////

double
NPlus_N_NPlus(double const & x)
{

    if(x <= 0.3)
        return 1.0;
    else if ( (x > 0.3) && (x < 0.7) )
        return 1.0e-4;
    else
        return 1.0;

}




/////////////////////////////////////////////////////////////////////
// n-type Photoelectrochemical Cell
/////////////////////////////////////////////////////////////////////

double
ddpPEC_DonorDopingProfile(double const & x)
{
//	if(x < -0.5)
    return 2.0;
//	else
//		return 0.0;
}

double
ddpPEC_AcceptorDopingProfile(double const & x)
{
//	if(x < -0.5)
//		return 0.0;
//	else
//		return 2.0;
    return 0.0;
}

double
ddpPEC_ReductantsIC(double const & x)
{
    return 30.0;
}

double
ddpPEC_OxidantsIC(double const & x)
{
    return 29.0;
}


