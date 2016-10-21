#include "../include/dopingprofile.hpp"


int
ddpDopingProfile_type::
setDopingProfile(ddpGrid_type const & grid,
                 ddpCarrier_type const & carrier,
                 double (* DopingProfileFunction) (double const & x) 	)

{

    ddpBijFlag_type BijFlag = DG;

    // set C(x) = N_D(x) or C(x) = -N_A(x)
    if(carrier.carrierProps.ChargeSign == Negative)
    {

        // set donor function
        this->setDonorFunction(DopingProfileFunction);

        // Project doping function onto basis
        ddpProjectFunction(DonorDopingProfileFunction,
                           grid,
                           carrier.Bijections.DGForward,
                           carrier.VandeMondeMatrices.globalVandeMondeDG,
                           carrier.Bijections.PTForward,
                           BijFlag,
                           carrier.weightsSparse,
                           carrier.carrierProps.MassU,
                           Dof);
        // C(x) = N_D

    }

    else
    {
        // set acceptor function
        this->setAcceptorFunction(DopingProfileFunction);

        // Project doping function onto basis
        ddpProjectFunction( AcceptorDopingProfileFunction,
                            grid,
                            carrier.Bijections.DGForward,
                            carrier.VandeMondeMatrices.globalVandeMondeDG,
                            carrier.Bijections.PTForward,
                            BijFlag,
                            carrier.weightsSparse,
                            carrier.carrierProps.MassU,
                            Dof);

        // switch signs so C(x) = -N_A
        Dof = -1.0*Dof;
    }


    // Set the point value array
    PTValue = (carrier.VandeMondeMatrices.globalVandeMondeDG)
              * Dof;

    return 0;
}

// Projects the doping function on to the basis and
// stores it in this structures Dof
int
ddpDopingProfile_type::
setDopingProfile(ddpGrid_type const & grid,
                 ddpCarrier_type const & carrier1,
                 ddpCarrier_type const & carrier2,
                 double (* DonorProfileFunction) (double const & x),
                 double (* AcceptorProfileFunction) (double const & x)
                )
{
    ddpBijFlag_type BijFlag = DG;

    // set donor function
    this->setDonorFunction(DonorProfileFunction);

    // set acceptor function
    this->setAcceptorFunction(AcceptorProfileFunction);

    ddpProjectFunction(DonorDopingProfileFunction,
                       grid,
                       carrier1.Bijections.DGForward,
                       carrier1.VandeMondeMatrices.globalVandeMondeDG,
                       carrier1.Bijections.PTForward,
                       BijFlag,
                       carrier1.weightsSparse,
                       carrier1.carrierProps.MassU,
                       DonorDof);

    // carrier2 = type for N_A
    ddpProjectFunction(AcceptorDopingProfileFunction,
                       grid,
                       carrier2.Bijections.DGForward,
                       carrier2.VandeMondeMatrices.globalVandeMondeDG,
                       carrier2.Bijections.PTForward,
                       BijFlag,
                       carrier2.weightsSparse,
                       carrier2.carrierProps.MassU,
                       AcceptorDof);

    // Set C(x) = N_D - N_A
    Dof = DonorDof - AcceptorDof;

    // Set the point value array
    // TODO:  Assumes carrier1 and carrier2 have same vandeMonde matrix
    PTValue = (carrier1.VandeMondeMatrices.globalVandeMondeDG)
              * Dof;

    return 0;
}

// set the donor doping function
int
ddpDopingProfile_type::
setDonorFunction(double (* DopingProfileFunction) (double const & x))
{
    DonorDopingProfileFunction = DopingProfileFunction;
    return 0;
}


// set the acceptor doping function
int
ddpDopingProfile_type::
setAcceptorFunction(double (* DopingProfileFunction) (double const & x))
{
    AcceptorDopingProfileFunction = DopingProfileFunction;
    return 0;
}



