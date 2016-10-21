#ifndef _DOPINGPROFILE_H_
#define _DOPINGPROFILE_H_

#include "main.hpp" // for stuff in SimulationConditions.cpp
#include "includes.hpp"
#include "carrier.hpp"
//#include "poisson.hpp"

typedef struct ddpDopingProfile_type
{

///////////////////////////////////////////////////////////////////////////////
//   Data Members
///////////////////////////////////////////////////////////////////////////////
    // Dofs
    ddpDenseVector_type Dof;
    ddpDenseVector_type PTValue;

    ddpDenseVector_type DonorDof;
    ddpDenseVector_type AcceptorDof;



    // Doping Functions
    double (* DonorDopingProfileFunction)(double const & x);
    double (* AcceptorDopingProfileFunction)(double const & x);


///////////////////////////////////////////////////////////////////////////////
//   Class Methods
///////////////////////////////////////////////////////////////////////////////

    // Projects the doping function on to the basis and
    // stores it in this structures Dof

    //Unipolar Model
    int setDopingProfile(ddpGrid_type const & grid,
                         ddpCarrier_type const & carrier,
                         double (* DopingProfileFunction) (double const & x	)  );

    // bipolar model
    int setDopingProfile(ddpGrid_type const & grid,
                         ddpCarrier_type const & carrier1,
                         ddpCarrier_type const & carrier2,
                         double (* DonorProfileFunction) (double const & x),
                         double (* AcceptorProfileFunction) (double const & x)  );

    // set the donor doping function
    int setDonorFunction(double (* DopingProfileFunction) (double const & x));

    // set the acceptor doping function
    int setAcceptorFunction(double (* DopingProfileFunction) (double const & x));


} ddpDopingProfile_type;

#endif
