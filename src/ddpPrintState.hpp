#ifndef _DDPPRINTSTATE_H_
#define _DDPPRINTSTATE_H_

#include "carrier.hpp"
#include "poisson.hpp"
#include <iomanip>


int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type const & carrier1,
              ddpPoisson_type const & poisson,
              int const & timeStamp);

int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type  & electrons,
              ddpCarrier_type  & holes,
              ddpPoisson_type const & poisson,
              int const & timeStamp,
              double const & tCurrent
             );

int
ddpPrintStateTest(ddpProblemInfo_type const & problem,
                  ddpCarrier_type  & carrier,
                  ddpCarrier_type  & carrier1,
                  ddpCarrier_type  & carrier2,
                  ddpPoisson_type const & poisson,
                  int const & timeStamp
                 );

int
ddpPrintState(ddpProblemInfo_type const & problem,
              ddpCarrier_type const & electrons,
              ddpCarrier_type const & holes,
              ddpCarrier_type const & reductants,
              ddpCarrier_type const & oxidants,
              ddpPoisson_type const & poisson,
              int const & timeStamp,
              double const & tCurrent);
int
ddpPrintFinalDofs(ddpCarrier_type const & electrons,
                  ddpCarrier_type const & holes,
                  ddpCarrier_type const & reductants,
                  ddpCarrier_type const & oxidants,
                  ddpPoisson_type const & poisson);

int
ddpPrintFinalDofs(ddpCarrier_type const & electrons,
                  ddpCarrier_type const & holes,
                  ddpPoisson_type const & poisson);

int ddpReadInElectronDOFS(ddpDenseVector_type & electrons_u,
                          ddpDenseVector_type & electrons_q);
int
ddpReadInFinalStates(ddpCarrier_type & electrons,
                     ddpCarrier_type & holes,
                     ddpPoisson_type & poisson);

int
ddpReadInFinalStates(ddpCarrier_type & electrons,
                     ddpCarrier_type & holes,
                     ddpCarrier_type & reductants,
                     ddpCarrier_type & oxidants,
                     ddpPoisson_type & poisson);

int
progressBar(const int & current, const int & total,
            const int & width);

#endif
