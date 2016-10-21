#ifndef _TEST_H_
#define _TEST_H_


#include "main.hpp"
#include "ddpPrintState.hpp"
#include "ddpComputeDeltaT.hpp"
#include "testUtilities.hpp"
#include "carrier.hpp"
#include "poisson.hpp"
#include "dopingprofile.hpp"

//
// This is the class for testing.  Each test will be a instantiation of
// this class.



typedef struct ddpTest_type
{

    ddpPoisson_type testPoisson;
    ddpDopingProfile_type DopingProfile;

    // for testing Poisson
    ddpPoisson_type truePoisson;
    ddpCarrier_type testCarrier1;
    ddpCarrier_type testCarrier2;


    double (* Carrier1ICFunction) (const double & x);
    double (* Carrier2ICFunction) (const double & x);



    double (* TrueSolutionPotential) (const double & x);
    double (* TrueSolutionElecField) (const double & x);

    double (* DopingProfileFunction1) (const double & x);
    double (* DopingProfileFunction2) (const double & x);

    double (* TrueDopingProfileFunction) (const double & x);


    ddpDenseVector_type TrueDopingDof;

    // for testing carriers
    ddpCarrier_type testCarrier;
    ddpCarrier_type trueCarrier;


    double (* TrueSolutionU) (const double & x);
    double (* TrueSolutionQ) (const double & x);


    // Carrier Initial And Boundary Conditions
    double (* CarrierICFunction) (const double & x);

    ddpBCType_type CarrierLeftBCType;
    ddpBCType_type CarrierRightBCType;

    ddpBCType_type Carrier1LeftBCType;
    ddpBCType_type Carrier1RightBCType;

    ddpBCType_type Carrier2LeftBCType;
    ddpBCType_type Carrier2RightBCType;

    double (* CarrierLeftBCFunction) (const double & x);
    double (* CarrierRightBCFunction) (const double & x);

    double (* Carrier1LeftBCFunction) (const double & x);
    double (* Carrier1RightBCFunction) (const double & x);

    double (* Carrier2LeftBCFunction) (const double & x);
    double (* Carrier2RightBCFunction) (const double & x);

    // Potential BC Functions

    double (* PotentialLeftBCFunction) (const double & x);
    double (* PotentialRightBCFunction) (const double & x);

    double (* DopingProfileFunction) (const double & x);



    ddpGrid_type testgrid;
    ddpGrid_type testgrid2;
    ddpProblemInfo_type testproblem;


    // Things to use to project funtion onto a basis
    ddpBijForw_type DGForward;
    ddpBijForw_type MXForward;

    ddpSparseMatrix_type globalVandeMondeDG;
    ddpSparseMatrix_type globalVandeMondeMX;

    ddpBijForw_type PTForward;

    ddpSparseVector_type sparseWeights;

    ddpSparseMatrix_type MassU;
    ddpSparseMatrix_type A00;

    ddpBijFlag_type BIJFlagDG;
    ddpBijFlag_type BIJFlagMX;

    ddpDenseVector_type EPTS;
    ddpDenseVector_type ElecFieldDof;

    ddpDenseVector_type EpsilonPTValues;
    ddpDenseVector_type EpsilonInvPTValues;

//////////////////////////////////////////////////////////////////////////
// Setting / Running Tests
/////////////////////////////////////////////////////////////////////////


    int initializeTest(const ddpTestChargeCarrier_type & carrierType,
                       const ddpGrid_type & inputgrid,
                       const ddpProblemInfo_type & inputproblem,
                       const ddpCarrierConstants_type & testConstants);

    int ReinitializeTestWithScaling(const ddpTestChargeCarrier_type & carrierType,
                                    const ddpGrid_type & inputgrid,
                                    const ddpProblemInfo_type & inputproblem,
                                    const ddpCarrierConstants_type & testConstants);

    // for poisson
    int ReinitializeTestWithInterface(const ddpGrid_type & left_grid,
                                      const ddpGrid_type & right_grid,
                                      const ddpProblemInfo_type & inputproblem,
                                      const ddpCarrierConstants_type & testConstants);

    // for carriers
    int initializeReactiveFluxTest(const ddpGrid_type & left_grid,
                                   const ddpGrid_type & right_grid,
                                   const ddpProblemInfo_type & inputproblem,
                                   const ddpCarrierConstants_type & testConstants);

    int setTestDoping();

    bool runTestDoping();

    int setTestPoisson();

    int setImplicitTestCarrier();

    int setSteadyStateTestCarrier();

    int setTestReactiveFluxes();

    bool runTestPoisson();

    void runTestPoissonL2Error(ddpDenseVector_type & errors);

    bool runTestPoissonWithInterface();

    bool runImplicitTestCarrier();

    bool runImplicitTestCarrierWithScaling();

    bool runSteadyStateTestCarrier();

    bool runSteadyStateShurTestCarrier();

    void runTestCarrierL2Error(ddpDenseVector_type & errors);

    void runImplicitTestCarrierL2Error(ddpDenseVector_type & errors);

    bool runImplicitTestReactiveFluxes();

    bool runImplicitExplicitTestReactiveFluxes();

    int outputTestCarrier(int testNumber);

    int outputTestPoisson(int testNumber);

    int	outputTestPoissonWithInterface(int testNumber);

} ddpTest_type;





#endif
