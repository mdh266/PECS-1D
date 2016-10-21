#ifndef _TESTUTILITIES_H_
#define _TESTUTILITIES_H_

#include "includes.hpp"
#include "main.hpp"
#include "carrier.hpp"
#include "poisson.hpp"
#include "dopingprofile.hpp"


typedef enum {PositiveCharge, NegativeCharge} ddpTestChargeCarrier_type;

double MinusSineFunction(const double & x);

double MinusCosineFunction(const double & x);

double SineFunction(const double & x);

double Sine2Function(double const & x);

double Sine2PiFunction(double const & x);

double CosineFunction(const double & x);

double CosinePiFunction(const double & x);

double Cosine2PiFunction(const double & x);

double Cosine2Function(const double & x);

double PiFunction(const double & x);

double ThirtyFunction(const double & x);

double TwoFunction(const double & x);

double MinusPiFunction(const double & x);

double OneFunction(const double & x);

double MinusOneFunction(const double & x);

double OneOverPiLineFunction(double const & x);

double ExponentialFunction(const double & x);

double PotentialInterfaceZeroRHS(double const & x);

double DisplElecFieldInterfaceZeroRHS(double const & x);

double PotentialInterfaceNonZeroRHS1(double const & x);

double DisplElecFieldInterfaceNonZeroRHS1(double const & x);

double PotentialInterfaceNonZeroRHS2(double const & x);

double DisplElecFieldInterfaceNonZeroRHS2(double const & x);

double MinusOneOverPiFunction(double const & x);

double DD_L2_RHS(double const & x);

double DD_L2_Density(double const & x);

double  DD_L2_Current(double const & x);

double DD_L2_BC(double const & x);
double Poisson_L2_RHS(double const & x);

double Poisson_L2_Potential(double const & x);

double  Poisson_L2_ElecField(double const & x);

double Doping_N_A(double const & x);

double Doping_N_D(double const & x);

double Doping_PN(double const & x);

double True_N_A(double const & x);

bool Are_Same(ddpCarrier_type const & carrier1,
              ddpCarrier_type const & carrier2);

bool Are_Same_Current(ddpCarrier_type const & carrier1,
                      ddpCarrier_type const & carrier2);

bool Dofs_Are_Same(ddpPoisson_type const & P1, ddpPoisson_type const & P2);

bool Are_Same(ddpDenseVector_type const & DP1, ddpDenseVector_type const & DP2);

int L2_Error(ddpPoisson_type const & P1,
             double (*TrueU)(const double & x),
             double (*TrueQ)(const double & x),
             ddpDenseVector_type & errors);


int L2_Error(ddpCarrier_type const & TestP,
             double (*TrueU)(const double & x),
             double (*TrueQ)(const double & x),
             ddpDenseVector_type & errors);


int print2File(ddpCarrier_type const & testCarrier,
               ddpCarrier_type const & trueCarrier,
               ddpPoisson_type const & testPoisson,
               int testNum);

///////////////////////////////////////////////////////////////////////////////
// Diffusion Problem Functions
///////////////////////////////////////////////////////////////////////////////

// Dirichlet Homogenous BC
double DiffusionSolutionHomoU(double const & x);

double DiffusionSolutionHomoQ(double const & x);


// Dirichlet Non-Homogenous BC

double DiffusionNonHomoLeftBC(double const & t);

double DiffusionNonHomoRightBC(double const & t);

double DiffusionNonHomoLeftBC_Scaling(double const & t);

double DiffusionNonHomoRightBC_Scaling(double const & t);

double DiffusionSolutionNonHomoU(double const & x);

double DiffusionSolutionNonHomoQ(double const & x);

double DiffusionSolutionNonHomoU_Scaling(double const & x);

double DiffusionSolutionNonHomoQ_Scaling(double const & x);


// Neumann BC
double DiffusionSolutionNeumannHomoU(double const & x);

double DiffusionSolutionNeumannHomoQ(double const & x);

double DiffusionNeumannNonHomoBC(double const & t);

double DiffusionSolutionNeumannNonHomoU(double const & x);

double DiffusionSolutionNeumannNonHomoQ(double const & x);

double DiffusionNeumannHomoBC_Scaling(double const & t);

double DiffusionSolutionNeumannHomoU_Scaling(double const & x);

double DiffusionSolutionNeumannHomoQ_Scaling(double const & x);

double DiffusionWithScalingInitial(const double & x);

double DiffusionWithScalingAndNeumannRightBC(double const & t);

double DiffusionWithScalingU(double const & x);

double DiffusionWithScalingQ(double const & x);

///////////////////////////////////////////////////////////////////////////////
// Electron Drift Diffusion Problem Functions
///////////////////////////////////////////////////////////////////////////////

// Dirichlet Homogenous BC
double EDD_Dirichlet_Homo_IC(double const & x);

double EDD_Dirichlet_Homo_SolutionU(double const & x);

double EDD_Dirichlet_Homo_SolutionnQ(double const & x);

double EDD_Dirichlet_Homo_IC(double const & x);

// Dirichlet Non-Homogenous BC

double EDD_Dirichlet_Non_Homo_IC(double const & x);

double EDD_Dirichlet_Non_Homo_LeftBC(double const & t);

double EDD_Dirichlet_Non_Homo_RightBC(double const & t);

double EDD_Dirichlet_Non_Homo_SolutionU(double const & x);

double EDD_Dirichlet_Non_Homo_SolutionnQ(double const & x);


///////////////////////////////////////////////////////////////////////////////
// Holes Drift Diffusion Problem Functions
///////////////////////////////////////////////////////////////////////////////

// Dirichlet Homogenous BC
double HDD_Dirichlet_Homo_IC(double const & x);

double HDD_Dirichlet_Homo_SolutionU(double const & x);

double HDD_Dirichlet_Homo_SolutionnQ(double const & x);

double HDD_Dirichlet_Homo_IC(double const & x);

// Dirichlet Non-Homogenous BC

double HDD_Dirichlet_Non_Homo_IC(double const & x);

double HDD_Dirichlet_Non_Homo_LeftBC(double const & t);

double HDD_Dirichlet_Non_Homo_RightBC(double const & t);

double HDD_Dirichlet_Non_Homo_SolutionU(double const & x);

double HDD_Dirichlet_Non_Homo_SolutionnQ(double const & x);



// MIXED BOUNDARY CONDITIONS

// LEFT HOMOGENOUS DIRICHLET
double EDD_Dir_Homo_Robin_Non_Homo_RightBC(double const & t);

// LEFT NON-HOMOGENOUS DIRICHLET
double EDD_Dir_Non_Homo_Robin_Non_Homo_RightBC(double const & t);

// RIGHT HOMOGENOUS DIRICHLET
double EDD_Robin_Non_Homo_Dir_Homo_LeftBC(double const & t);

// NON HOMOGENOUS DIRICHLET
double EDD_Robin_Non_Homo_Dir_Non_Homo_LeftBC(double const & t);

// HOLES

// LEFT NEUMANN
double HDD_Non_Homo_Neumann_Homo_Dir_LeftBC(double const & t);

// LEFT ROBIN
double HDD_Non_Homo_Robin_Non_Homo_Dir_LeftBC(double const & t);

// RIGHT NEUMANN
double HDD_Homo_Dir_Non_Homo_Neumann_RightBC(double const & t);

// RIGHT ROBIN
double HDD_Homo_Dir_Non_Homo_Robin_RightBC(double const & t);

double HDD_Right_Robin_IC(double const & x);

double HDD_Right_Robin_Solution_U(double const & x);

double HDD_Right_Robin_Solution_Q(double const & x);

//////////////////////////////////////////////////////////////////////
// STEADY STATE SOLVER
/////////////////////////////////////////////////////////////////////
double HDD_SteadyState_Dir_U(double const & x);

double HDD_SteadyState_Dir_Q(double const & x);

double Diffusion_SteadyState_Robin_U(double const & x);

double HDD_SteadyState_Robin_U(double const & x);

////////////////////////////////////////////////////////////////////
// REACTIVE FLUXES WITH FORWARD EULER
///////////////////////////////////////////////////////////////////

double ReactiveFluxesU(double const & x);

double ReactiveFluxesQ(double const & x);


#endif
