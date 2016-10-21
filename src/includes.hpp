// includes.hpp
// Contains included libraries, using directives, structures, and typedefs
//
///////////////////////////////////////////////////////////////////////////////


#ifndef _INCLUDES_H_
#define _INCLUDES_H_
///////////////////////////////////////////////////////////////////////////////
// INCLUDES
///////////////////////////////////////////////////////////////////////////////

//#define NDEBUG

#include<iostream>
#include<fstream>

// Import functions from the C++ STL.
#include<map>
#include<vector>
#define _USE_MATH_DEFINES
#include<math.h>
#include<limits>
#include<assert.h>
#include<cstdlib>
#include<string>
#include<sstream>
#include<algorithm>


// disable annoying warnings from the intel compiler
//#pragma warning disable 869 383 981 2196 279 2536

// Import functions from the Intel MKL.

// Import functions from Gnu scientific library.
#include<gsl/gsl_integration.h>
#include<gsl/gsl_sf_legendre.h>
#include<gsl/gsl_math.h>

/////////////////////////////////////////////////////////////////////////////
// GRVY
////////////////////////////////////////////////////////////////////////////

//#include<grvy.h>
#include<sys/time.h>
#include<time.h>
//#include<hdf5.h>

//////////////////////////////////////////////////////////////////////////////
// TO USE UMFPACK uncomment and switch to
/////////////////////////////////////////////////////////////////////////////
/*
#include<mkl.h>

extern "C" {
#include <umfpack.h>
#include <amd.h>
#include <SuiteSparse_config.h>
}
#include <Eigen/UmfPackSupport>
#include <unsupported/Eigen/SparseExtra>


*/

#include "quadrule.hpp"

// Import Sparse Eigen numerical linear algebra library.
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>



///////////////////////////////////////////////////////////////////////////////
// USING FROM THE STL
///////////////////////////////////////////////////////////////////////////////
//use what we need from the standard namespace
using
std::cout;
using
std::cin;
using
std::endl;
using
std::ofstream;
using
std::ifstream;

///////////////////////////////////////////////////////////////////////////////
//  TYPEDEFS
//////////////////////////////////////////////////////////////////////////////


// Matrix and Vector type
typedef Eigen::SparseMatrix<double> ddpSparseMatrix_type;
typedef Eigen::DiagonalMatrix<double,Eigen::Dynamic> ddpDiagonalMatrix_type;
typedef Eigen::VectorXd ddpDenseVector_type;
typedef Eigen::SparseVector<double> ddpSparseVector_type;


//Linear solver uses UMFPACK
//typedef Eigen::UmfPackLU<Eigen::SparseMatrix< double > > ddpSolver_type;

// Linear solver uses EIGEN
typedef Eigen::SparseLU<Eigen::SparseMatrix<double> > ddpSolver_type;


// Number Of Derivatives
typedef enum {zero,one} ddpNumDeriv_type;

//TODO: Inconsistency between MX and Mixed, fix this
typedef enum {DG, Mixed, Points} ddpBijFlag_type;

// Direction type
typedef enum {Plus, Minus} ddpDirection_type;

// Charge Sign type
typedef enum {Positive, Negative} ddpChargeSign_type;

// Map to allow chargeSign type to be read in and mapped from string to enum
typedef std::map<std::string, ddpChargeSign_type> ddpChargeSignMap_type;

// Bijections for DOF to point values
typedef std::map < std::pair <int, int >, int > ddpBijForw_type;
typedef std::map < int, std::pair <int, int > > ddpBijBack_type;

// Status of whether coupling of drift diffusion to poisson is on or off
typedef enum {On, Off} ddpElecFieldCoupling_type;

// Map to allow couplingToPoisson to be on or off from string to enum for read in
typedef std::map<std::string, ddpElecFieldCoupling_type> ddpElecFieldCouplingMap;

// Status of whether the device is illuminated or not
typedef enum {Illuminated, Dark} ddpIlluminationStatus_type;

// Map to allow couplingToPoisson to be on or off from string to enum for read in
typedef std::map<std::string, ddpIlluminationStatus_type> ddpIlluminationMap;

// Boundary Condition Type
typedef enum {Dirichlet, Robin} ddpBCType_type;

// MATERIAL TYPE STILL NEED 9.23.2014
typedef enum {Semiconductor, Electrolyte} ddpMaterial_type;

typedef enum {Original, Continuation} ddpIVP_Type;

typedef std::map<std::string, ddpIVP_Type> ddpIVPMap;

//////////////////////////////////////////////////////////////////////////////
// STRUCTURES
/////////////////////////////////////////////////////////////////////////////


/** \brief Element type contains all the information about each element.*/
typedef struct ddpElement_type
{
    double Left;
    double Right;
    double Delta;
    int OrderMX;
    int OrderDG;
    int GaussLegendreNumPoints;
    double * GaussLegendrePoints;
    double * GaussLegendreWeights;

} ddpElement_type;

/** \brief Holds the end points of the domain.*/
typedef struct ddpDomain_type
{
    double LeftEndPoint;
    double RightEndPoint;

} ddpDomain_type;


/** \brief Holds the grid, with a list of all the elements and mesh info.*/
typedef struct ddpGrid_type
{
    double DeltaxMin;
    int OrderDGMax;
    ddpDomain_type Domain;
    int NumElementsNoGhost;
    int NumElementsWithGhost;
    int NumBoundaryElements;
    double BoundaryLayerWidth;

    /// array of elements.
    ddpElement_type * ElementList;

} ddpGrid_type;

/** Contains all the probelm info and will be passed around to
	* pick up coefficents. Will be read in from the input file and
	* stored here before assigned permanently. */
typedef struct ddpProblemInfo_type
{
    int MaxOrderDG;
    int MaxOrderMX;
    double TimeInitial;
    double TimeFinal;
    int GaussLegendreNumPoints;
    int NumFrames;
    double temperature;
    double electronCharge;
    double vacuumPermittivity;
    double semiCondRelativePerm;
    double electrolyteRelativePerm;
    double BoltzmannConst;
    double thermalVoltage;
    double characteristicLength;
    double characteristicTime;
    double characteristicDensity;
    double intrinsicDensity;
    ddpElecFieldCoupling_type ElecFieldCouplingStatus;
    double diffInEndPTS;
    double PhysicalStartEndPT;
    double appliedBias;
    double Absorption_Coeff;
    double Photon_Flux;
    double BuiltInBias;
    ddpIlluminationStatus_type IlluminationStatus;
    ddpIVP_Type IVP_Type;
    double increase_time_step_factor;

} ddpProblemInfo_type;


/** Defines the basis functions: elem number, order and DG vs Mixed. */
typedef struct ddpBasisFunction_type
{
    int element;
    int order;
    enum {psi, upsilon} family;

} ddpBasisFunction_type;


/** \brief Bijections from local DOF (elem, order) to global DOF index. */
typedef struct ddpBijection_type
{
    ddpBijForw_type DGForward;
    ddpBijBack_type DGBackwrd;
    ddpBijForw_type MXForward;
    ddpBijBack_type MXBackwrd;
    ddpBijForw_type PTForward;
    ddpBijBack_type PTBackwrd;
} ddpBijection_type;

/** \brief Vandemonde Matrices used to take vector of DOFs
 * 	 of a functions to point values of that function. */
typedef struct ddpVandeMondeMatrices_type
{
    // regular versions
    ddpSparseMatrix_type globalVandeMondeDG;
    ddpSparseMatrix_type globalVandeMondeDGPrime;
    ddpSparseMatrix_type globalVandeMondeMX;
    ddpSparseMatrix_type globalVandeMondeMXPrime;
    ddpSparseMatrix_type globalVandeMondeFluxMX;

    // transposed Versions
    ddpSparseMatrix_type globalVandeMondeDGTransposed;
    ddpSparseMatrix_type globalVandeMondeDGPrimeTransposed;

} ddpVandeMondeMatrices_type;

/** \brief Flux Matrices used build up fluxes in LDG method,
 * see ddpMakeDiffusiveFluxProperties. */
typedef struct ddpDGFluxMatrices_type
{
    // regular versions
    ddpSparseMatrix_type plusplus;
    ddpSparseMatrix_type plusminus;
    ddpSparseMatrix_type minusplus;
    ddpSparseMatrix_type minusminus;
    // transposed versions
    ddpSparseMatrix_type plusplusTransposed;
    ddpSparseMatrix_type plusminusTransposed;
    ddpSparseMatrix_type minusplusTransposed;
    ddpSparseMatrix_type minusminusTransposed;

} ddpDGFluxMatrices_type;


/** \brief Structure that contains the information and function
 * of each boundary condition. */
typedef struct ddpBoundaryCondition_type
{
    // Dirichlet, Neumann, or Robin BC Type
    ddpBCType_type BCType;

    // Ohmic BC Value
    double OhmicValue;

    // Allow For Robin BC
    double RobinValue;

    // Testing BC Function Pointers
    double (* TestDirichletValue)(const double & t);
    double (* TestRobinValue)(const double & t);

    // This to get over the x = -1 or x = 0 problem
    double xLeftEndPoint;

} ddpBoundaryCondition_type;

/** \brief  structure that contains the information and function of
 * the generation function. */
typedef struct ddpGenerationFunction_type
{
    double AbsorpCoeff;
    double PhotonFlux;
    double ScaledAbsorpCoeff;

    double GenFun(const double & x)
    {
        return AbsorpCoeff * PhotonFlux * exp(ScaledAbsorpCoeff * (x) );
    }

} ddpGenerationFunction_type;


/** \brief All the properties of the carriers that do not change
 *  with time.*/
/** This includes the matrices, vectors and paramater values.
* The matrix will be of the form,
*   \f[ \left[ \begin{matrix}
*			 \mu^{-1} A & B_{1} \\
*			 B_{2} & \frac{1}{\Delta t} M + C
*			 \end{matrix} \right]  \f]
*
* where,
*
*
*	 \f[ A(\textbf{p},\textbf{q} ) \; = \;
*				\int_{\Omega} \ \textbf{p} \ \cdot \textbf{q} \ dx \f]
*	 \f[ M(v,u ) \; = \;
*				\int_{\Omega} \ v \ u \ dx \f]
*
*	 \f[ B_{1}(v,\textbf{q} ) \; = \;
*				\int_{\Omega} \ \nabla \ v  \ \textbf{q} \ dx
*				\ + \
*				\int_{\partial \Omega_{D}} v  \  \textbf{q}
*				\ \cdot \boldsymbol \eta \ ds	\f]
*
*	 \f[ B_{2}(\textbf{p},u ) \; = \;
*				\int_{\Omega} \ \nabla \ \cdot \ \textbf{p} \ u \ dx
*				\ + \
*				\int_{\partial \Omega_{N}} \textbf{p}
*							\cdot \boldsymbol \eta \ u \ ds	\f]
*
*/
typedef struct ddpCarrierProperties_type
{
    // Mass, Stiffness and flux Matrices
    ddpSparseMatrix_type MassU;
    ddpDiagonalMatrix_type InvMassU;

    ddpSparseMatrix_type MassQ;
    ddpDiagonalMatrix_type InvMassQ;

    // LDG Matrices
    ddpSparseMatrix_type FluxRightUFromQ;
    ddpSparseMatrix_type FluxLeftUFromQ;
    ddpSparseMatrix_type FluxRightQFromU;
    ddpSparseMatrix_type FluxLeftQFromU;

    ddpSparseMatrix_type TotalFluxFromBCRHS;
    ddpSparseMatrix_type TotalFluxFromBCPenalty;

    // forward Euler matrices for U
    ddpSparseMatrix_type StiffUFromQ;
    ddpSparseMatrix_type TotalUFromQRHS;

    // forward Euler matrices for Q
    ddpSparseMatrix_type StiffQFromU;
    ddpSparseMatrix_type TotalQFromURHS;

    // steady state calculations matrices
    ddpSparseMatrix_type MassUToSystem;
    ddpSparseMatrix_type B1;
    ddpSparseMatrix_type B2;
    ddpSparseMatrix_type C;
    ddpSparseMatrix_type Dir_RHS;
    ddpSparseMatrix_type Robin_RHS;
    ddpSparseMatrix_type Penalty_RHS;
    ddpSparseMatrix_type LeftFlux_Inside;
    ddpSparseMatrix_type RightFlux_Inside;
    ddpSparseMatrix_type RecomboToSystem;
    ddpSparseMatrix_type DriftToSystem;

    // For doing fully implict robin boundary
    ddpSparseMatrix_type LeftBoundary_LookInside;
    ddpSparseMatrix_type RightBoundary_LookInside;

    // Boundary condition structs
    ddpBoundaryCondition_type BCLeft;
    ddpBoundaryCondition_type BCRight;

    // Carrier's mobility and diffusivity
    double Mobility;
    double Diffusivity;
    double MaxMobility;
    double MaxDiffusivity;


    // TODO: DO WE NEED THESE?
    double Scale4Diffusivity;
    double Scale4Mobility;

    double CurrentScale;
    double RecomboScale;
    double LengthScale;
    double TimeScale;

    // Other properties
    ddpChargeSign_type ChargeSign;
    double Sign4Force;
    double Sign4Poisson;
    double RecombinationTime;
    double TransferRate;
    double EquilibriumDensity;
    double AbsorptionCoeff;
    double PhotonFlux;
    double IntrisincDensity;
    double Scale4Recombo;
    double Beta; // for LDG fluxes
    double Tau; // for LDG fluxes
    double RecombinationVelocity;
    ddpMaterial_type Material;
} ddpCarrierProperties_type;


/** \brief  Carrier properties that do change with time. */
/** These properties include DOFS and matrices/vectors to set up
	*  du/dt = F(u,x,t) that get updated at every time step. */
typedef struct ddpCarrierState_type
{

    // Dofs
    ddpDenseVector_type uDof;
    ddpDenseVector_type qDof;
    ddpDenseVector_type qDof_prev;

    // Work vectors to be used in makeDot to make uDotDof
    ddpDenseVector_type BC_Dir_Input;
    ddpDenseVector_type BC_Robin_Input;
    ddpDenseVector_type QRHS;
    ddpDenseVector_type RHSFromQ;
    ddpDenseVector_type RHSFromElec;
    ddpDenseVector_type RHSFromGeneration;
    ddpDenseVector_type RHSFromRecombination;
    ddpDenseVector_type Generation_RHS;
    ddpDenseVector_type RHSTotal;

    // Steady state matrices
    ddpSparseMatrix_type PE;
    ddpSparseMatrix_type LDG_ABig;
    ddpDenseVector_type LDG_RHS;
    ddpDenseVector_type BigSolution;
    ddpDenseVector_type OldBigSolution;

    ddpDenseVector_type ElecPTS;
    ddpDenseVector_type ElecPTSTimesWeights;
    ddpSparseVector_type SparseElecPTSTimesWeights;
    ddpDenseVector_type DriftTerm;

    // Work Matrices for drift componenet
    ddpDiagonalMatrix_type ElecTimesWeightsMatrix;


    ddpSolver_type Solver_LDG;

} ddpCarrierState_type;


/** \brief This structure helps us store the read in values
	* and then attach them to the appropriate structures later on. */
typedef struct ddpCarrierConstants_type
{
    // Electron constants
    double electron_Mobility;
    ddpChargeSign_type electron_ChargeSign;
    double electron_RecombinationTime;
    double electron_TransferRate;
    double electron_RecombinationVelocity;

    // hole constants
    double hole_Mobility;
    ddpChargeSign_type hole_ChargeSign;
    double hole_RecombinationTime;
    double hole_TransferRate;
    double hole_RecombinationVelocity;

    // reductant constants
    double reductant_Mobility;
    ddpChargeSign_type reductant_ChargeSign;

    // oxidant constants
    double oxidant_Mobility;
    ddpChargeSign_type oxidant_ChargeSign;

} ddpCarrierConstants_type;


/** \brief Properties of the potential, electric field structure that do not change in time.*/
/** 	The matrix will be of the form,
			*   \f[ \left[ \begin{matrix}
			*			 A & B \\
			*			 B^{T} & 0
			*			 \end{matrix} \right]  \f]
			*
			* where,
			*
			*	 \f[ A(\textbf{p},\textbf{D}) \; = \;
			*				\int_{\Omega} \  \textbf{p} \  \cdot \  \textbf{D} \ dx \f]
			*
			*	 \f[ B(\textbf{p},\Phi) \; = \;
			*				\int_{\Omega} \ \nabla \ \cdot \ \textbf{p} \  \Phi  \ dx \f]
			*
			*	 \f[ B^{T}(v,\textbf{D} ) \; = \;
			*				\int_{\Omega} \ \nabla v \ \cdot \ \textbf{D} \ dx \f]
			*
			*
			*/
typedef struct ddpPoissonProperties_type
{

    // Matrices
    ddpSparseMatrix_type ABig;
    ddpSparseMatrix_type A00;
    ddpSparseMatrix_type A01;
    ddpSparseMatrix_type A10;
    ddpSparseMatrix_type C;
    ddpSparseMatrix_type VBig;
    ddpSparseMatrix_type VRHS;
    ddpSparseMatrix_type CRHS;

    // Constants
    double SemiconductorPerm;
    double ElectrolytePerm;
    double Lambda;
    double coupledOrNot;
    ddpDiagonalMatrix_type EpsilonInverse;

} ddpPoissonProp_type;


/** \brief Propeties, constants, matrices vectors that
	* are for potential and electric field that do not change	with time.*/
typedef struct ddpPoissonState_type
{
// Data Members

    // Dofs
    ddpDenseVector_type elecDof;
    ddpDenseVector_type potDof;

    // work vectors to be fed into create the actual dofs
    ddpDenseVector_type carrierDof;
    ddpDenseVector_type BCInput;

    ddpDenseVector_type RHSFromBC;
    ddpDenseVector_type RHSFromCandU;

    ddpDenseVector_type RHSTop;
    ddpDenseVector_type RHSBottom;
    ddpDenseVector_type RHSTotal;
    ddpDenseVector_type Soln;

    ddpSolver_type solverABig;

} ddpPoissonState_type;

#endif

