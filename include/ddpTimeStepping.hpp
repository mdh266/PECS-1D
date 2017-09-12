#ifndef _ddpTimeStepping_H_
#define _ddpTimeStepping_H_

#include "carrier.hpp"
#include "poisson.hpp"
#include "includes.hpp"
#include "ddpComputeDeltaT.hpp"
#include "ddpRecombination.hpp"

/** \brief Time stepping routines for the ENTIRE SYSTEMS. */
typedef class ddpTimeStepping_type
{
	protected:
		double electrons_e, holes_e; 	// equilibrium densities
		double k_et, k_ht; 				// carrier transfer rates
		double k_f, k_b;				// forward and backward transfer rates

		// Carrier densites at the interface
		double electrons_interface,
					 holes_interface,
					 reductants_interface,
					 oxidants_interface,
					 electron_steady_state,
					 hole_steady_state,
					 left_value;
		
		// Workspace vectors
		ddpDenseVector_type SemiconductorElecFieldDof;
		ddpDenseVector_type ElectrolyteElecFieldDof;
		
		ddpRecombination_type SRH_Recombination;

		ddpDenseVector_type electron_q_values;
		ddpDenseVector_type hole_q_values;
		ddpDenseVector_type q_values;
 
	public:
		/** \brief Constructor for all the time stepping methods.*/
		ddpTimeStepping_type(ddpCarrier_type & electrons,
												 ddpCarrier_type & holes,
												 ddpCarrier_type & reductants,
												 ddpCarrier_type & oxidants,
												 ddpPoisson_type & Poisson,
												 const ddpDopingProfile_type & DopingProfile,
												 const ddpProblemInfo_type & problem,
												 const double & tCurrent);
		
		/** \brief Default destructor. */					 
		~ddpTimeStepping_type();
		

		/** \brief Assembles all the LDG matrices which dont change in time. */
		int AssembleLDGMatrices(ddpCarrier_type & electrons,
														 ddpCarrier_type & holes,
														 ddpCarrier_type & reductants,
														 ddpCarrier_type & oxidants,
														 const ddpProblemInfo_type & problem,
														 const ddpGrid_type & semiconductor_grid,
														 const ddpGrid_type & electrolyte_grid,
						 								 const double & tDelta);

		int AssembleLDGMatrices(ddpCarrier_type & electrons,
									 ddpCarrier_type & holes,
									 ddpCarrier_type & reductants,
									 ddpCarrier_type & oxidants,
									 const ddpProblemInfo_type & problem,
									 const ddpGrid_type & semiconductor_grid,
									 const ddpGrid_type & electrolyte_grid,
	 								 const double & tDelta_semi,
									 const double & tDelta_elec);

		/** \brief Alternating-IMEX method.*/
		/** The drift component is always solved using explicit density values.  The density 
 * values on the itnerface are treated in an alternating way. Using implicit values for the
 * carriers own density values on the interface and all the other densities values using explicit values.
 * In this way we cyle through the algorithm  \n
 * 0.) Update electric field and potential using all densities previous time steps values. \n
 * 1.) Update electrons and holes using the other carries as well as redox previous time steps densities. \n
 * 2.) Update the reductants using the new values of electrons and holes and the previous time steps density
 * values for the oxidants. \n
 * 3.) Update the oxidants using the new values of electrons, holes and reductants.*/
		int IMIMEX_Solve(ddpCarrier_type & electrons,
											ddpCarrier_type & holes,
											ddpCarrier_type & reductants,
											ddpCarrier_type & oxidants,
											ddpPoisson_type & Poisson,
											const ddpGrid_type & semiconductor_grid,
											const ddpGrid_type & electrolyte_grid,
											const ddpProblemInfo_type & problem,
											const double & tCurrent,
											const double & tDelta);

		int IMEXEX_Solve(ddpCarrier_type & electrons,
											ddpCarrier_type & holes,
											ddpCarrier_type & reductants,
											ddpCarrier_type & oxidants,
											ddpPoisson_type & Poisson,
											const ddpGrid_type & semiconductor_grid,
											const ddpGrid_type & electrolyte_grid,
											const ddpProblemInfo_type & problem,
											const double & tCurrent,
											const double & tDelta);


		int	AssembleBackwardEulerRHS(ddpCarrier_type & electrons,
														 ddpCarrier_type & holes,
														 ddpCarrier_type & reductants,
														 ddpCarrier_type & oxidants,
														 const ddpGrid_type & semiconductor_grid,
														 const double & tCurrent,
						 								 const double & tDelta);

		int	AssembleBackwardEulerSemiconductorRHS(ddpCarrier_type & electrons,
														 ddpCarrier_type & holes,
														 ddpCarrier_type & reductants,
														 ddpCarrier_type & oxidants,
														 const ddpGrid_type & semiconductor_grid,
														 const double & tCurrent,
						 								 const double & tDelta);

		int	AssembleBackwardEulerElectrolyteRHS(ddpCarrier_type & electrons,
														 ddpCarrier_type & holes,
														 ddpCarrier_type & reductants,
														 ddpCarrier_type & oxidants,
														 const double & tCurrent,
						 								 const double & tDelta);

	 // NOTE: Needs to be called before AssembleLDGMatrices
		int UpdateImplicitDriftTerms(ddpCarrier_type & electrons,
										ddpCarrier_type & holes,
										ddpCarrier_type & reductants,
										ddpCarrier_type & oxidants,
										ddpPoisson_type & Poisson,
										const ddpProblemInfo_type & problem,
										const ddpGrid_type & grid);
	

	 // NOTE: Needs to be called before AssembleBackwardEulerRHS
		int UpdateExplicitDriftTerms(ddpCarrier_type & electrons,
										ddpCarrier_type & holes,
										ddpCarrier_type & reductants,
										ddpCarrier_type & oxidants,
										ddpPoisson_type & Poisson);

		int	UpdateExplicitSemiconductorDriftTerms(ddpCarrier_type & electrons,
									ddpCarrier_type & holes,
									ddpPoisson_type & Poisson);
	
	 // NOTE: Needs to be called before AssembleBackwardEulerElectrolyteRHS
		int UpdateExplicitElectrolyteDriftTerms(ddpCarrier_type & reductants,
									ddpCarrier_type & oxidants,
									ddpPoisson_type & Poisson);
	

		int BackwardEulerSolve(ddpCarrier_type & carrier1,
													 ddpCarrier_type & carrier2);



		int BackwardEulerSolve(ddpCarrier_type & electrons,
												ddpCarrier_type & holes,
												ddpCarrier_type & reductants,
												ddpCarrier_type & oxidants);

		bool check_converged(ddpCarrier_type & electrons,
												 ddpCarrier_type & holes,
												 const double & tol);

} ddpTimeStepper_type;
#endif