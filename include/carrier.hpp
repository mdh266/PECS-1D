#ifndef _CARRIER_H_
#define _CARRIER_H_

#include "includes.hpp"
#include "main.hpp"

typedef struct ddpCarrier_type
{
    /** This class will build portions of the LDG local cell matrix for the general
    	* (non-dimensional) drift diffusion equation:
    	*
    	*	\f[ \begin{align}
    	* 	u_{t} \  - \  \boldsymbol \nabla \ \cdot
    	* 	\ \mu \left( s \boldsymbol \nabla \Phi u \ + \ \boldsymbol \nabla u \ \right)
    	* 	\; &= \;
    	* 	R(u) + G	&& \text{in} \;  \Omega   \\
    	* 	u \; &= \; u_{D} &&  \text{on} \;  \Omega_{D}     \\
    	* 	- \mu \left(s \boldsymbol \nabla \Phi \ u \ + \ \boldsymbol \nabla u \ \right)
    	* 	\ \cdot \ \boldsymbol \eta
    	*	\;  &= \; K (u) && \text{on} \; \partial \Omega_{N}
    	* 	\end{align} \f]
    	*
    	* 	We rewrite this in mixed form:
    	*
    	* \f[ \begin{align}
    	*			u_{t} \ + \ \nabla \ \textbf{q} \ &= \ R(u) \ + G && \text{in} \ \Omega \\
    	*			\mu^{-1} \ \textbf{q} \ & =
    	*								 \ -s \nabla \Phi \ u \ - \nabla u && \text{in} \ \Omega \\
    	*			\mu^{-1} \ \textbf{q} \ \cdot \boldsymbol \eta &=
    	*									\ K(u) && \text{on} \ \partial \ \Omega_{N} \\
    	*			u \ &= \ u_{D} && \text{on} \ \partial \Omega_{D}
    	*	\end{align} \f]
    	*
    	* 	The weak formulation for IMEX will be:
    	*
    	*	Find \f$(u, \textbf{q}) \in W \times [t^{k-1}, t^{k}]
    	*			\times \textbf{W}^{d} \times[ t^{k-1}, t^{k}] \f$ such that,
    	*
    	*  \f[ \begin{align}
    	*	\frac{1}{\Delta t} \left(  v , u^{k}  \right)
    	* 	-
    	*	\tau \langle  [[ \ v \ ]] ,
    	*				 [[ u^{k} ]] \rangle_{\mathcal{E}_{h}^{0}}
    	*	-
    	*	\left( \boldsymbol \nabla v  ,  \textbf{q}^{k}  \right)
    	*	+
    	*	 \langle [[ \ v \ ]] ,
    	*	\{ \textbf{q}^{k} \} \rangle_{\mathcal{E}_{h}^{0} \cap \partial \Omega_{D}}
    	*	\ &= \
    	*	\left( v ,  R(u^{k-1})  + G \right) -
    	*	\langle   v, K( u^{k-1})    \rangle_{\Sigma}   \nonumber \\
    	*	 - \left(  \boldsymbol \nabla \cdot \textbf{p} ,   u^{k-1} \right)
    	*  \ - \
    	*	\langle  [[ \,  \textbf{p} \, ]] ,
    	*	\{  u^{k} \}  \rangle_{\mathcal{E}^{0}_{h} \cap \partial \Omega_{N}}
    	*	\ + \
    	*	 \left( \textbf{p} , \textbf{q}^{k} \right)
    	*	\ &= \
    	*	 +
    	*	\left( s \textbf{P}  \cdot \boldsymbol \nabla \Phi , u^{k-1} \right)
    	*	-
    	*	\langle  \textbf{p}   ,  u_{D}  \rangle_{ \partial \Omega_{N} }
    	*	\end{align} \f]
    	*
    	*  For all \f$(v,\textbf{p}) \in W \times \textbf{W}^{d})\f$.
    	*
    	* 	The corresponding matrix will be of the form,
    	*   \f[ \left[ \begin{matrix}
    	*		\mu^{-1} A & B_{1} + F_{1} \\
    	*		 B_{2} + F_{2} & \frac{1}{\Delta t} M + C
    	*	 \end{matrix} \right]  \f]
    	* 	This matrix will assembled once in initialize, setSover.  It will be stored in ddpCarrierProperties carrierProps.
    	*  	The corresponding right hand side vector will assembled at every time step
    	* 	in AssembleRHS and stored in ddpCarrierState_type carrierState.
    	*
    	* 	\note We use IMEX time stepping so all non-linear terms and drift terms are
    	*				time lagged and therefore on the right hand side. While they are
    	*				are built in parallel, this place takes place outside of this class and
    	*				in ddpTimeStepping_type.
    	*	\note Only the recombination term will be assembled in parallel using OpenMP.
    	*
    	*	\note The LDG flux matrices \f$F_{1}\f$ and \f$F_{2}\f$.
    	*				are built sequentially and this occurs in a loop
    	*  			outside this class, but calls the functions ddpMakeCarrierProperties and ddpMakeDiffusiveFluxProperties.
    	*

    	*/
    //////////////////////////////////////////////////////////////////////////////
    // Data Members
    ///////////////////////////////////////////////////////////////////////////////

    /** \brief Carrier vectors and properties that do change in time. */
    ddpCarrierState_type carrierState;

    /** \brief Carrier Matrices and coefficients dont change in time.*/
    ddpCarrierProperties_type carrierProps;

    //TODO: These and others could be move out if they become an issue.
    /** \brief Quadrature weights and points in sparse vectors.*/
    ddpSparseVector_type
    weightsSparse,
    PTSSparse;

    /** \brief Quadrature weights and points in dense vectors.*/
    ddpDenseVector_type
    PTSDense,
    weightsDense;

    /** \brief Bijections from local DOFs to global DOFs.*/
    ddpBijection_type Bijections;

    // VandeMondeMatrices, dont need if memory beocomes issue
    ddpVandeMondeMatrices_type VandeMondeMatrices;
    ddpDGFluxMatrices_type DGFluxMatrices;

    // Workspace vectors
    ddpDenseVector_type uNext; //used in forwardEuler to do time stepping

///////////////////////////////////////////////////////////////////////////////
// Data Function Pointers
///////////////////////////////////////////////////////////////////////////////

    // Initial Condition Function Pointers
    double (* InitialConditionsU)(double const & x);

///////////////////////////////////////////////////////////////////////////////
// Generation And Recombination Objects
//////////////////////////////////////////////////////////////////////////////
    ddpGenerationFunction_type GenerationFunction;

///////////////////////////////////////////////////////////////////////////////
// Class Methods
///////////////////////////////////////////////////////////////////////////////

    // initialization
    int initialize(const ddpGrid_type & grid,
                   const ddpProblemInfo_type & problem,
                   const ddpCarrierConstants_type & carrierConstants,
                   char const * nameOfCarrier);


    int setSolver(ddpGrid_type const & grid,
                  ddpProblemInfo_type const & problem);


    int setTestGenerationRHS(ddpGrid_type const & grid,
                             double (*GenerationFunction)(const double & x) );

    int setGenerationRHS(ddpGrid_type const & grid,
                         const ddpIlluminationStatus_type & IlluminationStatus);

    int assembleLDGMatrices(ddpGrid_type const & grid,
                            ddpProblemInfo_type const & problem);

    int setInitialConditions(ddpGrid_type & grid,
                             double (* InitialFunction) (const double & t) );


    // setting boundary conditions
    int setLeft_DirBC(double (* DopingProfile)(const double & x) );

    int setRight_DirBC(double (* DopingProfile)(const double & x) );

    int setLeft_RobincBC();

    int setRight_RobinBC();

    int	setLeft_RobinValue(const double & r);

    int	setRight_RobinValue(const double & r);

    int setLeftTestDirichletBC(double (* BCLeft_Input) (const double & t) );

    int setRightTestDirichletBC(double (* BCRight_Input) (const double & t) );

    int setLeftTestRobincBC(double (* BCLeft_Input) (const double & t) );

    int setRightTestRobinBC(double (* BCRight_Input) (const double & t) );

    // others
    int setEquilibriumDensity(double (* DopingProfile)(const double & x) );

    double getEquilibriumDensity();

    double getTransferRate();

    double getInterfaceDensity();


    ///////////////////////////
    // Steady state solvers	//
    //////////////////////////

    // Monolithic solver
    int	Solve(ddpDenseVector_type const & ElecFieldDof,
              ddpProblemInfo_type const & problem,
              ddpGrid_type const & grid,
              double const & timeCurrent);

    // Uses Shur Complement
    int	ShurSolve(ddpDenseVector_type const & ElecFieldDof,
                  ddpProblemInfo_type const & problem,
                  ddpGrid_type const & grid,
                  double const & timeCurrent);

    ///////////////////////////
    // Time Steppers			 	//
    //////////////////////////

    int AssembleLDGSystem(ddpProblemInfo_type const & problem,
                          ddpGrid_type const & grid,
                          double const & deltaT);

    int	AssembleLDGSystem(ddpProblemInfo_type const & problem,
                          ddpGrid_type const & grid,
                          double const & interface_value,
                          double const & deltaT);

    // NOTE: should be used before assembleLDGSystem
    int UpdateImplicitDriftTerm(ddpDenseVector_type & ElecFieldDof,
                                const ddpProblemInfo_type & problem,
                                const ddpGrid_type & grid);

    // NOTE: updateDriftTerm should be used before assembleRHS
    int UpdateExplicitDriftTerm(ddpDenseVector_type & ElecFieldDof);

    int	AssembleRHS(double const & timeCurrent,
                    double const & deltaT);

    int	BackwardEuler();

    bool check_converged(const double & tol);

} ddpCarrier_type;

#endif
