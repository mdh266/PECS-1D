#ifndef _POISSON_H_
#define _POISSON_H_

#include "includes.hpp"
#include "main.hpp"
#include "dopingprofile.hpp"
#include "carrier.hpp"


typedef struct ddpPoisson_type
{
    /** Assembles each cells local left hand side matrix for a constant
    *		Debeye length. Also has can compute the error of the approximation
    * 	for testing. The data is then stored in ddpPoissonProperties_type PoissonProps.
    *
    *   The general problem is,
    * 	\f[ \begin{align}
    *		\epsilon_{r}^{-1} \ \textbf{D} \ + \ \nabla \Phi \ &= \ 0  && \text{in} \; \Omega \\
    *	  \ \nabla \ \cdot \ \textbf{D} \  &= \frac{1}{\lambda^{2}}
    *		f(\textbf{x})
    *						&& \text{in} \; \Omega \\
    *		\textbf{D} \ \cdot \ \boldsymbol \eta \ &= \ 0
    *						&& \text{on} \; \partial \Omega_{N} \\
    *		\Phi \ &=  \ \Phi_{D}
    *							 && \text{on} \; \partial \Omega_{D}
    *   \end{align} \f]
    *
    * For \f$\lambda^{2} \ = \ \frac{\Phi^{*} \epsilon }{L^{2} q C^{*}}\f$.
    *
    *	This becomes the problem in the weak formulation:
    *
    * Find
    *	\f$( \ \Phi \ , \ \textbf{D} \ ) \ \in
    *	\left( \ \text{W} \ \times \ [ \ 0 , \ T \ ] ,
    *	\ \textbf{V}^{d} \ \times \
    *  	[ \ 0 , \ T \ ] \ \right) \f$ such that:
    *
    * 	\f[ \begin{align}
    * \ \left(  \textbf{p} \ , \ \epsilon_{r}^{-1} \ \textbf{D} \right)_{\Omega}
    *	\  - \
    * \left( \ \boldsymbol \nabla \cdot \ \textbf{p}  \ , \  \Phi  \right)_{\Omega}
    *  \; &= \;
    *  - \langle \textbf{p} \ , \ \Phi_{D} \rangle_{\Gamma_{D}}   \\
    * -  \left(  v \ , \ \boldsymbol \nabla  \cdot \textbf{D} \right)_{\Omega} \;
    *  &= \; -
    *	\frac{1}{\lambda^{2}} \
    * \left(  v,\ f(\textbf{x})  \right)_{\Omega}
    *	\end{align} \f]
    *
    *
    *  For all \f$( v \  , \ \textbf{p}  ) \, \in \, W \, \times\, \textbf{V}^{d}\f$.
    *
    *
    *	We obtain the electric field \f$-\boldsymbol \nabla \Phi\f$ by setting,
    *
    *	\f[
    * - \boldsymbol \nabla \Phi \ = \ \epsilon_{r}^{-1}  \ \textbf{D}
    *	\f]
    *
    *
    *  This method only assembles the left hand side of the weak formulation.
    */

    //////////////////////////////////////////////////////////////////////////////
    //DATA MEMBERS
    //////////////////////////////////////////////////////////////////////////////
    ddpPoissonProperties_type PoissonProps;
    ddpPoissonState_type PoissonState;

    ddpDomain_type Poisson_domain;
    ddpGrid_type Poisson_grid;

    // Weights and Points Vectors
    ddpSparseVector_type
    weightsSparse,
    PTSSparse;

    ddpDenseVector_type
    PTSDense,
    weightsDense;

    // Bijections
    ddpBijection_type Bijections;

    // VandeMondeMatrices, dont need if memory beocomes issue
    ddpVandeMondeMatrices_type VandeMondeMatrices;
    ddpDGFluxMatrices_type DGFluxMatrices;

    /// Matrices to turn carrier DOF vectors into 2 x size
    ddpSparseMatrix_type SemicondcutorToTotal;
    ddpSparseMatrix_type ElectrolyteToTotal;

    // Matrices to turn electric field DOF to 1/2 size
    ddpSparseMatrix_type EFTotalToSemiconductor;
    ddpSparseMatrix_type EFTotalToElectrolyte;

    // Matrices to turn electric field DOF to 1/2 size
    ddpSparseMatrix_type PotTotalToSemiconductor;
    ddpSparseMatrix_type PotTotalToElectrolyte;


    //////////////////////////////////////////////////////////////////////////////
    // Data Functions
    //////////////////////////////////////////////////////////////////////////////

    // Boundary Condition Pointers
    double (* testBCDirLeft) (double const & x);
    double (* testBCDirRight) (double const & x);
    double leftVoltage;
    double rightVoltage;


    //////////////////////////////////////////////////////////////////////////////
    // Structure Methods
    //////////////////////////////////////////////////////////////////////////////

    int
    initialize(const ddpGrid_type & grid, const ddpProblemInfo_type & problem,
               const ddpCarrierConstants_type & carrierConstants);

    int
    initializeWithInterface(const ddpGrid_type & left_grid,
                            const ddpGrid_type & right_grid,
                            const ddpProblemInfo_type & problem,
                            const ddpCarrierConstants_type & carrierConstants);

    int
    remakeProperties(const ddpProblemInfo_type & problem,
                     const ddpGrid_type& grid);

    int
    setDirichletBoundaryConditions(double (* testBCLeft) (const double & x),
                                   double (* testBCRight) (const double & x) );

    int
    setBias(const ddpProblemInfo_type & problem);

    int
    getElecFieldVals(ddpDenseVector_type & EPTS) const;

    int
    getElectrolyteElecFieldVals(ddpDenseVector_type & EPTS) const;

    int
    getSemiconductorElecFieldDOFS(ddpDenseVector_type & EPTS) const;

    int
    getElectrolyteElecFieldDOFS(ddpDenseVector_type & EPTS) const;

    int
    getSemiconductorPotentialDOFS(ddpDenseVector_type & EPTS) const;

    int
    getElectrolytePotentialDOFS(ddpDenseVector_type & EPTS) const;

    ///////////////////////////////////////////////////////////////////
    // BIPOLAR MODEL
    ///////////////////////////////////////////////////////////////////


    int
    solveSystem(ddpCarrier_type const & carrier1,
                ddpCarrier_type const & carrier2,
                ddpDopingProfile_type const & dopingProfile,
                const double & time,
                ddpProblemInfo_type const & problem );

    ///////////////////////////////////////////////////////////////////
    // UNIPOLAR MODEL
    ///////////////////////////////////////////////////////////////////

    int
    solveSystem(ddpCarrier_type const & carrier1,
                ddpDopingProfile_type const & dopingProfile,
                const double & time,
                ddpProblemInfo_type const & problem );

    ///////////////////////////////////////////////////////////////////
    // BOTH MODELS
    ///////////////////////////////////////////////////////////////////

    int
    solveSystem(ddpDenseVector_type const & carrier1Dof,
                ddpDopingProfile_type const & dopingProfile,
                const double & time,
                ddpProblemInfo_type const & problem );


    ///////////////////////////////////////////////////////////////////
    // SEMICONDUCTOR-ELECTROLYTE MODEL
    ///////////////////////////////////////////////////////////////////


    int
    solveSystem(ddpCarrier_type const & electrons,
                ddpCarrier_type const & holes,
                ddpCarrier_type const & reductants,
                ddpCarrier_type const & oxidants,
                ddpDopingProfile_type const & dopingProfile,
                const double & time,
                ddpProblemInfo_type const & problem);



}  ddpPoisson_type;

#endif

