#include "../include/ddpTimeStepping.hpp"


ddpTimeStepping_type::
ddpTimeStepping_type(ddpCarrier_type & electrons,
                     ddpCarrier_type & holes,
                     ddpCarrier_type & reductants,
                     ddpCarrier_type & oxidants,
                     ddpPoisson_type & Poisson,
                     const ddpDopingProfile_type & DopingProfile,
                     const ddpProblemInfo_type & problem,
                     const double & tCurrent)
    : SRH_Recombination(electrons, holes)
{
    // set the equilibrium values
    electrons_e = electrons.carrierProps.EquilibriumDensity;
    holes_e	=	holes.carrierProps.EquilibriumDensity;

    // set the transfer rates
    k_et = electrons.getTransferRate();
    k_ht = holes.getTransferRate();

    int MXSubdomainDOF = Poisson.EFTotalToSemiconductor.rows();

    // set the sizes of the electric field dofs
    SemiconductorElecFieldDof = ddpDenseVector_type::Zero(MXSubdomainDOF);
    ElectrolyteElecFieldDof = ddpDenseVector_type::Zero(MXSubdomainDOF);
}

ddpTimeStepping_type::
~ddpTimeStepping_type()
{
    // Default
}


int
ddpTimeStepping_type::
AssembleLDGMatrices(ddpCarrier_type & electrons,
                    ddpCarrier_type & holes,
                    ddpCarrier_type & reductants,
                    ddpCarrier_type & oxidants,
                    const ddpProblemInfo_type & problem,
                    const ddpGrid_type & semiconductor_grid,
                    const ddpGrid_type & electrolyte_grid,
                    const double & tDelta)
{
    electrons.AssembleLDGSystem(problem, semiconductor_grid, tDelta);
    holes.AssembleLDGSystem(problem, semiconductor_grid, tDelta);
    reductants.AssembleLDGSystem(problem, electrolyte_grid, tDelta);
    oxidants.AssembleLDGSystem(problem, electrolyte_grid, tDelta);
    return 0;
}

int
ddpTimeStepping_type::
AssembleLDGMatrices(ddpCarrier_type & electrons,
                    ddpCarrier_type & holes,
                    ddpCarrier_type & reductants,
                    ddpCarrier_type & oxidants,
                    const ddpProblemInfo_type & problem,
                    const ddpGrid_type & semiconductor_grid,
                    const ddpGrid_type & electrolyte_grid,
                    const double & tDelta_semi,
                    const double & tDelta_elec)
{
    electrons.AssembleLDGSystem(problem, semiconductor_grid, tDelta_semi);
    holes.AssembleLDGSystem(problem, semiconductor_grid, tDelta_semi);

    reductants.AssembleLDGSystem(problem, electrolyte_grid, tDelta_elec);
    oxidants.AssembleLDGSystem(problem, electrolyte_grid, tDelta_elec);
    return 0;
}

int
ddpTimeStepping_type::
AssembleBackwardEulerRHS(ddpCarrier_type & electrons,
                         ddpCarrier_type & holes,
                         ddpCarrier_type & reductants,
                         ddpCarrier_type & oxidants,
                         const ddpGrid_type & semiconductor_grid,
                         const double & tCurrent,
                         const double & tDelta)
{
    ////////////////////////////////////////////////////////////////////
    // Compute the interface conditions
    ///////////////////////////////////////////////////////////////////

    // get interface density values
    electrons_interface = electrons.getInterfaceDensity();
    holes_interface = holes.getInterfaceDensity();
    reductants_interface = reductants.getInterfaceDensity();
    oxidants_interface = oxidants.getInterfaceDensity();


    // get the forward and backward transfer rates at this time step
    k_f = k_et *	(electrons_interface - electrons_e);
    k_b = k_ht *  (holes_interface - holes_e);

    // set the Robin boundary conditions at this time step
    electrons.setRight_RobinValue(k_f * oxidants_interface );
    holes.setRight_RobinValue( k_b * reductants_interface );

    reductants.setLeft_RobinValue(
        (k_f * oxidants_interface - k_b * reductants_interface ));

    oxidants.setLeft_RobinValue(
        (k_b * reductants_interface - k_f *oxidants_interface) );

    ////////////////////////////////////////////////////////////////////
    // Update Recombination
    ///////////////////////////////////////////////////////////////////
//    grvy_timer_begin("Assemble Recombination Term");
    SRH_Recombination.updateRecombination(semiconductor_grid, electrons, holes);
//    grvy_timer_end("Assemble Recombination Term");

    ////////////////////////////////////////////////////////////////////
    // Assemble RHSes
    ///////////////////////////////////////////////////////////////////

    // Time stepping for Semiconductor
    electrons.AssembleRHS(tCurrent,
                          tDelta);

    holes.AssembleRHS(tCurrent,
                      tDelta);

    reductants.AssembleRHS(tCurrent,
                           tDelta);

    oxidants.AssembleRHS(tCurrent,
                         tDelta);

    return 0;
}

int
ddpTimeStepping_type::
AssembleBackwardEulerSemiconductorRHS(ddpCarrier_type & electrons,
                                      ddpCarrier_type & holes,
                                      ddpCarrier_type & reductants,
                                      ddpCarrier_type & oxidants,
                                      const ddpGrid_type & semiconductor_grid,
                                      const double & tCurrent,
                                      const double & tDelta)
{
    ////////////////////////////////////////////////////////////////////
    // Compute the interface conditions
    ///////////////////////////////////////////////////////////////////

    // get interface density values
    electrons_interface = electrons.getInterfaceDensity();
    holes_interface = holes.getInterfaceDensity();
    reductants_interface = reductants.getInterfaceDensity();
    oxidants_interface = oxidants.getInterfaceDensity();


    // get the forward and backward transfer rates at this time step
    k_f = k_et *	(electrons_interface - electrons_e);
    k_b = k_ht *  (holes_interface - holes_e);

    // set the Robin boundary conditions at this time step
    electrons.setRight_RobinValue(k_f * oxidants_interface );
    holes.setRight_RobinValue( k_b * reductants_interface );
    ////////////////////////////////////////////////////////////////////
    // Update Recombination
    ///////////////////////////////////////////////////////////////////
//    grvy_timer_begin("Assemble Recombination Term");
    SRH_Recombination.updateRecombination(semiconductor_grid, electrons, holes);
//    grvy_timer_end("Assemble Recombination Term");

    ////////////////////////////////////////////////////////////////////
    // Assemble RHSes
    ///////////////////////////////////////////////////////////////////

    // Time stepping for Semiconductor
    electrons.AssembleRHS(tCurrent,
                          tDelta);

    holes.AssembleRHS(tCurrent,
                      tDelta);

    return 0;
}

int
ddpTimeStepping_type::
AssembleBackwardEulerElectrolyteRHS(ddpCarrier_type & electrons,
                                    ddpCarrier_type & holes,
                                    ddpCarrier_type & reductants,
                                    ddpCarrier_type & oxidants,
                                    const double & tCurrent,
                                    const double & tDelta)
{
    // get the forward and backward transfer rates at this time step
    k_f = k_et *	(electrons_interface - electrons_e);
    k_b = k_ht *  (holes_interface - holes_e);

    reductants.setLeft_RobinValue(
        (k_f * oxidants_interface - k_b * reductants_interface ));

    oxidants.setLeft_RobinValue(
        (k_b * reductants_interface - k_f *oxidants_interface) );


    ////////////////////////////////////////////////////////////////////
    // Assemble RHSes for just redox couples
    ///////////////////////////////////////////////////////////////////

    reductants.AssembleRHS(tCurrent,
                           tDelta);

    oxidants.AssembleRHS(tCurrent,
                         tDelta);

    return 0;
}

int
ddpTimeStepping_type::
UpdateImplicitDriftTerms(ddpCarrier_type & electrons,
                         ddpCarrier_type & holes,
                         ddpCarrier_type & reductants,
                         ddpCarrier_type & oxidants,
                         ddpPoisson_type & Poisson,
                         const ddpProblemInfo_type & problem,
                         const ddpGrid_type & grid)
{   // Get the electric field values for the carriers
    Poisson.getSemiconductorElecFieldDOFS(SemiconductorElecFieldDof);
    Poisson.getElectrolyteElecFieldDOFS(ElectrolyteElecFieldDof);

    electrons.UpdateImplicitDriftTerm(SemiconductorElecFieldDof, problem, grid);
    holes.UpdateImplicitDriftTerm(SemiconductorElecFieldDof, problem, grid);
    reductants.UpdateImplicitDriftTerm(ElectrolyteElecFieldDof, problem, grid);
    oxidants.UpdateImplicitDriftTerm(ElectrolyteElecFieldDof, problem, grid);

    return 0;
}

int
ddpTimeStepping_type::
UpdateExplicitDriftTerms(ddpCarrier_type & electrons,
                         ddpCarrier_type & holes,
                         ddpCarrier_type & reductants,
                         ddpCarrier_type & oxidants,
                         ddpPoisson_type & Poisson)

{   // Get the electric field values for the carriers
    Poisson.getSemiconductorElecFieldDOFS(SemiconductorElecFieldDof);
    Poisson.getElectrolyteElecFieldDOFS(ElectrolyteElecFieldDof);

    electrons.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);
    holes.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);
    reductants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);
    oxidants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);

    return 0;
}


int
ddpTimeStepping_type::
UpdateExplicitSemiconductorDriftTerms(ddpCarrier_type & electrons,
                                      ddpCarrier_type & holes,
                                      ddpPoisson_type & Poisson)

{
    // Get the electric field values for the carriers
    Poisson.getSemiconductorElecFieldDOFS(SemiconductorElecFieldDof);
    electrons.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);
    holes.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);
    return 0;
}


int
ddpTimeStepping_type::
UpdateExplicitElectrolyteDriftTerms(	ddpCarrier_type & reductants,
                                        ddpCarrier_type & oxidants,
                                        ddpPoisson_type & Poisson)

{   // Get the electric field values for the carriers
    Poisson.getElectrolyteElecFieldDOFS(ElectrolyteElecFieldDof);

    reductants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);
    oxidants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);

    return 0;
}

int
ddpTimeStepping_type::
BackwardEulerSolve(ddpCarrier_type & electrons,
                   ddpCarrier_type & holes,
                   ddpCarrier_type & reductants,
                   ddpCarrier_type & oxidants)
{
    electrons.BackwardEuler();
    holes.BackwardEuler();
    reductants.BackwardEuler();
    oxidants.BackwardEuler();

    return 0;
}

int
ddpTimeStepping_type::
BackwardEulerSolve(ddpCarrier_type & carrier1,
                   ddpCarrier_type & carrier2)
{
    carrier1.BackwardEuler();
    carrier2.BackwardEuler();

    return 0;
}

int
ddpTimeStepping_type::
IMIMEX_Solve(ddpCarrier_type & electrons,
             ddpCarrier_type & holes,
             ddpCarrier_type & reductants,
             ddpCarrier_type & oxidants,
             ddpPoisson_type & Poisson,
             const ddpGrid_type & semiconductor_grid,
             const ddpGrid_type & electrolyte_grid,
             const ddpProblemInfo_type & problem,
             const double & tCurrent,
             const double & tDelta)
{
    // Get the electric field values for the carriers
    Poisson.getSemiconductorElecFieldDOFS(SemiconductorElecFieldDof);
    Poisson.getElectrolyteElecFieldDOFS(ElectrolyteElecFieldDof);

    // set the drift terms using explicit data
    holes.UpdateImplicitDriftTerm(SemiconductorElecFieldDof,
                                  problem,
                                  semiconductor_grid);
    electrons.UpdateImplicitDriftTerm(SemiconductorElecFieldDof,
                                      problem,
                                      semiconductor_grid);


    //////////////////////////////////////////////////////////////////////////////
    // Semiconductor Solve
    //////////////////////////////////////////////////////////////////////////////

    // get interface density values
    reductants_interface = reductants.getInterfaceDensity();
    oxidants_interface = oxidants.getInterfaceDensity();

    double implicit_electron_value = k_et * oxidants_interface;// * electrons_interface;
    double explicit_electron_value = -k_et * electrons_e * oxidants_interface;
    double implicit_hole_value 		 = k_ht * reductants_interface;// * holes_interface;;
    double explicit_hole_vlaue 		 = -k_ht * holes_e * reductants_interface;

    // assemble matrices with implicit value
    electrons.AssembleLDGSystem(problem,
                                semiconductor_grid,
                                implicit_electron_value,
                                tDelta);

    holes.AssembleLDGSystem(problem,
                            semiconductor_grid,
                            implicit_hole_value,
                            tDelta);

    // set the Robin boundary conditions at this time step
    // using explicit values
    electrons.setRight_RobinValue(explicit_electron_value);
    holes.setRight_RobinValue(explicit_hole_vlaue);

    // Update recombination
 //   grvy_timer_begin("Assemble Recombination Term");
    SRH_Recombination.updateRecombination(semiconductor_grid, electrons, holes);
 //   grvy_timer_end("Assemble Recombination Term");

    // Set RHS
    electrons.AssembleRHS(tCurrent,
                          tDelta);


    holes.AssembleRHS(tCurrent,
                      tDelta);
    // solve and update
    electrons.BackwardEuler();
    holes.BackwardEuler();


    //////////////////////////////////////////////////////////////////////////////
    // Reductant Solve
    //////////////////////////////////////////////////////////////////////////////
    electrons_interface = electrons.getInterfaceDensity();
    holes_interface 		= holes.getInterfaceDensity();

    // get interface density values
    double implicit_reductant_value = k_ht * (holes_interface - holes_e);
    double explicit_reductant_value = -k_et * (electrons_interface - electrons_e)
                                      * oxidants_interface;
    // update the drift term
    reductants.UpdateImplicitDriftTerm(ElectrolyteElecFieldDof,
                                       problem,
                                       electrolyte_grid);

    // assemble matrices with implicit value
    reductants.AssembleLDGSystem(problem,
                                 electrolyte_grid,
                                 implicit_reductant_value,
                                 tDelta);

    // set the robin boundary condition using explicit values of rho_o, implicit rho_n
    reductants.setLeft_RobinValue(explicit_reductant_value);

    // assmble rhs
    reductants.AssembleRHS(tCurrent,
                           tDelta);

    // solve and update
    reductants.BackwardEuler();

    //////////////////////////////////////////////////////////////////////////////
    // Oxidant Solve
    //////////////////////////////////////////////////////////////////////////////

    // get interface density values
    reductants_interface = reductants.getInterfaceDensity();

    double implicit_oxidant_value = k_et * (electrons_interface - electrons_e);
    double explicit_oxidant_value =	-k_ht * (holes_interface - holes_e)
                                    * reductants_interface;

    // update the drift term
    oxidants.UpdateImplicitDriftTerm(ElectrolyteElecFieldDof,
                                     problem,
                                     electrolyte_grid);

    // assemble matrices with implicit value
    oxidants.AssembleLDGSystem(problem,
                               electrolyte_grid,
                               implicit_oxidant_value,
                               tDelta);

    // set the robin boundary condition using explicit values of rho_o, implicit rho_n
    oxidants.setLeft_RobinValue(explicit_oxidant_value);


    // assmble rhs
    oxidants.AssembleRHS(tCurrent,
                         tDelta);

    // solve and update
    oxidants.BackwardEuler();

    return 0;
}


int
ddpTimeStepping_type::
IMEXEX_Solve(ddpCarrier_type & electrons,
             ddpCarrier_type & holes,
             ddpCarrier_type & reductants,
             ddpCarrier_type & oxidants,
             ddpPoisson_type & Poisson,
             const ddpGrid_type & semiconductor_grid,
             const ddpGrid_type & electrolyte_grid,
             const ddpProblemInfo_type & problem,
             const double & tCurrent,
             const double & tDelta)

{
    // Get the electric field values for the carriers
    Poisson.getSemiconductorElecFieldDOFS(SemiconductorElecFieldDof);
    Poisson.getElectrolyteElecFieldDOFS(ElectrolyteElecFieldDof);

    // set the drift terms using explicit data
    holes.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);
    electrons.UpdateExplicitDriftTerm(SemiconductorElecFieldDof);

    //////////////////////////////////////////////////////////////////////////////
    // Semiconductor Solve
    //////////////////////////////////////////////////////////////////////////////

    // get interface density values
    reductants_interface = reductants.getInterfaceDensity();
    oxidants_interface = oxidants.getInterfaceDensity();

    double implicit_electron_value = k_et * oxidants_interface;// * electrons_interface;
    double explicit_electron_value = -k_et * electrons_e * oxidants_interface;
    double implicit_hole_value 		 = k_ht * reductants_interface;// * holes_interface;;
    double explicit_hole_vlaue 		 = -k_ht * holes_e * reductants_interface;

    // assemble matrices with implicit value
    electrons.AssembleLDGSystem(problem,
                                semiconductor_grid,
                                implicit_electron_value,
                                tDelta);

    holes.AssembleLDGSystem(problem,
                            semiconductor_grid,
                            implicit_hole_value,
                            tDelta);

    // set the Robin boundary conditions at this time step
    // using explicit values
    electrons.setRight_RobinValue(explicit_electron_value);
    holes.setRight_RobinValue(explicit_hole_vlaue);

    // Update recombination
  //  grvy_timer_begin("Assemble Recombination Term");
    SRH_Recombination.updateRecombination(semiconductor_grid, electrons, holes);
  //  grvy_timer_end("Assemble Recombination Term");

    // Set RHS
    electrons.AssembleRHS(tCurrent,
                          tDelta);


    holes.AssembleRHS(tCurrent,
                      tDelta);

    // solve and update
    electrons.BackwardEuler();
    holes.BackwardEuler();


    //////////////////////////////////////////////////////////////////////////////
    // Reductant Solve
    //////////////////////////////////////////////////////////////////////////////
    electrons_interface = electrons.getInterfaceDensity();
    holes_interface 		= holes.getInterfaceDensity();

    // get interface density values
    double implicit_reductant_value = k_ht * (holes_interface - holes_e);
    double explicit_reductant_value = -k_et * (electrons_interface - electrons_e)
                                      * oxidants_interface;

    // assemble matrices with implicit value
    reductants.AssembleLDGSystem(problem,
                                 electrolyte_grid,
                                 implicit_reductant_value,
                                 tDelta);

    // set the robin boundary condition using explicit values of rho_o, implicit rho_n
    reductants.setLeft_RobinValue(explicit_reductant_value);

    // update the drift term
//	reductants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);

    // assmble rhs
    reductants.AssembleRHS(tCurrent,
                           tDelta);

    // solve and update
    reductants.BackwardEuler();

    //////////////////////////////////////////////////////////////////////////////
    // Oxidant Solve
    //////////////////////////////////////////////////////////////////////////////

    // get interface density values
    reductants_interface = reductants.getInterfaceDensity();

    double implicit_oxidant_value = k_et * (electrons_interface - electrons_e);
    double explicit_oxidant_value =	-k_ht * (holes_interface - holes_e)
                                    * reductants_interface;

    // assemble matrices with implicit value
    oxidants.AssembleLDGSystem(problem,
                               electrolyte_grid,
                               implicit_oxidant_value,
                               tDelta);

    // set the robin boundary condition using explicit values of rho_o, implicit rho_n
    oxidants.setLeft_RobinValue(explicit_oxidant_value);

    // update the drift term
//	oxidants.UpdateExplicitDriftTerm(ElectrolyteElecFieldDof);

    // assmble rhs
    oxidants.AssembleRHS(tCurrent,
                         tDelta);

    // solve and update
    oxidants.BackwardEuler();

    return 0;
}

bool
ddpTimeStepping_type::
check_converged(ddpCarrier_type & electrons,
                ddpCarrier_type & holes,
                const double & tol)
{
    electron_q_values =  electrons.VandeMondeMatrices.globalVandeMondeDG *
                         electrons.carrierState.qDof;

    hole_q_values = holes.VandeMondeMatrices.globalVandeMondeDG *
                    holes.carrierState.qDof;

    q_values = electron_q_values - hole_q_values;

    double current_value = q_values(0);

    for(int i = 1; i < q_values.size(); i++)
    {
        if( fabs(current_value - q_values(i)) > tol)
            return false;
    }

    // now its true
    left_value = current_value;
    return true;
}
