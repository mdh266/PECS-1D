#include "../include/testUtilities.hpp"

// ZeroFunction is delcared in Utilities.cpp
// so it can be used in carrier.cpp
// no longer needed here


// Functions that will be used in test_System.

double MinusSineFunction(const double & x)
{
    return -sin(x);
}

double MinusCosineFunction(const double & x)
{
    return -cos(x);
}
double SineFunction(const double & x)
{
    return sin(x);
}

double Sine2Function(double const & x)
{
    return sin(2*x);
}

double Sine2PiFunction(double const & x)
{
    return sin(2*M_PI*x);
}

double CosineFunction(const double & x)
{
    return cos(x);
}

double CosinePiFunction(const double & x)
{
    return cos(M_PI*x);
}

double Cosine2Function(double const & x)
{
    return cos(2*x);
}

double Cosine2PiFunction(double const & x)
{
    return cos(2*M_PI*x);
}


double PiFunction(const double & x)
{
    return M_PI;
}

double ThirtyFunction(const double & x)
{
    return 30.0;
}

double TwoFunction(const double & x)
{
    return 2.0;
}

double MinusPiFunction(const double & x)
{
    return -M_PI;
}

double OneFunction(const double & x)
{
    return 1.0;
}

double MinusOneFunction(const double & x)
{
    return -1.0;
}

double OneOverPiLineFunction(double const & x)
{
    return (1.0/M_PI) * x;
}

double ExponentialFunction(const double & x)
{
    return exp(x);
}


double MinusOneOverPiFunction(double const & x)
{
    return -1.0/M_PI;
}

double PotentialInterfaceZeroRHS(double const & x)
{
    if(x < 0.0)
        return (2.0/3.0) * (x + 1.0);
    else
        return (1.0/3.0) * (x + 2.0);
}

double DisplElecFieldInterfaceZeroRHS(double const & x)
{
    return (-2.0/3.0)*M_PI;
}

double PotentialInterfaceNonZeroRHS1(double const & x)
{
    double Phi_C = 2.0/3.0 + 1.0/(6.0*M_PI);

    if(x < 0.0)
        return -(0.5/M_PI)*x*x + ((2.0/3.0 - 1.0/(3.0*M_PI))*x) + Phi_C;
    else
        return ( 1.0/3.0 -1.0/(6.0*M_PI))*x + Phi_C;
}

double DisplElecFieldInterfaceNonZeroRHS1(double const & x)
{
    if(x < 0.0)
        return x - (2.0/3.0)*M_PI + 1.0/3.0;
    else
        return -(2.0/3.0)*M_PI + 1.0/3.0;
}

double PotentialInterfaceNonZeroRHS2(double const & x)
{
    double Phi_C = 2.0/3.0 + 1.0/(6.0*M_PI);

    if(x < 0.0)
        return (2.0/3.0 + 1.0/(6.0*M_PI) )*x + Phi_C;
    else
        return -(1.0/(4*M_PI) )*x*x + ( 1.0/3.0 + 1.0/(12.0*M_PI) )*x + Phi_C;
}

double DisplElecFieldInterfaceNonZeroRHS2(double const & x)
{
    double Phi_C = -(2.0*M_PI)/3.0 - 1.0/(6.0);

    if(x < 0.0)
        return Phi_C;
    else
        return x + Phi_C;
}


double Doping_N_A(double const & x)
{
    if(x <= 0.5)
        return 0.0;
    else
        return 1.0;
}

double Doping_N_D(double const & x)
{
    if(x <= 0.5)
        return 1.0;
    else
        return 0.0;
}

double Doping_PN(double const & x)
{
    if(x <= 0.5)
        return 1.0;
    else
        return -1.0;
}

double True_N_A(double const & x)
{
    if(x < 0.5)
        return 0.0;
    else
        return -1.0;
}

double DD_L2_RHS(double const & x)
{
    double N = 5;
    return -N*(N-1) * std::pow(x,N-2); // + N * std::pow(x,N-1);
}

double DD_L2_Density(double const & x)
{
    double N = 5;
    return std::pow(x,N);
}

double  DD_L2_Current(double const & x)
{
    double N = 5;
    return -(N) * std::pow(x, N-1); // + std::pow(x,N);
}

double DD_L2_BC(double const & x)
{
    return DD_L2_Current(1);
}

double Poisson_L2_RHS(double const & x)
{
    return -cos(x);
//	double N = 5;
//  return -N * N*(N-1) * std::pow(x,N-2);
}

double Poisson_L2_Potential(double const & x)
{
    return cos(x);
//	double N = 5;
//	return N*std::pow(x,N) - (N-1) * x;
}

double  Poisson_L2_ElecField(double const & x)
{
    return -sin(x);
//	double N = 5;
//	return -N * N * std::pow(x,N-1) + (N-1);
}

// True solutions are evaluated at T = 1.0

// Solution to diffusion equation with homogenous boundary
// conditions

double DiffusionSolutionHomoU(double const & x)
{
    return  exp(-1.0)*sin(x);
}

double DiffusionSolutionHomoQ(double const & x)
{
    // q = -du/dx
    return -exp(-1.0)*cos(x);
}

// The boundary condition functions to test
// diffusion equation with non-homogenous boundary conditions

double DiffusionNonHomoLeftBC(double const & t)
{
    return exp(-t);
}

double DiffusionNonHomoRightBC(double const & t)
{
    return -exp(-t);
}


// Solution to diffusion equation with non-homogenous boundary
// conditions

double DiffusionSolutionNonHomoU(double const & x)
{
    return  exp(-1.0)*cos(x);
}

double DiffusionSolutionNonHomoQ(double const & x)
{
    // q = -du/dx
    return exp(-1.0)*sin(x);
}

double DiffusionNonHomoLeftBC_Scaling(double const & t)
{
    // u at x = 0
    double tau = 2.0;
    double C 	 = 1.0;
    double L 	 = M_PI; // length scale
    double mu  = 1.0;
    return /*(C*tau /(L*L)) */ exp(-t * tau/mu);
}

double DiffusionNonHomoRightBC_Scaling(double const & t)
{
    // u at x = \Pi
    double tau = 2.0;	 // time scale
    double C 	 = 5.0;
    double L 	 = M_PI; // length scale
    double mu  = 1.0;
    return - /*(tau /(L*L)) */ exp(-t * tau /mu );
}

double DiffusionSolutionNonHomoU_Scaling(double const & x)
{
    // u at t = 1
    double tau = 2.0;	 // time scale
    double L 	=  M_PI; // length scale
    double mu  = 1.0;  // Mobility
    double C = 5.0; // doping scale
    return  C*exp(-1.0 * tau / mu)*cos(L*x);
}

double DiffusionSolutionNonHomoQ_Scaling(double const & x)
{
    double tau = 2.0;	 // time scale
    double L 	=  M_PI; // length scale
    double mu  =  1.0;  // Mobility
    double C = 5.0; // doping scale
    // q = -du/dx at t = 1
    return C*exp(-1.0 * tau / mu)*sin(L*x);
}
// Solution to diffusion equation with Neumann boundary
// conditions

double DiffusionSolutionNeumannHomoU(double const & x)
{
    return  exp(-4)*cos(2*x);
}

double DiffusionSolutionNeumannHomoQ(double const & x)
{
    // q = -du/dx
    return 2*exp(-4.0)*sin(2*x);
}

// The boundary condition functions to test
// diffusion equation with Neumann boundary conditions

double DiffusionNeumannNonHomoBC(double const & t)
{
    return -2*exp(-4*t);
}

// Solution to diffusion equation with Neumann boundary
// conditions

double DiffusionSolutionNeumannNonHomoU(double const & x)
{
    return  exp(-4)*sin(2*x);
}

double DiffusionSolutionNeumannNonHomoQ(double const & x)
{
    // q = -du/dx
    return -2*exp(-4.0)*cos(2*x);
}

double DiffusionSolutionNeumannHomoU_Scaling(double const & x)
{
    double C 	 = 5.0;
    double L 	 = M_PI; // length scale
    double mu  = 1.0;
    double tau = 2.0;

    return  C*exp(-4*tau)*cos(2*x*L);
}

double DiffusionSolutionNeumannHomoQ_Scaling(double const & x)
{
    double C 	 = 5.0;
    double L 	 = M_PI; // length scale
    double mu  = 1.0;
    double tau = 2.0;

    // q = -du/dx
    return 2*C*exp(-4.0*tau)*sin(2*x*L);
}

double DiffusionWithScalingInitial(const double & x)
{
    return Sine2PiFunction(x);
}

double DiffusionWithScalingAndNeumannRightBC(double const & t)
{
    double C = 5.0;
    double L = M_PI;
    double tau = 2.0;

    return (tau/L) * (-2*exp(-4*t*tau) );
}

double DiffusionWithScalingU(double const & x)
{
    double C = 5.0;
    double L = M_PI;
    double tau = 2.0;
    return  C *exp(-4*tau)*sin(2*L*x);
}

double DiffusionWithScalingQ(double const & x)
{
    double C = 5.0;
    double L = M_PI;
    double tau = 2.0;

    return -2*C*exp(-4.0*tau)*cos(2*L*x);
}



/// DRIFT DIFFUSION EQUATION
//  HOMOGENOUS DIRICHLET BC


// ELECTRONS Initial conditions and true solutions

double EDD_Dirichlet_Homo_IC(double const & x)
{
    // Drift Diffusion electrons initial conditions

    return exp(-0.5*x) * sin(x);
}

double EDD_Dirichlet_Homo_SolutionU(double const & x)
{
    // Drift Diffusion solution for electrons solution at T = 1.0
    return exp(-1.25) * exp(-0.5*x) * sin(x);
}

double EDD_Dirichlet_Homo_SolutionnQ(double const & x)
{
    // q = -u -du/dx
    // derivative of solution to Drift Diffusion for
    // electrons solution at T = 1.0
    double Eu = EDD_Dirichlet_Homo_SolutionU(x);
    double du = exp(-1.25) * exp(-0.5*x) * ( cos(x) - 0.5*sin(x) );

    return -Eu - du;
}

// HOLES Initial conditions and true solutions

double HDD_Dirichlet_Homo_IC(double const & x)
{
    // Drift Diffusion electrons initial conditions

    return exp(0.5*x) * sin(x);
}

double HDD_Dirichlet_Homo_SolutionU(double const & x)
{
    // Drift Diffusion solution for electrons solution at T = 1.0
    return exp(-1.25) * exp(0.5*x) * sin(x);
}

double HDD_Dirichlet_Homo_SolutionnQ(double const & x)
{
    // q = -du/dx
    // derivative of solution to Drift Diffusion for
    // electrons solution at T = 1.0

    double du = exp(-1.25) * exp(0.5*x) * ( cos(x) + 0.5*sin(x) );
    double Eu = HDD_Dirichlet_Homo_SolutionU(x);
    return Eu - du;
}


// NON HOMOGENOUS DIRICHLET BC

//ELECTRONS
double EDD_Dirichlet_Non_Homo_IC(double const & x)
{
    // Drift Diffusion electrons initial conditions

    return exp(-0.5*x) * cos(x);
}

double EDD_Dirichlet_Non_Homo_LeftBC(double const & t)
{
    return exp(-1.25*t);
}

double EDD_Dirichlet_Non_Homo_RightBC(double const & t)
{
    return -exp(-1.25*t)*exp(-0.5*M_PI);
}

double EDD_Dirichlet_Non_Homo_SolutionU(double const & x)
{
    // Drift Diffusion solution for electrons solution at T = 1.0
    return exp(-1.25) * exp(-0.5*x) * cos(x);
}

double EDD_Dirichlet_Non_Homo_SolutionnQ(double const & x)
{
    // q = -du/dx
    // derivative of solution to Drift Diffusion for
    // electrons solution at T = 1.0
    double Eu = EDD_Dirichlet_Non_Homo_SolutionU(x);
    double du =  - exp(-1.25) * exp(-0.5*x) * ( 0.5*cos(x) + sin(x) );
    return -Eu - du;
}

// HOLES


double HDD_Dirichlet_Non_Homo_IC(double const & x)
{
    // Drift Diffusion electrons initial conditions

    return exp(+0.5*x) * cos(x);
}

double HDD_Dirichlet_Non_Homo_LeftBC(double const & t)
{
    return exp(-1.25*t);
}

double HDD_Dirichlet_Non_Homo_RightBC(double const & t)
{
    return -exp(-1.25*t)*exp(0.5*M_PI);
}

double HDD_Dirichlet_Non_Homo_SolutionU(double const & x)
{
    // Drift Diffusion solution for electrons solution at T = 1.0
    return exp(-1.25) * exp(0.5*x) * cos(x);
}

double HDD_Dirichlet_Non_Homo_SolutionnQ(double const & x)
{
    // q = -du/dx
    // derivative of solution to Drift Diffusion for
    // electrons solution at T = 1.0
    double Eu = HDD_Dirichlet_Non_Homo_SolutionU(x);
    double du = exp(-1.25) * exp(0.5*x) * ( 0.5*cos(x) - sin(x) );
    return Eu-du;

}
//////////////////////////////////////////////////////////////////////////////
// Mixed BOUNDARY CONDITIONS
//////////////////////////////////////////////////////////////////////////////



//LEFT HOMOGENOUS DIRICHLET , ROBIN RIGHT
double EDD_Dir_Homo_Robin_Non_Homo_RightBC(double const & t)
{
    return +exp(-1.25*t)*exp(-0.5*M_PI);
}


// LEFT NON-HOMOGENOUS DIRICHLET, RIGHT NON-HOMOGENOUS ROBIN
double EDD_Dir_Non_Homo_Robin_Non_Homo_RightBC(double const & t)
{
    return 0.5*exp(-1.25*t)*exp(-0.5*M_PI);
}

// LEFT NEUMANN NON-HOMOGENOUS, RIGHT NON-HOMOGENOUS DIRICHLET
double EDD_Robin_Non_Homo_Dir_Non_Homo_LeftBC(double const & t)
{
    return 1.5*exp(-1.25*t);
}

// LEFT ROBIN NON-HOMOGENOUS, RIGHT HOMOGENOUS DIRICHLET
double EDD_Robin_Non_Homo_Dir_Homo_LeftBC(double const & t)
{
    return -exp(-1.25*t);
}

// HOLES

// LEFT NEUMANN
double HDD_Non_Homo_Neumann_Homo_Dir_LeftBC(double const & t)
{
    return -exp(-1.25*t);
}

// LEFT ROBIN
double HDD_Non_Homo_Robin_Non_Homo_Dir_LeftBC(double const & t)
{
    return 0.5*exp(-1.25*t);
}

// RIGHT NEUMANN
double HDD_Homo_Dir_Non_Homo_Neumann_RightBC(double const & t)
{
    return exp(0.5*M_PI) * exp(-1.25*t);
}

//RIGHT ROBIN
double HDD_Homo_Dir_Non_Homo_Robin_RightBC(double const & t)
{
    return (5.0/(4.0*sqrt(2.0)) ) * exp(0.5*M_PI) * exp(-13.0*t/16.0);
}

double HDD_Right_Robin_IC(double const & x)
{
    return exp(0.5*x)*sin(0.75*x);
}

double HDD_Right_Robin_Solution_U(double const & x)
{
    return exp(0.5*x - 13.0/16.0)*sin(0.75*x);
}

double HDD_Right_Robin_Solution_Q(double const & x)
{
    double Eu = HDD_Right_Robin_Solution_U(x);
    double du = 0.25*exp(0.5*x - 13.0/16.0)*(3.0*cos(0.75*x) + 2.0*sin(0.75*x));
    return Eu - du;
}

///////////////////////////////////////////////////////////////////////
// STEADY STATE SOLUTIONS
///////////////////////////////////////////////////////////////////////

double HDD_SteadyState_Dir_U(double const & x)
{
    double C = 1.0/(1.0 - exp(-1.0) );
    return C * ( 1.0 - exp( -x / M_PI) );
}

double HDD_SteadyState_Dir_Q(double const & x)
{
    double C = 1.0/(1.0 - exp(-1.0) );
    return -1.0/M_PI * C;
}

double Diffusion_SteadyState_Robin_U(double const & x)
{
    return 1.0 - M_PI*x;
}

double HDD_SteadyState_Robin_U(double const & x)
{
    return M_PI * ( exp( -(x - M_PI)/M_PI ) - 1.0 );
}

////////////////////////////////////////////////////////////////////////
// REACTIVE FLUXES WITH FORWARD EULER
///////////////////////////////////////////////////////////////////////
double ReactiveFluxesU(double const & x)
{
    return 1.0 - ( x / M_PI);
}

double ReactiveFluxesQ(double const & x)
{
    return 0.15816;
}

// tests whether the numerical solution and true solution are the same
// within a numerical tolerance.  Tests L_{\infinty} convergence.
bool Are_Same(ddpCarrier_type const & carrier1,
              ddpCarrier_type const & carrier2)
{
    double tol = 1.0e-2;

    double u1, u2, q1, q2;

    if(carrier1.carrierState.uDof.size() != carrier2.carrierState.uDof.size())
    {
        cout << "uDofs size not same" << endl;
        return false;
    }

    if(carrier1.carrierState.qDof.size() != carrier2.carrierState.qDof.size())
    {
        cout << "qDofs size not same" << endl;
        return false;

    }

    ddpDenseVector_type
    U1PTS =
        carrier1.VandeMondeMatrices.globalVandeMondeDG
        * carrier1.carrierState.uDof,
        Q1PTS =
            carrier1.VandeMondeMatrices.globalVandeMondeDG
            * carrier1.carrierState.qDof,
            U2PTS =
                carrier2.VandeMondeMatrices.globalVandeMondeDG
                * carrier2.carrierState.uDof,
                Q2PTS =
                    carrier2.VandeMondeMatrices.globalVandeMondeDG
                    * carrier2.carrierState.qDof;


    for(int i = 0; i < U1PTS.size(); i++)
    {
        u1 = U1PTS(i);
        u2 = U2PTS(i);
        q1 = Q1PTS(i);
        q2 = Q2PTS(i);

        if(!(fabs(u1) >= 0))
        {
            cout << "u1 = nan" << endl;
            return false;
        }
        if(!(fabs(q1) >= 0))
        {
            cout << "q1 = nan" << endl;
            return false;
        }
        if( fabs( u1 - u2  ) > tol )
        {
            cout << "U values not the same at " << i << endl;
            cout << "u1 = " << u1 << endl;
            cout << "u2 = " << u2 << endl;
            return false;
        }

        if( fabs( q1 - q2 ) > tol)
        {
            cout << "Q values not the same at " << i << endl;
            cout << "q1 = " << q1 << endl;
            cout << "q2 = " << q2 << endl;
            return false;
        }

    }

    return true;
}


// tests whether the numerical solution and true solution are the same
// within a numerical tolerance.  Tests L_{\infinty} convergence.
bool Are_Same_Current(ddpCarrier_type const & carrier1,
                      ddpCarrier_type const & carrier2)
{
    double tol = 1.0e-4;

    double u1, u2, q1, q2;

    if(carrier1.carrierState.qDof.size() != carrier2.carrierState.qDof.size())
    {
        cout << "qDofs size not same" << endl;
        return false;

    }

    ddpDenseVector_type
    Q1PTS =
        carrier1.VandeMondeMatrices.globalVandeMondeDG
        * carrier1.carrierState.qDof,
        Q2PTS =
            carrier2.VandeMondeMatrices.globalVandeMondeDG
            * carrier2.carrierState.qDof;


    for(int i = 0; i < Q1PTS.size(); i++)
    {
        q1 = Q1PTS(i);
        q2 = Q2PTS(i);

        if( fabs( q1 - q2 ) > tol)
        {
            cout << "Q values not the same at " << i << endl;
            cout << "q1 = " << q1 << endl;
            cout << "q2 = " << q2 << endl;
            return false;
        }

    }

    return true;
}

bool Dofs_Are_Same(ddpPoisson_type const & P1, ddpPoisson_type const & P2)
{
    double tol = 1.0e-8;

    if(P2.PoissonState.elecDof.size() != P1.PoissonState.elecDof.size())
    {
        cout << "elecDofs size not same" << endl;
        return false;
    }

    if(P2.PoissonState.potDof.size() != P1.PoissonState.potDof.size())
    {
        cout << "potDofs size not same" << endl;
        return false;
    }

    for(int i = 0; i < P1.PoissonState.potDof.size(); i++)
    {
        if( fabs(P1.PoissonState.potDof(i)-P2.PoissonState.potDof(i)) > tol )
        {
            cout << "potDofs not same at entry " << i << endl;
            return false;
        }
        if( fabs(P1.PoissonState.elecDof(i) - P2.PoissonState.elecDof(i)) > tol )
        {
            cout << "elecDofs not same at entry " << i << endl;
            return false;
        }

    }
    return true;
}

bool Are_Same(ddpDenseVector_type const & DP1, ddpDenseVector_type const & DP2)
{
    double tol = 1.0e-8;

    if(DP2.size() != DP1.size())
    {
        cout << "PTValues vector size not same" << endl;
        return false;
    }

    for(int i = 0; i < DP1.size(); i++)
    {
        if( fabs(DP1(i)-DP2(i)) > tol )
        {
            cout << "PTValues different at entry " << i << endl;
            return false;
        }
    }

    return true;
}

int L2_Error(ddpPoisson_type const & TestP,
             double (*TrueU)(const double & x),
             double (*TrueQ)(const double & x),
             ddpDenseVector_type & errors)
{
    double PotentialError;
    double ElectricFieldError;

    ddpDenseVector_type Pot1PTVals = TestP.VandeMondeMatrices.globalVandeMondeDG
                                     * TestP.PoissonState.potDof;

    ddpDenseVector_type Elec1PTVals = TestP.VandeMondeMatrices.globalVandeMondeMX
                                      * TestP.PoissonState.elecDof;

    ddpDenseVector_type Weights = TestP.weightsDense;
    int n = Weights.size();
    double x = 0;

    for(int i = 0; i < n; i++)
    {
        x = TestP.PTSDense(i);
        PotentialError += (Pot1PTVals(i) - TrueU(x))
                          * (Pot1PTVals(i) - TrueU(x) ) * Weights(i);

        ElectricFieldError += (Elec1PTVals(i) - TrueQ(x) )
                              * (Elec1PTVals(i) - TrueQ(x)) * Weights(i);
    }

    PotentialError = sqrt(PotentialError);
    ElectricFieldError = sqrt(ElectricFieldError);

    errors(0) = PotentialError;
    errors(1) = ElectricFieldError;


//	cout << "PotentialError = " << PotentialError << endl;
//	cout << "ElectricFieldError = " << ElectricFieldError << endl;

    return 0;
}
int L2_Error(ddpCarrier_type const & TestP,
             double (*TrueU)(const double & x),
             double (*TrueQ)(const double & x),
             ddpDenseVector_type & errors)
{
    double PotentialError;
    double ElectricFieldError;


    ddpDenseVector_type Pot1PTVals = TestP.VandeMondeMatrices.globalVandeMondeDG
                                     * TestP.carrierState.uDof;

    ddpDenseVector_type Elec1PTVals = TestP.VandeMondeMatrices.globalVandeMondeDG
                                      * TestP.carrierState.qDof;

    ddpDenseVector_type Weights = TestP.weightsDense;
    int n = Weights.size();
    double x = 0;

    for(int i = 0; i < n; i++)
    {
        x = TestP.PTSDense(i);
        PotentialError += (Pot1PTVals(i) - TrueU(x))
                          * (Pot1PTVals(i) - TrueU(x) ) * Weights(i);

        ElectricFieldError += (Elec1PTVals(i) - TrueQ(x) )
                              * (Elec1PTVals(i) - TrueQ(x)) * Weights(i);
    }

    PotentialError = sqrt(PotentialError);
    ElectricFieldError = sqrt(ElectricFieldError);

    errors(0) = PotentialError;
    errors(1) = ElectricFieldError;
//	cout << "PotentialError = " << PotentialError << endl;
//	cout << "ElectricFieldError = " << ElectricFieldError << endl;

    return 0;
}




// used for printing the point values of soultions out into a format
// that can be read and converted into a movie via python script
int print2File(ddpCarrier_type const & testCarrier,
               ddpCarrier_type const & trueCarrier,
               ddpPoisson_type const & testPoisson,
               int testNum)
{
    ofstream prt;
    std::string prefix = "MovieCarrier";
    std::string extension = ".dat";
    std::string filename;
    std::stringstream ss;

    ss << prefix << testNum << extension << '\0';
    filename = ss.str();
    char str[filename.size()];

    for(unsigned int i = 0; i < filename.size(); i++)
        str[i] = filename[i];

    prt.open(str);

    ddpDenseVector_type
    UPTS_test = testCarrier.VandeMondeMatrices.globalVandeMondeDG
                * testCarrier.carrierState.uDof,
                UPTS_true = trueCarrier.VandeMondeMatrices.globalVandeMondeDG
                            * trueCarrier.carrierState.uDof;

    ddpDenseVector_type
    QPTS_test = testCarrier.VandeMondeMatrices.globalVandeMondeDG
                * testCarrier.carrierState.qDof,
                QPTS_true = trueCarrier.VandeMondeMatrices.globalVandeMondeDG
                            * trueCarrier.carrierState.qDof;

    ddpDenseVector_type
    ElecPTS = testPoisson.VandeMondeMatrices.globalVandeMondeMX
              * testPoisson.PoissonState.elecDof,
              POTPTS = testPoisson.VandeMondeMatrices.globalVandeMondeDG
                       * testPoisson.PoissonState.potDof;


    for(int i = 0; i < UPTS_test.size(); i++)
    {
        prt << testCarrier.PTSDense(i) << " "
            << UPTS_test(i) << " "
            << QPTS_test(i) << " "
            << UPTS_true(i) << " "
            << QPTS_true(i) << " "
            << ElecPTS(i) << " "
            << POTPTS(i) << endl;
    }

    prt.close();
    return 0;
}



