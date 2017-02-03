///	\mainpage  

\section Overview

\subsection Introduction
This software is designed to numerically simulate a coupled system of
reaction diffusion equations that describe the macroscopic dynamics of 
charge transport in photoelectrochemical (PEC) solar cell.
The main objective is to accurately capture the reactive dynamics of the 
semiconductor-electrolyte interface. The underlying spatial descritizations
are based on the local discontinuous Galerkin on each carrier (ddpCarrier_type) 
and mixed finite element method on Poisson's equation (ddpPoisson_type).  The use of specific tailored implicit-explicit (IMEX) methods to capture transient solutions of the
time dependent nonlinear systems of PDES (ddpTimeStepping_type).


\subsection Requirements 
The dependency requirements for this software are, \n
1.) <a href="https://eigen.tuxfamily.org"> Eigen sparse linear algebra library </a> \n
2.) <a href="https://www.boost.org"> boost C++ library </a> \n
3.) <a href="https://www.gnu.org/software/gsl"> GNU Scientific Library </a> \n
4.) <a href="https://www.cmake.org"> CMake </a> \n
5.) <a href="https://www.openmp.org"> OpenMP </a>  (Optional) 


\subsection Model
We focus on the reactive dynamics of the semiconductor-electrolyte interface.
Therfore our domain is the half cell, whose abstract representation is, 
\image html DomainDecomposition.png
	

In the semiconductor component (\f$ \Omega_{S} \f$) the transport of electrons 
\f$(\rho_{n})\f$ and holes \f$(\rho_{p})\f$ is governed by 
the drift-diffusion-Poisson system of equations for electrons
\f[ \begin{align} 
\frac{\partial \rho_{n}}{\partial t}
\, + \, 
\boldsymbol \nabla  \cdot  \,\left( \, \mu_{n}  \boldsymbol \nabla \Phi \, \rho_{n} 
\, - \,
D_{n} \  \boldsymbol \nabla \, \rho_{n}  \, \right) 
\; &= \; 
-R(\rho_{n}, \rho_{p})
\, + \, 
G && \text{in} \ (0,T] \ \times \ \Omega_{S} \nonumber \\
\frac{\partial \rho_{p}}{\partial t} 
\, + \, 
\boldsymbol \nabla \cdot  \, \left( - \mu_{p} \boldsymbol \nabla \Phi \, \rho_{p} 
\, - \, 
D_{p} \  \boldsymbol \nabla \, \rho_{p} \, \right) 
\; &= \;
-R(\rho_{n}, \rho_{p})
\, + \, 
G
&& \text{in} \ (0,T] \ \times \ \Omega_{S} \\
 -\boldsymbol \nabla \cdot \left( \, \epsilon_{r}^{S} \, \boldsymbol \nabla \Phi \right)
\; &= \;  
\frac{q}{\epsilon_{0}} \left[ 
\rho_{p}^{e} - 
\rho_{n}^{e}  - (\rho_{n} -\rho_{p}) \right] 
  && \text{in} \ (0,T] \ \times \ \Omega_{S}  \nonumber .
 \end{align} \f]


Where \f$\mu_{i}\f$ and \f$D_{i}\f$ are the mobility and diffusivity of carrier \f$i = n,p\f$
respectively.  The functions \f$\rho_{n}^{e}\f$
and \f$\rho_{p}^{e}\f$ are the equilibrium electron and hole densities respectively.  
The constants
\f$q\f$ and \f$\epsilon_{0}\f$ charge of the electron and permittivity of free space. 
The material permittivity is \f$\epsilon_{r}\f$ is assumed to be an invertible tensor.


We use Shockley-Reed-Hall recombination as our sink functional,
\f[ \begin{equation}
R(\rho_{n}, \rho_{p})
\; = \;
\frac{ \rho_{n} \, \rho_{p} \, - \, \rho_{i}^{2} }
{\tau_{n} \, (\rho_{n} \, + \, \rho_{i}) \, + \, \tau_{p} \, ( \rho_{p} \, + \, \rho_{i})}.
\end{equation}
\f]


The term \f$\rho_{i}\f$ is the intrinsic electron density and can be spatially varying. 
\f$\tau_{n}, \ \tau_{p} \f$ are constants called the electron and hole lifetimes.  
The generation of electrons and holes is modeled using a macroscopic source function,



\f[
\begin{equation}
G(\textbf{x}) \; = \; \left\{ 
\begin{array}{lr}
\alpha (\textbf{x}) 
\, G_{0} \, e^{- \, \int_{0}^{s} \, \alpha 
(\, \textbf{x}_{0} \, + \, s' \, \boldsymbol \theta_{0} \, ) \, ds'} 
\qquad & \qquad \text{if} \quad  \textbf{x} 
\; = \; \textbf{x}_{0} \, + \, s \, \boldsymbol \theta_{0} \\
0 \qquad & \qquad \text{otherwise}
\end{array} \right. 
\end{equation} \f]



The point \f$\textbf{x}_{0}\f$ is the photon's incident location and 
\f$\boldsymbol \theta_{0}\f$ is the incident direction. The absorption coefficient 
\f$\alpha(\textbf{x}) \f$ has been averaged over all 
energy values of light that generate free carriers.  The term 
\f$G(\textbf{x}_{0}) \; [ \, \text{cm}^{-2} \, \text{s}^{-1} \, ]\f$ 
represents the surface photon flux at the point \f$\textbf{x}_{0}\f$. 


The portion of the boundary of the semiconductor \f$\Gamma_{S}\f$ is 
an Ohmic metal contact where the charge densities take on their equilibrium values,

\f[ \begin{align}
\rho_{n}   \; &= \;  
\rho_{n}^{e}  && \text{on} \ (0,T] \ \times \ \Gamma_{S}  \\
\rho_{p} \; &= \; 
\rho_{p}^{e}  && \text{on} \ (0,T] \ \times \ \Gamma_{S}
\end{align} 
\f]


The potential on an Ohmic contact is the sum of the applied voltage 
\f$\Phi_{\text{app}} \f$ to and the so-called ``built-in'' potential 
\f$\Phi_{\text{bi}}\f$ ,

\f[
\begin{align} 
\Phi  \;  = \; \Phi_{\text{bi}} \, + \, \Phi_{\text{app.}}  
&& \text{on} \ (0,T] \ \times \ \Gamma_{S}.
\end{align}
\f]

In the electrolyte component \f$(\Omega_{E})\f$ the transport of reductants 
(\f$ \rho_{r} \f$) and oxidants (\f$ \rho_{o} \f$ ) is governed by a 
similar drift-diffusion-Poisson system,

\f[
\begin{align}
\frac{\partial \rho_{r}}{\partial t}
\, + \, 
\boldsymbol \nabla  \cdot \, \left(  \mu_{r}  \, \boldsymbol \nabla \Phi \, \rho_{r} 
\, - \, D_{r} \,
 \boldsymbol \nabla \, \rho_{r}  \, \right) 
\; &= \; 0 , && \text{in} \ (0,T] \ \times \ \Omega_{E}  \nonumber \\
\frac{\partial \rho_{o}}{\partial t} 
\, + \, 
\boldsymbol \nabla \cdot   \, \left( - \mu_{o} \, \boldsymbol \nabla \Phi \, \rho_{o} 
\, - \, D_{o} \,
\boldsymbol \nabla \, \rho_{o} \, \right) 
\; &= \; 0, && \text{in} \ (0,T] \ \times \ \Omega_{E}  \\
-\boldsymbol \nabla \cdot \left( \, \epsilon^{E}_{r} \, \boldsymbol \nabla \Phi \right)
\; &= \;
- \, \frac{q}{\epsilon_{0}} \, ( \rho_{r} - \rho_{o}). && \text{in} \ (0,T] \ \times \ \Omega_{E}   \nonumber
\end{align}
\f]

The reductant and hole mobilites are \f$\mu_{r}, \ \mu_{o} \  \f$ and the electrolyte permittivity is
\f$\epsilon_{r}^{E}\f$.  The lack of doping profile in our electrolyte reflects the fact 
that it is charge neutral.  Our model does not include any generation or recombination 
mechanisms in the electrolyte domain since we are only considering so-called 
``heterogenous reactions." That is chemical reactions can only occur at the interface 
and not within the bulk of the electrolyte.  

We assume the interface is isolated  and that the electrolyte variables take on their 
bulk values on boundary \f$\Gamma_{E}\f$,

\f[
\begin{align}
\rho_{r} \vert_{\Gamma_{E}}  \; = \;  \rho_{r}^{\infty}   ,
&&
\rho_{o} \vert_{\Gamma_{E}} \; = \; \rho_{o}^{\infty}  ,
&&
\Phi \vert_{\Gamma_{E}} \; = \; \Phi^{\infty} , 
&& \text{on} \ (0,T] \ \times \ \Gamma_{E} 
\end{align}
\f]


The potential and displacement electric field are continuous across the interface,

\f[
\begin{align}
\Phi \, \vert_{\Sigma^{-}} \, = \, \Phi \, \vert_{\Sigma^{+}} 
&&  -\epsilon_{r}^{S} \, \boldsymbol \nabla \, \Phi 
\ \cdot \ \boldsymbol \eta_{\Sigma} \ \vert_{\Sigma^{-}} 
\, = \,
-\epsilon^{E}_{r} \, \boldsymbol \nabla \, \Phi \,
\cdot \ \boldsymbol \eta_{\Sigma} \ \vert_{\Sigma^{+}} 
&& \text{on} \ (0,T] \ \times \ \Sigma
\end{align}
\f]

The chemical reactions of the charge carriers on the interface are modeled using 
the following boundary conditions on the currents,

\f[
\begin{align}
\textbf{J}_{n} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left(  \mu_{n}  \,  \,\boldsymbol \nabla \Phi \, \rho_{n} 
\, - \,  D_{n}  \, 
\boldsymbol \nabla \, \rho_{n}  \, \right)  \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
k_{et} \, ( \, \rho_{n} \, - \, \rho_{n}^{e} \, ) \, \rho_{o} , \\ 
\textbf{J}_{p} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left( - \mu_{p}  \, \boldsymbol \nabla \Phi \, \rho_{p} 
\, - \, D_{p}  \,
\boldsymbol \nabla \, \rho_{p} \, \right) \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
k_{ht} \, ( \, \rho_{p} \, - \, \rho_{p}^{e} \, ) \, \rho_{r},  \\ 
\textbf{J}_{r} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left(  \, \mu_{r}  \,\boldsymbol \nabla \Phi \, \rho_{r} 
\, - \, 
D_{r}  \,
\boldsymbol \nabla \, \rho_{r}  \, \right) \cdot \boldsymbol \eta_{\Sigma} \;  
&=
\;  \, k_{ht} \,  (\rho_{p} \, - \, \rho_{p}^{e}) \rho_{r} \, - \, 
k_{et} \, (\rho_{n} \, - \, \rho_{n}^{e}) \, \rho_{o}  , \\
\textbf{J}_{o} \cdot \boldsymbol \eta_{\Sigma}
\; = \; 
\left( - \mu_{o} \  \boldsymbol \nabla \Phi \, \rho_{o} 
\, - \, 
D_{o} \ 
\boldsymbol \nabla \, \rho_{o} \, \right) \cdot \boldsymbol \eta_{\Sigma} 
\; &= \; 
- k_{ht} \, (\rho_{p} \, - \, \rho_{p}^{e}) \, \rho_{r} \, + \, 
k_{et} \, (\rho_{n} \, - \, \rho_{n}^{e}) \, \, \rho_{o}.
\end{align}
\f]


The initial conditions are taken to be,
\f[ \begin{align}
\rho_{n}   \; &= \;  
\rho_{n}^{e}  && \text{on} \ \{0\} \ \times \ \Omega_{S}  \\
\rho_{p} \; &= \; 
\rho_{p}^{e}  && \text{on} \ \{0\} \ \times \ \Omega_{S} \\
\rho_{r}   \; &= \;  \rho_{r}^{\infty}  && \text{on} \ \{0\} \ \times \ \Omega_{E}  \\
\rho_{o} \; &= \; \rho_{o}^{\infty}  && \text{on} \ \{0\} \ \times \ \Omega_{E}
\end{align} 
\f]

\note We use Einstein's relations \f$ D \ = \ U_{T} \ \mu \f$ and singular perturbation scaling, for more information see the
<a href="http://www.ma.utexas.edu/users/gamba/papers/Semi-Elec-Interf.pdf">
paper</a> on this model.


The output of these simulations will be the calculations of potential, electric field, charge
densities and the current.  The although in our model we have eliminated the charge of an electron, \f$q\f$, the output current through the device will involve the charge of the electron and is defined to be:

\f[
\textbf{J}(\textbf{x}) \; = \;
\left\{
\begin{array}{cc}
q \, \textbf{J}_{n} \, - \, q \, \textbf{J}_{p}, & 
\textbf{x} \in \Omega_{S} \\
-q \, \textbf{J}_{r}  & 
\textbf{x} \in \Omega_{E} \\
\end{array}
\right.
\f]

\note This definition of the current is continuous across the interface and \f$\textbf{J}_{r} = -\textbf{J}_{o}\f$  at the interface and in the electrolyte.


\subsection Numerics

The overall stragy is to create a domain decomposition where in each subdomain we have
a reaction-drift-diffusion system of equations for a pair of charge carriers.  The
potential and electric field lives in a superset of these two domains and information is
passed between the two subdomains as well as superset domain.


The drift-diffusion transport equations are numerically approximated by 
the a local discontinuous Galerkin
method in space, see ddpCarrier_type  and the classes 
within this namesapace for more details.
Poisson's equation for the potential and electric field is approximated using a 
mixed finite element method see ddpPoisson_type for more details.


Time stepping is handeled in a specific way such that nonlinear terms are
linearized by time lagging.  Poisson is updated using implicit density values, 
while the charge carriers use an IMEX strategy. This is mostly handled by 
ddpTimeStepping_type and the specific IMEX strategy is written into each C++ driver
program. 


The overall time stepping strategy in the production driver 
is termed a ``parallel Gummel-Schwarz method."  The steps in this algorithm are
presented in the flow chart, \image html parallel_gummel_schwarz.png
For more deail see . 

All other device details are set in the input file: ddp-input.ini.



\subsection Resulting Output
The resulting output of this simulation will be a collection of .dat files. 
The number of files is controlled by the parameter: numTimeStamps in ddp-input.ini.
The format of the output files will be,

x point | electron density | hole density | elec field | potent. | current | G(x) | time | 
\n
.. \n
.. \n
interface  \n
.. \n
.. \n
x point | reductant density | oxidant density | elec field | potent. | current | G(x) | time | \n

\author Michael Harmon

