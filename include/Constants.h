#ifndef TBG_LMC_CONSTANTS_H
#define TBG_LMC_CONSTANTS_H

#include <complex>

// PARAMETERS ##########################################################

// intrinsic parameters
// electronic
#define THETA_RESCALED 1.0501208797943464								// Theta by rescaling
#define NATOM     11908                   							    // # atoms (dimension of Hamiltonian)
#define a0  	  2.46                                              	// lattice constant (Angstroem)                                        
#define	qq0       7.2													// hopping renormalization 
#define	a_NN      1.411621												// intralayer nearest-neigbour distance	
#define	c0        3.35                                             		// interlayer distance (Angstroem)
#define	t0        -2.7                                                  // hopping parameter of pz-pi (eV)
#define TT        1.3e-5 //1.3e-4                       				// inverse temperature (1/eV)
#define VV        0.0//0.00005                                          // symmetric top-gate/back-gate potential (eV) 
#define dgap      0.0//0.00050                                          / symmetric monolayer on-site potential
#define FILLING	  0	                                                    // filling factor (0 at CNP, -1 for v=-2 like in Aleski's paper)
#define UU		  -0.3	

// numerical paramters
#define mu_init   0.02											     	// initial guess for chemical potenital -> not necessarily needed for Floquet (eV)
#define dev_mu    1e-7                   					            // exit value for while loop in mu_SP() 
#define DELTA     1e-8												    // correction prefactor for chemical potential in groundstate() 
#define dev_order 1e-6                   					            // exit deviation for while loop in mu_SP() 
#define max_iter  30                                                    // max itereration for selfeconsitent calculation of order parameter
   
#define PI 3.14159265359

// truncation of Tylor matrices
#define dim_MAT	  3000	                                                // dimension of stored matrices
#define dim_new   512                                                   // target dimension
#define COUNT 0														    // dummy intex for multiple claculations		
#define eta 0.0001                                                      // imaginary width of frequency 

// Orders of Taylor propagation 
#define FIRST 1.0                                                       // 1st order (1.0 on / 0.0 off)
#define SECOND 1.0                                                      // 2nd order (1.0 on / 0.0 off)
#define THIRD 0.0                                                       // 3rd order (1.0 on / 0.0 off)

// Peierls driving
#define w_peierls      1.5                                              // Frequency of Applied Field (eV)
#define Ax_peierls     0.01                                             // Amplitude of Applied Field in x-direction
#define Ay_peierls     0.01                                             // Amplitude of Applied Field in y-direction
#define Az_peierls     0.0                                              // Amplitude of Applied Field in z-direction
 
// FLOQUET
#define m_max 1                                                         // order of truncation: m in {-m,...,0,...+m} 
#define n_max 1                                                         // order of truncation: n in {-n,...,0,...+n} (m_max == n_max!) 
#define timesteps_F 1e2                                                 // # of steps to perform integration over one period T=2pi/Omega                             

using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

constexpr cdouble II(0,1);

#endif //TBG_LMC_CONSTANTS_H
