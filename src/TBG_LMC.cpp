/** 
 *  C++ transfer of ALEKSI's CODE
 *	TIGHT-BINDING MODEL FOR TWISTED BILAYER GRAPHENE (TBG)
 *  Copyright (C) 2019, Gabriel E. Topp
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *  02111-1307, USA.
 */
 
// MAIN file 
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>
#include <time.h>


// PARAMETERS ##########################################################
#include "Constants.h"                          

// CALCULATION OPTIONS #################################################
#ifndef NO_MPI                                                                                                     
    #include <mpi.h>
#endif

#ifndef NO_OMP                                                          
    #include <omp.h>                                                    
#endif

//#define NO_IC                                                         
//#define NN_ONLY                                                       

#include "StandardFunctions.h"
#include "Functions.h"

//  main() function #####################################################
int main(int argc, char * argv[])
{
    //************** MPI INIT ***************************
  	int numprocs=1, myrank=0, namelen;
    //if(myrank==0){cout << "HALLO0"<< endl;}   
#ifndef NO_MPI
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Get_processor_name(processor_name, &namelen);
    
	//cout << "Process " << myrank << " on " << processor_name << " out of " << numprocs << " says hello." << endl;
	MPI_Barrier(MPI_COMM_WORLD);
    
#endif
	if(myrank==0) cout << "\n\tProgram running on " << numprocs << " processors." << endl;

	//************** OPEN_MP INIT **************************************
#ifndef NO_OMP 	  
	//cout << "# of processes " << omp_get_num_procs() << endl;
#pragma omp parallel 
	//cout << "Thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " says hello!" << endl;     
#endif
	//******************************************************************
    //if(myrank==0){cout << "HALLO1"<< endl;}   
	// DECLARATIONS AND INTITALIZATIONS
	const int m = 31;
	const int n = m+1;
	 
	const double THETA = acos(0.5*(double(m)*double(m) + double(n)*double(n) + 4.*double(m)*double(n))/(double(m)*double(m) + double(n)*double(n) + double(m)*double(n)));
	if(myrank==0){ 
		cout << "THETA: " << THETA*180./PI << endl;
		cout << "THETA_RESCALED: " << THETA_RESCALED << endl;
	}
	
	const double lambda = sin(THETA/2.)/(sin(THETA_RESCALED*(PI/180.)/2.));
	if(myrank==0) cout << "Lambda: " << lambda << endl;
	
	// Read in basis vectors
	vector<dvec> lvec;
	ReadIn(lvec, "L_VECS.dat");
	if(myrank==0){
		cout << "A0[0] = " <<  lvec[0][0] << " " << "A0[1] = " << lvec[0][1] << endl;
		cout << "A1[0] = " <<  lvec[1][0] << " " << "A1[1] = " << lvec[1][1] << endl;
	}
		
	// chemical potential
	double mu = 0.0;
	
	//Read in atomic positions
	vector<dvec> UNIT_CELL;
	ReadIn(UNIT_CELL, "Unit_Cell.dat");
	if(myrank==0) cout << "Unit_Cell.dat --> " <<  UNIT_CELL.size() << " points" << endl;
	if(NATOM != UNIT_CELL.size())
	{
		cout << "WRONG ATOMNUMBER!!! ---------------------------------------------------------------------------------------------" << endl;
		return 0;
	}
	
	//Read in vector of k-points
	vector<dvec> K_PATH;
	ReadIn(K_PATH, "k_path.dat");
	if(myrank==0) cout << "high-symmetry path --> " << K_PATH.size() << " points" << endl;
	int num_kpoints_PATH = K_PATH.size();
	
	//Read in vector of k-points (full path)
	vector<dvec> K_PATH_FULL;
	ReadIn(K_PATH_FULL, "k_path_full.dat");
	if(myrank==0) cout << "high-symmetry path (full) --> " << K_PATH_FULL.size() << " points" << endl;
	int num_kpoints_PATH_full = K_PATH_FULL.size();
	
	//~ // full BZ
	//~ //vector of weights
	//~ vector<dvec> kweights_full;
	//~ ReadIn(kweights_full, "k_weights_full.dat");
			
	//~ //vector of BZ vectors
	//~ vector<dvec> BZ_FULL;
	//~ ReadIn(BZ_FULL, "k_BZ_full.dat");
	//~ if(myrank==0) cout << "full BZ --> " << BZ_FULL.size() << " points" << endl;
	//~ int num_kpoints_BZ_full = BZ_FULL.size();
	
	//~ // reduced BZ
	//~ //vector of weights
	//~ vector<dvec> kweights_irr;
	//~ ReadIn(kweights_irr, "k_weights_irr.dat");
			
	//~ //vector of BZ vectors
	//~ vector<dvec> BZ_IRR;
	//~ ReadIn(BZ_IRR, "k_BZ_irr.dat");
	//~ if(myrank==0) cout << "irr. BZ --> " << BZ_IRR.size() << " points" << endl;
	//~ int num_kpoints_BZ_irr = BZ_IRR.size();
	
	// Single Particle (SP) MEMORY
	
	// SP vector for Hamiltonian Hk
	cvec *Hk = new cvec(NATOM*NATOM);
	
	// SP eigenvalues
	dvec evals(NATOM);
	
	// SP bands 
	dvec BANDS(num_kpoints_PATH_full*NATOM);
	
	// Order Parameter
	cvec ORDER(3*NATOM,0.);
	
	// Linear conductivity
	cvec COND_O1(2*4,0.);
	
	// Quadraticconductivity
	cvec COND_O2(27,0.);
	
	double DOS;
		
	// FLOQUET MEMORY
	// vector to store Floquet matrix
	cvec *Hk_FLOQUET_DOWN = new cvec((2*m_max+1)*(2*n_max+1)*dim_new*dim_new); 
    //cvec *Hk_FLOQUET = new cvec((2*m_max+1)*(2*n_max+1)*NATOM*NATOM);    
	
	// vector for eigenvalues
	dvec *evals_FLOQUET_DOWN = new dvec(dim_new*(2*n_max+1));
	//dvec *evals_FLOQUET = new dvec(NATOM*(2*n_max+1));
	
	// Floquet
	dvec *BANDS_FLOQUET_DOWN = new dvec(num_kpoints_PATH*dim_new*(2*n_max+1));
	dvec *OVERLAP_FLOQUET_DOWN = new dvec(num_kpoints_PATH*dim_new*(2*n_max+1));	
	
	//dvec *BANDS_FLOQUET = new dvec(num_kpoints_PATH*NATOM*(2*n_max+1));
	//dvec *OVERLAP_FLOQUET = new dvec(num_kpoints_PATH*NATOM*(2*n_max+1));	
			
	// Berry Curvature
	vector<dvec> bands_BCs_DOWN(8,dvec(dim_new));
	vector<dvec> bands_BCs_FLOQUET_DOWN(8,dvec(dim_new*(2*m_max+1)));	
	
	vector<dvec> bands_BCs(8,dvec(NATOM));
	vector<dvec> bands_BCs_FLOQUET(8,dvec(NATOM*(2*m_max+1)));	

	
	// Matrices for (truncated) Taylored matrices (new dimension dim_MAT!)
	vector<cvec*> Hk_Taylor(7); 
	vector<cvec*> Hk_DOWN_LIST_dim_new(7);                                        
	vector<cvec*> Hk_DOWN_LIST_dim_MAT(7);
	for(int m=0; m<7; m++)
	{
		Hk_Taylor[m] = new cvec(NATOM*NATOM);	
		Hk_DOWN_LIST_dim_new[m] = new cvec(dim_new*dim_new);   
		Hk_DOWN_LIST_dim_MAT[m] = new cvec(dim_MAT*dim_MAT);   
	}
	
	// CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	const clock_t begin_time = clock();                                 // time summed over all threads
#ifndef NO_OMP	 
	double dtime = omp_get_wtime();	                                    // time per core
#endif
	
	int count = COUNT;
	if(myrank==0){cout << "COUNT: " << count << endl;}
	
	//if(myrank==0){cout << "Start caluclation of chemical potential" << endl;}
	//mu_SP(lambda, Hk[0], evals, kweights_full, BZ_FULL, UNIT_CELL, lvec, mu, numprocs, myrank);

	//if(myrank==0){cout << "Start caluclation of equilibrium SP bands" << endl;}
	//if(count==0)
	//{
	//	Hk_bands_SP(lambda, BANDS, Hk[0], evals, K_PATH, UNIT_CELL, lvec, "DATA/bands_"+to_string(count)+".dat", numprocs, myrank);
	//}
	
	if(count==0)
	{
		if(myrank==0){cout << "Calculate truncated matrices for taylored Hamiltonian (PATH)" << endl;}
		Calc_List(lambda, evals, Hk[0], K_PATH, lvec, UNIT_CELL, Hk_Taylor, Hk_DOWN_LIST_dim_MAT, "PATH", numprocs, myrank);
	}
	
	//if(myrank==0){cout << "Calculate truncated matrices for  taylored Hamiltonian (BZ)" << endl;}
	//Calc_List(lambda, evals, Hk[0], BZ_FULL, lvec, UNIT_CELL, Hk_Taylor, Hk_DOWN_LIST_dim_MAT, "BZ", numprocs, myrank);
	
	if(myrank==0){cout << "Start caluclation of equilibrium bands" << endl;}
	Hk_bands_DOWN(K_PATH, Hk_DOWN_LIST_dim_new, "DATA/bands_DOWN_"+to_string(count)+".dat", numprocs, myrank);
	
	//if(myrank==0){cout << "Start caluclation of equilibrium Berry curvature along k-path" << endl;}	
	//EQ_BC_PATH(K_PATH, bands_BCs, Hk_DOWN_LIST_dim_new, "PATH", numprocs, myrank);
	
	//if(myrank==0){cout << "Start caluclation of equilibrium Berry curvature of BZ" << endl;}	
	//EQ_BC_PATH(BZ_FULL, bands_BCs, Hk_DOWN_LIST_dim_new, "BZ", numprocs, myrank);
	

	//if(myrank==0){cout << "Calculate Density of states" << endl;}
	//DOS_OMEGA(0.00001, -0.03, 0.03, mu, DOS, BZ_FULL, Hk_DOWN_LIST_dim_new, "DOS_"+to_string(count)+".dat", numprocs, myrank);
	
	//if(myrank==0){cout << "Calculate of linear conductivity" << endl;}
	//LIN_COND_OMEGA(0.00001, 0.0001, 0.010, mu, COND_O1, BZ_FULL, Hk_DOWN_LIST_dim_new, "SIGMA_"+to_string(count)+".dat", numprocs, myrank);
	//LIN_COND_OMEGA(0.00005, 5.0001, 5.010, mu, COND_O1, BZ_FULL, Hk_DOWN_LIST_dim_new, "HIGH_FREQU/SIGMA_"+to_string(count)+".dat", numprocs, myrank);
	
	//cdouble OMEGA1 = 0.005;
	//cdouble OMEGA2 = -0.005;
	//O2_COND(OMEGA1, OMEGA2, mu, COND_O1, BZ_FULL, Hk_DOWN_LIST_dim_new, "SIGMA_TEST.dat");
	
	//if(myrank==0){cout << "Calculate of quadratic conductivity" << endl;}
	//O2_COND_OMEGA(0.00001, 0.001, 0.020, -0.010, mu, COND_O2, BZ_FULL, Hk_DOWN_LIST_dim_new, "SIGMA_2_"+to_string(count)+".dat", numprocs, myrank);
	//O2_COND_OMEGA(0.00005, 5.0001, 5.010, -5.005, mu, COND_O2, BZ_FULL, Hk_DOWN_LIST_dim_new, "HIGH_FREQU/SIGMA_2_"+to_string(count)+".dat", numprocs, myrank);
	
	
	if(myrank==0){cout << "Start caluclation of effective Floquet bands along k-path (downfolded)" << endl;}		
	HF_EFF(Hk_DOWN_LIST_dim_new, K_PATH, "DATA/FLOQUET_EFF_"+to_string(count)+".dat", "PATH", numprocs, myrank);
	
	if(myrank==0){cout << "Start caluclation of Floqeut bands along k-path (truncated)" << endl;}		
	Hk_bands_Floquet_DOWN(BANDS_FLOQUET_DOWN[0], OVERLAP_FLOQUET_DOWN[0], Hk_FLOQUET_DOWN[0], evals_FLOQUET_DOWN[0], K_PATH, Hk_DOWN_LIST_dim_new, numprocs, myrank);

	//if(myrank==0){cout << "Start caluclation of Floquet Berry Curvature along kpath (truncated)" << endl;}	
	//FLOQUET_BC_PATH_DOWN(K_PATH, evals_FLOQUET_DOWN[0], bands_BCs_FLOQUET_DOWN, Hk_DOWN_LIST_dim_new, numprocs, myrank);
	
	//if(count==0)
	//{
	//	if(myrank==0){cout << "Start caluclation of Floqeut bands along k-path (full)" << endl;}		
	//	Hk_bands_Floquet(lambda, BANDS_FLOQUET[0], OVERLAP_FLOQUET[0], Hk_FLOQUET[0], evals_FLOQUET[0], K_PATH, UNIT_CELL, lvec, numprocs, myrank);
	//}

	//if(myrank==0){cout << "Start caluclation of Floquet Berry Curvature along kpath (full)" << endl;}	
	//FLOQUET_BC_PATH(lambda, K_PATH, evals_FLOQUET[0], bands_BCs_FLOQUET, UNIT_CELL, lvec, numprocs, myrank);


	if(myrank==0)
	{ 
		cout << "Total caluclation time (MPI): " << float(clock() - begin_time)/CLOCKS_PER_SEC << " seconds" << endl;
#ifndef NO_OMP	
		dtime = omp_get_wtime() - dtime;
		cout << "Total caluclation time (OMP): " << dtime << " seconds" << endl; 
#endif	
	}
	
#ifndef NO_MPI
	MPI_Finalize();
#endif	

// free memory	
	delete Hk;

	//delete Hk_FLOQUET;
	//delete evals_FLOQUET;
	//delete BANDS_FLOQUET;
	//delete OVERLAP_FLOQUET;
	
	delete Hk_FLOQUET_DOWN;
	delete evals_FLOQUET_DOWN;
	delete BANDS_FLOQUET_DOWN;
	delete OVERLAP_FLOQUET_DOWN;
	
	//delete Hk_BdG;
	//delete evals_BdG;
	//delete BANDS_BdG;

	
	for(int m=0; m<7; m++)
	{                            
		delete Hk_Taylor[m];
		delete Hk_DOWN_LIST_dim_new[m];
		delete Hk_DOWN_LIST_dim_MAT[m];
	}	
	
}



