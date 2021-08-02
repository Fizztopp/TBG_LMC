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
                                                    
#include "Constants.h"
#include "StandardFunctions.h"

#include <mpi.h>                                                         
#include <omp.h>   


void set_Hk0(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL)
/**
 * 	Set eq. single-particle Hamiltonian (without external field)
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 */
{
	const double t0_rs = t0/lambda;
	const double t1_rs = -0.11*t0;
	const double bb = a0*lambda/sqrt(3.);
	const double cc = 1.3618*a0*lambda;
	const double spatial_cut_off = 4.0*bb;
	
	const double kx = kvec[0];                                          
	const double ky = kvec[1];       
		
    // Bottom layer 
#ifndef NO_OMP    
    #pragma omp parallel												
	{
#endif		
	double d, rx, ry, rz;                                               
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 	
	for(int m=0; m<NATOM*NATOM; m++){
		Hk[m] = 0.0;
	}	
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 
	for(int i=0; i<NATOM/2; ++i)
	{
		Hk[fq(i,i,NATOM)] = VV/2.;                                     
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }     
		for(int j=i+1; j<NATOM/2; ++j)
		{	
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{   
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry));      				// [k] = 1/AA				          
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);
		}
	}	
	// Top layer 
#ifndef NO_OMP    	
	#pragma omp for
#endif	  
	for(int i=NATOM/2; i<NATOM; ++i)
	{
		Hk[fq(i,i,NATOM)] = -VV/2.;                                     
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }  
		for(int j=i+1; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif					
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry));      				// [k] = 1/AA					  
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
	// Inter-layer terms
#ifndef NO_IC
#ifndef NO_OMP    	
	#pragma omp for
#endif
	for(int i=0; i<NATOM/2; ++i)
	{
		for(int j=NATOM/2; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// == c0
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
					Hk[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry));      // [k] = 1/AA		
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
#endif
#ifndef NO_OMP 	
	}
#endif						
}


void set_Hk(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time)
/**
 * 	Set time-dependent single-particle Hamiltonian with external field
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */
{
	const double t0_rs = t0/lambda;
	const double t1_rs = -0.11*t0;
	const double bb = a0*lambda/sqrt(3.);
	const double cc = 1.3618*a0*lambda;
	const double spatial_cut_off = 4.0*bb;
	
	const double kx = kvec[0];                                          
	const double ky = kvec[1];       
		
    // Bottom layer 
#ifndef NO_OMP    
    #pragma omp parallel												
	{
#endif		
	double d, rx, ry, rz;                                               
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 	
	for(int m=0; m<NATOM*NATOM; m++){
		Hk[m] = 0.0;
	}	
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 
	for(int i=0; i<NATOM/2; ++i)
	{
		Hk[fq(i,i,NATOM)] = VV/2.;                                      
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }      
		for(int j=i+1; j<NATOM/2; ++j)
		{	
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{   
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz));      				// [k] = 1/AA				          
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);
		}
	}	
	// Top layer 
#ifndef NO_OMP    	
	#pragma omp for
#endif	  
	for(int i=NATOM/2; i<NATOM; ++i)
	{
		Hk[fq(i,i,NATOM)] = -VV/2.;                                     
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }  
		for(int j=i+1; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif					
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz));      				// [k] = 1/AA					  
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
	// Inter-layer terms
#ifndef NO_IC
#ifndef NO_OMP    	
	#pragma omp for
#endif
	for(int i=0; i<NATOM/2; ++i)
	{
		for(int j=NATOM/2; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// == c0
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
					Hk[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz));      // [k] = 1/AA		
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
#endif
#ifndef NO_OMP 	
	}
#endif						
}


void set_dHkdAx(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time)
/**
 * 	Set kx derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */
{
	const double t0_rs = t0/lambda;
	const double t1_rs = -0.11*t0;
	const double bb = a0*lambda/sqrt(3.);
	const double cc = 1.3618*a0*lambda;
	const double spatial_cut_off = 4.0*bb;
	
	const double kx = kvec[0];                                          
	const double ky = kvec[1];       
		
    // Bottom layer 
#ifndef NO_OMP    
    #pragma omp parallel												
	{
#endif		
	double d, rx, ry, rz;                                               
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 	
	for(int m=0; m<NATOM*NATOM; m++){
		Hk[m] = 0.0;
	}	
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 
	for(int i=0; i<NATOM/2; ++i)
	{    
		for(int j=i+1; j<NATOM/2; ++j)
		{	
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{   
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*rx);      				// [k] = 1/AA				          
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);
		}
	}	
	// Top layer 
#ifndef NO_OMP    	
	#pragma omp for
#endif	  
	for(int i=NATOM/2; i<NATOM; ++i)
	{
		Hk[fq(i,i,NATOM)] = -VV/2.;                                     
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }  
		for(int j=i+1; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif					
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*rx);      				// [k] = 1/AA					  
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
	// Inter-layer terms
#ifndef NO_IC
#ifndef NO_OMP    	
	#pragma omp for
#endif
	for(int i=0; i<NATOM/2; ++i)
	{
		for(int j=NATOM/2; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// == c0
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
					Hk[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*rx);      // [k] = 1/AA		
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
#endif
#ifndef NO_OMP 	
	}
#endif						
}


void set_dHkdAy(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time)
/**
 * 	Set ky derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */
{
	const double t0_rs = t0/lambda;
	const double t1_rs = -0.11*t0;
	const double bb = a0*lambda/sqrt(3.);
	const double cc = 1.3618*a0*lambda;
	const double spatial_cut_off = 4.0*bb;
	
	const double kx = kvec[0];                                         
	const double ky = kvec[1];       
		
    // Bottom layer 
#ifndef NO_OMP    
    #pragma omp parallel												
	{
#endif		
	double d, rx, ry, rz;                                               
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 	
	for(int m=0; m<NATOM*NATOM; m++){
		Hk[m] = 0.0;
	}	
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 
	for(int i=0; i<NATOM/2; ++i)
	{    
		for(int j=i+1; j<NATOM/2; ++j)
		{	
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{   
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*ry);      				// [k] = 1/AA				          
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);
		}
	}	
	// Top layer 
#ifndef NO_OMP    	
	#pragma omp for
#endif	  
	for(int i=NATOM/2; i<NATOM; ++i)
	{
		Hk[fq(i,i,NATOM)] = -VV/2.;                                     
		if(UNIT_CELL[i][3] < 0.9){
			Hk[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            Hk[fq(i,i,NATOM)] += dgap/2.;
        }  
		for(int j=i+1; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif					
					Hk[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*ry);      				// [k] = 1/AA					  
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
	// Inter-layer terms
#ifndef NO_IC
#ifndef NO_OMP    	
	#pragma omp for
#endif
	for(int i=0; i<NATOM/2; ++i)
	{
		for(int j=NATOM/2; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// == c0
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));    												// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
					Hk[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*exp(II*(Ax_t(time)*rx+Ay_t(time)*ry+Az_t(time)*rz))*(II*ry);      // [k] = 1/AA		
				}
			}
			Hk[fq(j,i,NATOM)] = conj(Hk[fq(i,j,NATOM)]);	
		}
	}
#endif
#ifndef NO_OMP 	
	}
#endif						
}


void set_Hk_Taylor(const double &lambda, dvec &kvec, vector<cvec*> Hk_Taylor, vector<dvec> &lvec, vector<dvec> &UNIT_CELL)
/**
 *	Taylor expansion of single-particle Hamiltonian for small fields
 * 	Set kx derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk_Taylor: Vector of complex matrices[10][NATOM*NATOM] to store Taylor matrices
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 */
{
	const double t0_rs = t0/lambda;
	const double t1_rs = -0.11*t0;
	const double bb = a0*lambda/sqrt(3.);
	const double cc = 1.3618*a0*lambda;
	const double spatial_cut_off = 4.0*bb;
	
	const double kx = kvec[0];                                          
	const double ky = kvec[1];       
	
    // bottom layer 
#ifndef NO_OMP    
    #pragma omp parallel												
	{
#endif		
	double d, rx, ry, rz;                                              
#ifndef NO_OMP    	
	#pragma omp for                                 
#endif		 	
	for(int m=0; m<NATOM*NATOM; m++){
		for(int n=0; n<7; n++)	{
			(*Hk_Taylor[n])[m] = 0.0;
		}	
	 }	
#ifndef NO_OMP    	
	#pragma omp for                                  					
#endif		 
	for(int i=0; i<NATOM/2; ++i)
	{
		//Back-gate voltage
		(*Hk_Taylor[0])[fq(i,i,NATOM)] = VV/2.;   
		if(UNIT_CELL[i][3] < 0.9){
			(*Hk_Taylor[0])[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            (*Hk_Taylor[0])[fq(i,i,NATOM)] += dgap/2.;
        }   
		for(int j=i+1; j<NATOM/2; ++j)
		{	
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));   
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif
					// 0th order
					(*Hk_Taylor[0])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry));     
					// 1st order
					(*Hk_Taylor[1])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*rx);     
					(*Hk_Taylor[2])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*ry);  
					(*Hk_Taylor[3])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*rz);  
					// 2nd order    	
					(*Hk_Taylor[4])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-rx*rx);      	
					(*Hk_Taylor[5])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-rx*ry);     	
					(*Hk_Taylor[6])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-ry*ry);    				
				}
			}
			for(int nn=0; nn<7; nn++)	{
				(*Hk_Taylor[nn])[fq(j,i,NATOM)]= conj((*Hk_Taylor[nn])[fq(i,j,NATOM)]);	
			}	
		}
	}	
	// Top layer 
#ifndef NO_OMP    	
	#pragma omp for
#endif	  
	for(int i=NATOM/2; i<NATOM; ++i)
	{
		(*Hk_Taylor[0])[fq(i,i,NATOM)] = -VV/2.;                                    
		if(UNIT_CELL[i][3] < 0.9){
			(*Hk_Taylor[0])[fq(i,i,NATOM)] += -dgap/2.;
        }
        else{
            (*Hk_Taylor[0])[fq(i,i,NATOM)] += dgap/2.;
        }  
		for(int j=i+1; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		// should be zero!
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));   											// [dd] = AA
	                if(d > spatial_cut_off){
						continue;
					}
#ifdef NN_ONLY  
					if(sqrt(rx*rx + ry*ry) >= 9./10.*lambda*a0){
						continue;
					}	
#endif					
					// 0th order
					(*Hk_Taylor[0])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry));     
					// 1st order
					(*Hk_Taylor[1])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*rx);     
					(*Hk_Taylor[2])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*ry);  
					(*Hk_Taylor[3])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(II*rz); 
					// 2nd order    	
					(*Hk_Taylor[4])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-rx*rx);      	
					(*Hk_Taylor[5])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-rx*ry);     	
					(*Hk_Taylor[6])[fq(i,j,NATOM)] += t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)/(d*d)*exp(II*(kx*rx+ky*ry))*(-ry*ry);    				
								           
				}
			}
			for(int nn=0; nn<7; nn++)	{
				(*Hk_Taylor[nn])[fq(j,i,NATOM)]= conj((*Hk_Taylor[nn])[fq(i,j,NATOM)]);	
			}	
		}
	}
	// Inter-layer terms 
#ifndef NO_IC
#ifndef NO_OMP    	
	#pragma omp for
#endif
	for(int i=0; i<NATOM/2; ++i)
	{
		for(int j=NATOM/2; j<NATOM; ++j)
		{
			for(int m=0; m<3; ++m)
			{
				for(int n=0; n<3; ++n)
				{
					rx = UNIT_CELL[i][0]-UNIT_CELL[j][0]+double(m-1)*lvec[0][0]+double(n-1)*lvec[1][0]; 		// [rx] = AA
					ry = double(m-1)*lvec[0][1]+UNIT_CELL[i][1]-UNIT_CELL[j][1]+double(n-1)*lvec[1][1]; 		// [ry] = AA
					rz = UNIT_CELL[i][2]-UNIT_CELL[j][2];	                                            		
	                d = sqrt(pow(rx,2.)+pow(ry,2.)+pow(rz,2.));   
	                if(d > spatial_cut_off){
						continue;
					}
						// 0th order   
						(*Hk_Taylor[0])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry)); 
						// 1st order    
						(*Hk_Taylor[1])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(II*rx);   
						(*Hk_Taylor[2])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(II*ry); 
						(*Hk_Taylor[3])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(II*rz); 
						// 2nd order	   
						(*Hk_Taylor[4])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(-rx*rx); 
						(*Hk_Taylor[5])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(-rx*ry); 
						(*Hk_Taylor[6])[fq(i,j,NATOM)] += 1./(d*d)*(t0_rs*exp(-qq0*abs(d-bb)/bb)*(rx*rx+ry*ry)+t1_rs*exp(-qq0*abs(d-cc)/bb)*rz*rz)*exp(II*(kx*rx+ky*ry))*(-ry*ry);          
	
				}
			}
			for(int nn=0; nn<7; nn++)	{
				(*Hk_Taylor[nn])[fq(j,i,NATOM)]= conj((*Hk_Taylor[nn])[fq(i,j,NATOM)]);	
			}
		}
	}
#endif
#ifndef NO_OMP 	
	}
#endif					
}


void set_Hk_DOWN_LIST(const double &lambda, dvec &kvec, dvec &evals, vector<cvec*> Hk_Taylor, vector<cvec*> Hk_DOWN_LIST, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, int &myrank)
/**
 * Transformes Taylor matrices to intital band basis & truncates matrices according to dimension dim_MAT. 
 * -lambda: renormalisation constant
 * -kvec: quasi momentum 
 * -evals: Real vector[NATOM] to store eigenvalues
 * -Hk_Taylor: Vector of complex matrices[10][NATOM*NATOM] to store Taylor matrices
 * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new*dim_new] to store truncated Taylor matrices in initial band basis
 * -lvec: superlattice bravis translational vectors (in lconst*Angstroem)
 * -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 * -myrank: Rank of process (MPI)
 */
 {
	if(myrank==0) cout << "Start set_Hk_DOWN_LIST() ---------------------------------------------------------------------------------------------------" << endl;
	int dimH = NATOM*NATOM;
	cvec *TEMP = new cvec(dimH); 
	cvec *SMAT = new cvec(dimH);
	
	double dtime1, dtime11;
	
	// Set band window with DIM_MAT
	vector<int> limits(2);
	limits[0] = NATOM/2-dim_MAT/2;
	limits[1] = NATOM/2+dim_MAT/2;
	
	if(myrank==0) cout << "set_Hk_DOWN_LIST(): Start diagonalization ----------------------------------------------------------------------------------" << endl; 
	set_Hk0(lambda, kvec, SMAT[0], lvec, UNIT_CELL);
	diagonalize(SMAT[0], evals);
	
	// Transform Tylor matrices to intital band basis
	for(int n=0; n<7; n++)	{
		dtime1 = omp_get_wtime();
		const clock_t begin_time1 = clock();
		times_nd(Hk_Taylor[n][0], SMAT[0], TEMP[0]);	
		dtime1 = omp_get_wtime() - dtime1;
		
		dtime11 = omp_get_wtime();
		const clock_t begin_time11 = clock();
		times(SMAT[0], TEMP[0], Hk_Taylor[n][0]);
		dtime11 = omp_get_wtime() - dtime11;
	}
	
	delete TEMP, SMAT;
	// Store truncated matrices
#ifndef NO_OMP 		
	#pragma omp parallel for
#endif 		
	for(int i=limits[0]; i<limits[1]; ++i)
	{
		for(int j=limits[0]; j<limits[1]; ++j)
		{
			for(int n=0; n<7; n++)	{
				(*Hk_DOWN_LIST[n])[fq(i-limits[0],j-limits[0],dim_MAT)] = (*Hk_Taylor[n])[fq(i,j,NATOM)]; 
			//if(myrank==0 && n==0)
			//{
			//	cout << (*Hk_Taylor[n])[fq(i,j,NATOM)] << endl;
			//}	 	
			}
		}		
	}	
}


void Calc_List(const double &lambda, dvec &evals, cvec &Hk, vector<dvec> &K_POINTS, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, vector<cvec*> Hk_Taylor, vector<cvec*> Hk_DOWN_LIST, const string& switcher, int &numprocs, int &myrank)
/**
 * Calculates truncated Taylor matrices along k-path in band basis and stores as file
 * -lambda: renormalisation constant
 * -evals: Real vector[NATOM] to store eigenvalues
 * -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 * -K_POINTS: vector[NPATH] of real vectors[3] to store k-points 
 * -lvec: superlattice bravis translational vectors (in lconst*Angstroem)
 * -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 * -Hk_Taylor: Vector of complex matrices[10][NATOM*NATOM] to store Taylor matrices
 * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new*dim_new] to store truncated Taylor matrices in initial band basis
 * -switcher: string to switch between k-path and BZ
 * -myrank: Rank of process (MPI)
 */
{
	if(myrank==0)cout << "Start Calc_List()" << " --------------------------------------------------------------------------------------------------------" << endl;
	int num_kpoints = K_POINTS.size();                     
	
	dvec kvec;	
	ofstream myfile;
	
	int count = COUNT;
	for(int k=myrank; k<num_kpoints; k+=numprocs)	
	{   
		// Set k-vector
		kvec = K_POINTS[k];
		// Sets Taylor expansion of Hk in k-orbital basis
		set_Hk_Taylor(lambda, kvec, Hk_Taylor, lvec, UNIT_CELL);
        // Declare Arrays needed for propagation and storage of observables                 		
		set_Hk_DOWN_LIST(lambda, kvec, evals, Hk_Taylor, Hk_DOWN_LIST, lvec, UNIT_CELL, myrank);
		// Write matrices to file (just one row)
		myfile.open("HK_TAYLOR/HK_TAYLOR_"+switcher+"_"+to_string(k)+".dat");
		if (myfile.is_open())
		{
			for(int m=0; m<7; m++)	
			{
				for(int i=0; i<dim_MAT*dim_MAT; ++i)
				{
					myfile << (*Hk_DOWN_LIST[m])[i] << " ";
				}
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
		if(myrank==0)
		{
			cout << "k_MAT #" << k << " stored!" << endl;
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}	
}


void ReadInMAT(vector<cvec*> Hk_DOWN_LIST, const string& filename)
{
/**
  * Read in taylore matrices from disc
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  *	-filename: String to define file
  */
	ifstream in(filename);
	string record;
	if(in.fail()){
		cout << "file" << filename << "could not be found!" << endl;
	}
	while (getline(in, record))
	{
		istringstream is( record );
		cvec row((istream_iterator<cdouble>(is)),	
		istream_iterator<cdouble>());
		//cout << row.size() << " " << dim_new << "---------------------------------------------------------------------" << endl;
		for(int m=0; m<7; ++m)	
		{	
			//~ for(int i=0; i<dim_new*dim_new; ++i)
			//~ {
				//~ (*Hk_DOWN_LIST[m])[i] = row[fq(m,i,dim_new*dim_new)];
				//~ //cout << row[fq(m,i,dim_new*dim_new)];
			//~ }
			for(int i=0; i<dim_new; ++i)
			{
				for(int j=0; j<dim_new; ++j)
				{
					(*Hk_DOWN_LIST[m])[fq(i,j,dim_new)] = row[fq(m,fq(i+dim_MAT/2-dim_new/2,j+dim_MAT/2-dim_new/2,dim_MAT),dim_MAT*dim_MAT)]; // cut out middle part!!
				}
			}
		}	

	}
	in.close();
}


void Hk_bands_DOWN(vector<dvec> &K_PATH, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank)
/**
 *	Calculate bands of downfolded Hamiltonian path K_PATH through BZ 
 *  -K_PATH: vector[NPATH] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int count = COUNT;
	const int num_kpoints_path = K_PATH.size();
	
	cvec Hk_DOWN(dim_new*dim_new);
	dvec evals_DOWN(dim_new);
	dvec BANDS(num_kpoints_path*dim_new);
	
	for(int k=myrank; k<num_kpoints_path; k+=numprocs)
	{
		// Read in tuncated taylor matrices from file
		ReadInMAT(Hk_DOWN_LIST, "HK_TAYLOR/HK_TAYLOR_PATH_"+to_string(k)+".dat");					
		for(int m=0; m<dim_new; m++)
			BANDS[fq(k, m, dim_new)] = real((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]);
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &BANDS[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void FUNC_DOS(cdouble &omega, double &mu, double &DOS, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename)
/**
 *	Calulates density-of-states for one frequency
 *  -omega: frequency
 *  -mu: chemical potential
 *  -DOS: denisty-of-states
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 */
{
	int count = COUNT;
	const int num_kpoints = BZ_FULL.size();                     // # of k-vectors from sampling of irreducible BZ
	DOS = 0.0;
	for(int k=0; k<num_kpoints; k++)
	{
		// Read in tuncated taylor matrices from file
		ReadInMAT(Hk_DOWN_LIST, "HK_TAYLOR/HK_TAYLOR_BZ_"+to_string(k)+".dat");
		for(int m=0; m<dim_new; m++)
		{
			DOS += -1./(double(num_kpoints)*PI)*imag(1./(omega+II*eta-(*Hk_DOWN_LIST[0])[fq(m,m,dim_new)])); 
		}			
	}
	if(filename!="no_file")
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{	
			myfile << DOS;
			myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void DOS_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, double &mu, double &DOS, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank)
/**
 *	Calulates density-of-states for frequency array
 *  -domega: frequency change
 *  -omega_min: minimum frequency
 *  -omega_max: maximum frequency
 *  -mu: chemical potential
 *  -DOS: denisty-of-states
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 *  -numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int omega_NN = int(real((omega_max-omega_min)/domega));
	cdouble omega;
	dvec DOS_MAT(omega_NN);
	for(int w=myrank; w<omega_NN; w+=numprocs)
	{
		omega = omega_min+mu + domega*double(w);                        // shift by mu to fit energy window with bands!
		FUNC_DOS(omega, mu, DOS, BZ_FULL, Hk_DOWN_LIST, "no_file");
		DOS_MAT[w] = DOS;
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &DOS_MAT[0], omega_NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int w=0; w<omega_NN; w++)
			{
				myfile << real(omega_min+mu + domega*double(w)) << " " << DOS_MAT[w] << endl;
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void LIN_COND(cdouble &omega, double &mu, cvec &COND_O1, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename)
/**
 *	Calulates linear conductivity for one frequency
 *  -omega: frequency
 *  -mu: chemical potential
 *  -COND_O1: linear conductivity
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 */
{
	const int num_kpoints = BZ_FULL.size();                     // # of k-vectors from sampling of irreducible BZ
	
	double f_nm;
	cdouble E_mn;

	int count = COUNT;
	for(int m=0; m<8; m++)
	{
		COND_O1[m] = 0.;
	}		
	for(int k=0; k<num_kpoints; k++)
	{
		// Read in tuncated taylor matrices from file
		ReadInMAT(Hk_DOWN_LIST, "HK_TAYLOR/HK_TAYLOR_BZ_"+to_string(k)+".dat");
		for(int lambda=0; lambda<2; lambda++)
		{
			for(int alpha=0; alpha<2; alpha++)
			{
				for(int m=dim_new/2-2; m<dim_new/2+2; m++) // ONLY FLAT BANDS!!
				{
					for(int n=0; n<dim_new; n++) // All bands
					{
						f_nm = fermi(real((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]), mu)-fermi(real((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]), mu);
						E_mn = (*Hk_DOWN_LIST[0])[fq(m,m,dim_new)] - (*Hk_DOWN_LIST[0])[fq(n,n,dim_new)];
						
						COND_O1[fq(lambda,alpha,2)] += II/(2.*double(num_kpoints))*1./(omega)*f_nm*(*Hk_DOWN_LIST[lambda+1])[fq(n,m,dim_new)]*(*Hk_DOWN_LIST[alpha+1])[fq(m,n,dim_new)]/(omega+II*eta-E_mn)*delta(m,n);        // diagonal part
						COND_O1[4+fq(lambda,alpha,2)] += II/(2.*double(num_kpoints))*1./(omega)*f_nm*(*Hk_DOWN_LIST[lambda+1])[fq(n,m,dim_new)]*(*Hk_DOWN_LIST[alpha+1])[fq(m,n,dim_new)]/(omega+II*eta-E_mn)*(1.-delta(m,n)); // off-diagonal part
					}		
				}	
				//COND_O1[fq(lambda,alpha,2)] += II/(2.*double(num_kpoints)*(omega+II*eta))*delta(lambda,alpha);
			}
		}			
	}	
	if(filename!="no_file")
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{	
			myfile << COND_O1[0] << COND_O1[1] << COND_O1[2] << COND_O1[3] << endl;
			myfile << COND_O1[4+0] << COND_O1[4+1] << COND_O1[4+2] << COND_O1[4+3];
			myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void LIN_COND_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, double &mu, cvec &COND_O1, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank)
/**
 *	Calulates density-of-states for frequency array
 *  -domega: frequency change
 *  -omega_min: minimum frequency
 *  -omega_max: maximum frequency
 *  -mu: chemical potential
 *  -COND_O1: linear conductivity
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 *  -numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int omega_NN = int(real((omega_max-omega_min)/domega));
	cvec COND_MAT(8*omega_NN);
	cdouble omega;
	for(int w=myrank; w<omega_NN; w+=numprocs)
	{
		omega = omega_min + domega*double(w);                           // only energy differences appear!
		LIN_COND(omega, mu, COND_O1, BZ_FULL, Hk_DOWN_LIST, "no_file");
		for(int i=0; i<8; i++)
			COND_MAT[fq(w,i,8)] = COND_O1[i];
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &COND_MAT[0], 8*omega_NN, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif		
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int w=0; w<omega_NN; w++)
			{
				myfile << real(omega_min + domega*double(w)) << " ";
				for(int m=0; m<8; m++)
				{
					myfile << real(COND_MAT[fq(w,m,8)]) << " ";
				}
				for(int m=0; m<8; m++)
				{
					myfile << imag(COND_MAT[fq(w,m,8)]) << " ";
				}
				myfile << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void O2_COND(cdouble &omega1, cdouble &omega2, double &mu, cvec &COND_O2, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename)
/**
 *	Calulates 2nd order conductivity for one frequency pair
 *  -omega1: frequency1
 *  -omega1: frequency2
 *  -mu: chemical potential
 *  -COND_O1: linear conductivity
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 */
{
	const int num_kpoints = BZ_FULL.size();                     // # of k-vectors from sampling of irreducible BZ
	
	double f_nl, f_lm;
	cdouble E_mn, E_ln, E_ml;

	int count = COUNT;
	for(int m=0; m<27; m++)
	{
		COND_O2[m] = 0.;
	}		
	for(int k=0; k<num_kpoints; k++)
	{
		// Read in tuncated taylor matrices from file
		ReadInMAT(Hk_DOWN_LIST, "HK_TAYLOR/HK_TAYLOR_BZ_"+to_string(k)+".dat");
		for(int lambda=0; lambda<3; lambda++)
		{
			for(int alpha=0; alpha<3; alpha++)
			{
				for(int beta=0; beta<3; beta++)
				{
					for(int m=0; m<dim_new; m++) // All bands
					{
						for(int n=0; n<dim_new; n++) // All bands
						{
							for(int l=0; l<dim_new; l++) // All bands
							{
								f_nl = fermi(real((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]), mu)-fermi(real((*Hk_DOWN_LIST[0])[fq(l,l,dim_new)]), mu);
								f_lm = fermi(real((*Hk_DOWN_LIST[0])[fq(l,l,dim_new)]), mu)-fermi(real((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]), mu);
								E_mn = (*Hk_DOWN_LIST[0])[fq(m,m,dim_new)] - (*Hk_DOWN_LIST[0])[fq(n,n,dim_new)];
								E_ln = (*Hk_DOWN_LIST[0])[fq(l,l,dim_new)] - (*Hk_DOWN_LIST[0])[fq(n,n,dim_new)];
								E_ml = (*Hk_DOWN_LIST[0])[fq(m,m,dim_new)] - (*Hk_DOWN_LIST[0])[fq(l,l,dim_new)];
								COND_O2[fq(lambda,alpha,3)+9*beta] += 1./(4.*double(num_kpoints)*(omega1)*(omega2))*(*Hk_DOWN_LIST[lambda+1])[fq(n,m,dim_new)]/(omega1+omega2+2.*II*eta-E_mn)*
																	 (f_nl*(*Hk_DOWN_LIST[alpha+1])[fq(l,n,dim_new)]*(*Hk_DOWN_LIST[beta+1])[fq(m,l,dim_new)]/(omega1+II*eta-E_ln) -
																	  f_lm*(*Hk_DOWN_LIST[alpha+1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[beta+1])[fq(l,n,dim_new)]/(omega1+II*eta-E_ml));        // diagonal part
							}
						}		
					}	
					//COND_O1[fq(lambda,alpha,2)] += II/(2.*double(num_kpoints)*(omega+II*eta))*delta(lambda,alpha);
				}	
			}
		}			
	}	
	if(filename!="no_file")
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{	
			for(int i=0; i<27; i++)
				myfile << COND_O2[i] << " ";
			myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void O2_COND_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, cdouble omega2, double &mu, cvec &COND_O2, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank)
/**
 *	Calulates 2nd order conductivity for frequency array
 *  -domega: frequency1 change
 *  -omega_min: minimum frequency1
 *  -omega_max: maximum frequency1
 *  -omega2: frequency2
 *  -mu: chemical potential
 *  -COND_O1: linear conductivity
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 *  -numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int omega_NN = int(real((omega_max-omega_min)/domega));
	cvec COND_MAT(27*omega_NN);
	cdouble omega;
	for(int w=myrank; w<omega_NN; w+=numprocs)
	{
		omega = omega_min + domega*double(w);                           // only energy differences appear!
		O2_COND(omega, omega2, mu, COND_O2, BZ_FULL, Hk_DOWN_LIST, "no_file");
		for(int i=0; i<27; i++)
			COND_MAT[fq(w,i,27)] = COND_O2[i];
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &COND_MAT[0], 27*omega_NN, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif		
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int w=0; w<omega_NN; w++)
			{
				myfile << real(omega_min + domega*double(w)) << " " << real(omega2) << " ";
				for(int m=0; m<27; m++)
				{
					myfile << real(COND_MAT[fq(w,m,27)]) << " ";
				}
				for(int m=0; m<27; m++)
				{
					myfile << imag(COND_MAT[fq(w,m,27)]) << " ";
				}
				myfile << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void EQ_BC(int k, vector<dvec> &OFF_BC, vector<dvec> &bands_BCs, vector<cvec*> Hk_DOWN_LIST, const string& filename, const string& switcher)
/**
 * 	Calculate quantum geometric tensor (QGT) in equlibrium for single k-point by derivatives of Hamiltonian
 *  -k: integer to pick k-vector from array
 *  -OFF_BC: matric to store interband terms
 *  -bands_BCs: matrix to store intraband terms
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -filename: string to name file
 *  -switcher: string to switch between k-path and BZ
 */
{
	int count = COUNT;
	ReadInMAT(Hk_DOWN_LIST,"HK_TAYLOR/HK_TAYLOR_"+switcher+"_"+to_string(k)+".dat");
	// Calculate Phase around loop
	for(int m=0; m<dim_new; m++)	
	{
		for(int j=0; j<8; j++){
			bands_BCs[j][m] = 0.;
		}	
		for(int n=0; n<dim_new; n++)
		{
			if(m==n)
				continue;
			// Metric (intraband)
			bands_BCs[0][m] += real((*Hk_DOWN_LIST[1])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[1][m] += real((*Hk_DOWN_LIST[1])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[2][m] += real((*Hk_DOWN_LIST[2])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[3][m] += real((*Hk_DOWN_LIST[2])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			
			// Berry curvature (intraband)
			bands_BCs[4][m] += -2.*imag((*Hk_DOWN_LIST[1])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[5][m] += -2.*imag((*Hk_DOWN_LIST[1])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[6][m] += -2.*imag((*Hk_DOWN_LIST[2])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			bands_BCs[7][m] += -2.*imag((*Hk_DOWN_LIST[2])[fq(m,n,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,m,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])*((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(n,n,dim_new)])));
			
			// Berry curvature (interband)
			for(int j=0; j<8; j++){
				OFF_BC[j][fq(m,n,dim_new)] = 0.;
			}	
			for(int l=0; l<dim_new; l++)
			{
				if(l==n)
					continue;
				if(l==m)
					continue;
				// Metric (interband)
				OFF_BC[0][fq(m,n,dim_new)] += real((*Hk_DOWN_LIST[1])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
				OFF_BC[1][fq(m,n,dim_new)] += real((*Hk_DOWN_LIST[1])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));	
				OFF_BC[2][fq(m,n,dim_new)] += real((*Hk_DOWN_LIST[2])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
				OFF_BC[3][fq(m,n,dim_new)] += real((*Hk_DOWN_LIST[2])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
			
				// BC (interband)
				OFF_BC[4][fq(m,n,dim_new)] += -2.*imag((*Hk_DOWN_LIST[1])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
				OFF_BC[5][fq(m,n,dim_new)] += -2.*imag((*Hk_DOWN_LIST[1])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));	
				OFF_BC[6][fq(m,n,dim_new)] += -2.*imag((*Hk_DOWN_LIST[2])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[1])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
				OFF_BC[7][fq(m,n,dim_new)] += -2.*imag((*Hk_DOWN_LIST[2])[fq(l,m,dim_new)]*(*Hk_DOWN_LIST[2])[fq(n,l,dim_new)]/(((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*((*Hk_DOWN_LIST[0])[fq(n,n,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])));
			}
		}			
	}	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0)
	{	
		if(filename!="no_file")
		{	
			ofstream myfile (filename);
			if (myfile.is_open())
			{
				for(int n=0; n<dim_new; ++n) 
				{
					// Berry curvature equal to phase diveded by area of loop
					for(int j=0; j<8; j++){
						myfile << bands_BCs[j][n] << " ";
					}	
					myfile << endl;	
				}	
				myfile.close();
			}
			else cout << "Unable to open file" << endl;	
		}
	}		
}	


void EQ_BC_PATH(vector<dvec> &K_POINTS, vector<dvec> &bands_BCs, vector<cvec*> Hk_DOWN_LIST, const string& switcher, int &numprocs, int &myrank)
/**
 * 	Calculate quantum geometric tensor (QGT) in equlibrium for momuntum array
 *  -K_POINTS: vector[NBZ] of real vectors[3] to store k-points 
 *  -bands_BCs: matrix to store intraband terms
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -switcher: string to switch between k-path and BZ
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int num_kpoints = K_POINTS.size();
	dvec GE_ARRAY(4*num_kpoints*dim_new);  
	dvec BC_ARRAY(4*num_kpoints*dim_new);   
	
	vector<dvec> OFF_BC(8, dvec(dim_new*dim_new));      
	dvec OFF_GE_ARRAY(4*num_kpoints*dim_new*dim_new);
	dvec OFF_BC_ARRAY(4*num_kpoints*dim_new*dim_new);
                                       
	
	for(int k=myrank; k<num_kpoints; k+=numprocs)
	{
		EQ_BC(k, OFF_BC, bands_BCs, Hk_DOWN_LIST, "no_file", switcher);
		for(int n=0; n<dim_new; ++n) 
		{
			for(int j=0; j<4; j++){
				GE_ARRAY[fq(k,n+j*dim_new,4*dim_new)] = bands_BCs[j][n];
				BC_ARRAY[fq(k,n+j*dim_new,4*dim_new)] = bands_BCs[j+4][n];			
			}	
			for(int m=0; m<dim_new; ++m) 
			{
				for(int j=0; j<4; j++){
					OFF_GE_ARRAY[fq(k,fq(m,n,dim_new)+j*dim_new*dim_new,4*dim_new*dim_new)] = OFF_BC[j][fq(m,n,dim_new)];
					OFF_BC_ARRAY[fq(k,fq(m,n,dim_new)+j*dim_new*dim_new,4*dim_new*dim_new)] = OFF_BC[j+4][fq(m,n,dim_new)];			
				}	
			}	
		}	
	}	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &GE_ARRAY[0], 4*dim_new*num_kpoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &BC_ARRAY[0], 4*dim_new*num_kpoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &OFF_GE_ARRAY[0], 4*num_kpoints*dim_new*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &OFF_BC_ARRAY[0], 4*num_kpoints*dim_new*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	int count = COUNT;
	if(myrank==0)
	{	
		ofstream myfile ("DATA/EQ_GE_LOOP_"+switcher+"_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*dim_new; ++n) 
				{
					myfile << GE_ARRAY[fq(k,n,4*dim_new)] << " ";	
				}	
				myfile << endl;
			}
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile ("DATA/EQ_BC_LOOP_"+switcher+"_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*dim_new; ++n) 
				{
					myfile << BC_ARRAY[fq(k,n,4*dim_new)] << " ";
				}	
				myfile << endl;
			}
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile ("DATA/EQ_GE_OFF_"+switcher+"_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*dim_new*dim_new; ++n) 
				{
					myfile << OFF_BC_ARRAY[fq(k,n,4*dim_new*dim_new)] << " ";
				}	
				myfile << endl;
			}
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile ("DATA/EQ_GE_OFF_"+switcher+"_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*dim_new*dim_new; ++n) 
				{
					myfile << OFF_GE_ARRAY[fq(k,n,4*dim_new*dim_new)] << " ";
				}	
				myfile << endl;
			}
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}	
}	


void HF_EFF(vector<cvec*> Hk_DOWN_LIST, vector<dvec> &K_PATH, const string& filename, const string& switcher, int &numprocs, int &myrank)
/**
 * 	Calculate the effective Floquethamiltonian in the high frequency expansion
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -K_PATH: vector[NBZ] of real vectors[3] to store k-points 
 *  -filename: string to name file
 *  -switcher: string to switch between k-path and BZ
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	int count = COUNT;
	const int num_kpoints_path = K_PATH.size();
	
	cvec TEMP(dim_new*dim_new);
	cvec Hk_DOWN(dim_new*dim_new);
	cvec Hk_DOWN0(dim_new*dim_new);
	cvec Hk_DOWN1(dim_new*dim_new);
	cvec Hk_DOWN2(dim_new*dim_new);
	cvec Hk_DOWN3(dim_new*dim_new);
	cvec Hk_DOWN4(dim_new*dim_new);
	cvec Hk_DOWN5(dim_new*dim_new);
	cvec Hk_DOWN6(dim_new*dim_new);
	cvec Hk_DOWN7(dim_new*dim_new);
	cvec Hk_DOWN8(dim_new*dim_new);
	cvec Hk_DOWN9(dim_new*dim_new);
	
	dvec evals_DOWN(dim_new);
	dvec BANDS(num_kpoints_path*dim_new);
	dvec BANDS0(num_kpoints_path*dim_new);
	dvec BANDS1(num_kpoints_path*dim_new);
	dvec BANDS2(num_kpoints_path*dim_new);
	dvec BANDS3(num_kpoints_path*dim_new);
	dvec BANDS4(num_kpoints_path*dim_new);
	dvec BANDS5(num_kpoints_path*dim_new);
	dvec BANDS6(num_kpoints_path*dim_new);
	dvec BANDS7(num_kpoints_path*dim_new);
	dvec BANDS8(num_kpoints_path*dim_new);
	dvec BANDS9(num_kpoints_path*dim_new);
	
		
	for(int k=myrank; k<num_kpoints_path; k+=numprocs)
	{
		ReadInMAT(Hk_DOWN_LIST,"HK_TAYLOR/HK_TAYLOR_"+switcher+"_"+to_string(k)+".dat");
		// Calculate Phase around loop
#ifndef NO_OMP 		
	#pragma omp parallel for
#endif 			
		for(int m=0; m<dim_new; m++)	
		{
			for(int n=0; n<dim_new; n++)
			{	
				Hk_DOWN[fq(m,n,dim_new)] = real((*Hk_DOWN_LIST[0])[fq(m,n,dim_new)])*delta(m,n);
				Hk_DOWN0[fq(m,n,dim_new)] = real((*Hk_DOWN_LIST[0])[fq(m,n,dim_new)])*delta(m,n);
				Hk_DOWN1[fq(m,n,dim_new)] = 0.;
				Hk_DOWN2[fq(m,n,dim_new)] = 0.;
				Hk_DOWN3[fq(m,n,dim_new)] = 0.;
				Hk_DOWN4[fq(m,n,dim_new)] = 0.;
				Hk_DOWN5[fq(m,n,dim_new)] = 0.;
				Hk_DOWN6[fq(m,n,dim_new)] = 0.;
				Hk_DOWN7[fq(m,n,dim_new)] = 0.;
				Hk_DOWN8[fq(m,n,dim_new)] = 0.;
				Hk_DOWN9[fq(m,n,dim_new)] = 0.;
				
				// Quadratic coupling H_AA
				Hk_DOWN[fq(m,n,dim_new)]  += 0.25*Ax_peierls*Ay_peierls*((*Hk_DOWN_LIST[4])[fq(m,n,dim_new)]+(*Hk_DOWN_LIST[6])[fq(m,n,dim_new)]);
				// H_AA^00 diag
				Hk_DOWN3[fq(m,n,dim_new)]  += 0.25*Ax_peierls*Ay_peierls*((*Hk_DOWN_LIST[4])[fq(m,n,dim_new)]+(*Hk_DOWN_LIST[6])[fq(m,n,dim_new)])*delta(m,n);	
				// H_AA^00 off-diag
				Hk_DOWN4[fq(m,n,dim_new)]  += 0.25*Ax_peierls*Ay_peierls*((*Hk_DOWN_LIST[4])[fq(m,n,dim_new)]+(*Hk_DOWN_LIST[6])[fq(m,n,dim_new)])*(1.-delta(m,n));		
				
				for(int l=0; l<dim_new; l++)
				{	
					// Linear coupling H_A
					Hk_DOWN[fq(m,n,dim_new)]  += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)]); 
					// H_A^0+/-1 diag 
					Hk_DOWN1[fq(m,n,dim_new)] += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)])*delta(m,n); 
					// H_A^0+/-1 off-diag
					Hk_DOWN2[fq(m,n,dim_new)] += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)])*(1.-delta(m,n)); 
					
					// Quadratic coupling H_AA
					Hk_DOWN[fq(m,n,dim_new)]  += 1./16.*Ax_peierls*Ax_peierls*Ay_peierls*Ay_peierls/w_peierls*II*
											     ((*Hk_DOWN_LIST[4])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[4])[fq(l,n,dim_new)] 
											     +(*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[6])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[6])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)]); 						
					// H_AA^0+/-2 diag
					Hk_DOWN5[fq(m,n,dim_new)]  += 1./16.*Ax_peierls*Ax_peierls*Ay_peierls*Ay_peierls/w_peierls*II*
											     ((*Hk_DOWN_LIST[4])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[4])[fq(l,n,dim_new)] 
											     +(*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[6])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[6])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)])*delta(m,n); 				
					// H_AA^0+/-2 off-diag
					Hk_DOWN6[fq(m,n,dim_new)]  += 1./16.*Ax_peierls*Ax_peierls*Ay_peierls*Ay_peierls/w_peierls*II*
											     ((*Hk_DOWN_LIST[4])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[4])[fq(l,n,dim_new)] 
											     +(*Hk_DOWN_LIST[5])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[6])[fq(l,n,dim_new)] - (*Hk_DOWN_LIST[6])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[5])[fq(l,n,dim_new)])*(1.-delta(m,n)); 			     
					
					// H_AA^00 diag nontrivial
					if(l==m)
						continue;
					if(abs((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])<1e-6)
						continue;
					Hk_DOWN9[fq(m,n,dim_new)]  += -0.5*Ax_peierls*Ay_peierls*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,m,dim_new)]+(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,m,dim_new)])/((*Hk_DOWN_LIST[0])[fq(m,m,dim_new)]-(*Hk_DOWN_LIST[0])[fq(l,l,dim_new)])*delta(m,n);	
											     
				}
			}		
		}
		// Calculate linear off diagonal part which comes only from flat-flat coupling
		for(int m=dim_new/2-2; m<dim_new/2+2; m++)	
		{
			for(int n=dim_new/2-2; n<dim_new/2+2; n++)
			{	
				for(int l=0; l<dim_new; l++)
				{		
					Hk_DOWN7[fq(m,n,dim_new)] += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)])*(1.-delta(m,n)); 
				}
			}
		}		
		// Calculate linear off diagonal part which comes only from flat-dispersive coupling
		for(int m=dim_new/2-2; m<dim_new/2+2; m++)	
		{
			for(int n=0; n<dim_new/2-2; n++)
			{	
				for(int l=0; l<dim_new; l++)
				{		
					Hk_DOWN8[fq(m,n,dim_new)] += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)])*(1.-delta(m,n)); 
				}
				Hk_DOWN8[fq(n,m,dim_new)] += conj(Hk_DOWN8[fq(m,n,dim_new)]); 
			}
			for(int n=dim_new/2+2; n<dim_new; n++)
			{	
				for(int l=0; l<dim_new; l++)
				{		
					Hk_DOWN8[fq(m,n,dim_new)] += 0.5*Ax_peierls*Ay_peierls/w_peierls*II*((*Hk_DOWN_LIST[1])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[2])[fq(l,n,dim_new)]-(*Hk_DOWN_LIST[2])[fq(m,l,dim_new)]*(*Hk_DOWN_LIST[1])[fq(l,n,dim_new)])*(1.-delta(m,n)); 
				}
				Hk_DOWN8[fq(n,m,dim_new)] += conj(Hk_DOWN8[fq(m,n,dim_new)]);
			}
		}	
					
		diagonalize_DOWN(Hk_DOWN, evals_DOWN);	
		
		//Transform matrices into Floquet band basis
		times_nd(Hk_DOWN0, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN0);
		
		times_nd(Hk_DOWN1, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN1);
		
		times_nd(Hk_DOWN2, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN2);
		
		times_nd(Hk_DOWN3, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN3);
		
		times_nd(Hk_DOWN4, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN4);
		
		times_nd(Hk_DOWN5, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN5);
		
		times_nd(Hk_DOWN6, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN6);
		
		times_nd(Hk_DOWN7, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN7);
		
		times_nd(Hk_DOWN8, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN8);
		
		times_nd(Hk_DOWN9, Hk_DOWN, TEMP);	
		times(Hk_DOWN, TEMP, Hk_DOWN9);
	
		for(int i=0; i<dim_new;i++)
		{
			BANDS[fq(k, i, dim_new)] = evals_DOWN[i];
			BANDS0[fq(k, i, dim_new)] = real(Hk_DOWN0[fq(i,i,dim_new)]);
			BANDS1[fq(k, i, dim_new)] = real(Hk_DOWN1[fq(i,i,dim_new)]);
			BANDS2[fq(k, i, dim_new)] = real(Hk_DOWN2[fq(i,i,dim_new)]);
			BANDS3[fq(k, i, dim_new)] = real(Hk_DOWN3[fq(i,i,dim_new)]);
			BANDS4[fq(k, i, dim_new)] = real(Hk_DOWN4[fq(i,i,dim_new)]);
			BANDS5[fq(k, i, dim_new)] = real(Hk_DOWN5[fq(i,i,dim_new)]);
			BANDS6[fq(k, i, dim_new)] = real(Hk_DOWN6[fq(i,i,dim_new)]);
			BANDS7[fq(k, i, dim_new)] = real(Hk_DOWN7[fq(i,i,dim_new)]);
			BANDS8[fq(k, i, dim_new)] = real(Hk_DOWN8[fq(i,i,dim_new)]);
			BANDS9[fq(k, i, dim_new)] = real(Hk_DOWN9[fq(i,i,dim_new)]);
		}	
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &BANDS[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS0[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS1[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS2[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS3[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS4[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS5[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS6[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);		
	MPI_Allreduce(MPI_IN_PLACE, &BANDS7[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &BANDS8[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
	MPI_Allreduce(MPI_IN_PLACE, &BANDS9[0], num_kpoints_path*dim_new, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	
#endif	
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
			if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF0_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS0[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
			if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF1_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS1[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
			if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF2_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS2[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
			if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF3_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS3[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
				if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF4_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS4[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
	if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF5_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS5[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
	if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF6_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS6[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
	if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF7_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS7[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}		
	if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF8_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS8[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}		
	if(myrank==0)
	{
		ofstream myfile ("DATA/FLOQUET_EFF9_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<dim_new; m++)
				{
					myfile << BANDS9[fq(k, m, dim_new)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}			
}	


void Hk_bands_SP(const double &lambda, dvec &BANDS, cvec &Hk, dvec &evals, vector<dvec> &K_PATH, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, const string& filename, int &numprocs, int &myrank)
/**
 *	Calculate bands of Hk0(k) for path K_PATH through BZ 
 *  -lambda: renormalisation constant
 *  -BANDS: Vector to store eigenvalues of all k-points 
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -evals: Vector to store eigenvalues of diagonalization
 *  -K_PATH: Vector of high-symmetry path vectors
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real Vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 * 	-filename: String to store data
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{
	const int num_kpoints_path = K_PATH.size();

	for(int k=myrank; k<num_kpoints_path; k+=numprocs)
	{
		if(myrank==0) cout << k << endl;
		set_Hk0(lambda, K_PATH[k], Hk, lvec, UNIT_CELL);
		diagonalize_eig(Hk, evals);
		for(int m=0; m<NATOM; m++)
			BANDS[fq(k, m, NATOM)] = evals[m];
	}
#ifndef NO_MPI	
	MPI_Allreduce(MPI_IN_PLACE, &BANDS[0], num_kpoints_path*NATOM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	if(myrank==0)
	{
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<NATOM; m++)
				{
					myfile << BANDS[fq(k, m, NATOM)] << " " ;
				}
				myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}	
}


void mu_SP(const double &lambda, cvec &Hk, dvec &evals, vector<dvec> &kweights, vector<dvec> &BZ, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, double &mu, int &numprocs, int &myrank)
/**
 * 	Calculate initial chemical potential
 * -lambda: renormalisation constant
 *  -Hk: Complex vector[NATOM*NATOM] to store Hamiltonian
 *  -evals: Real vector[NATOM] to store eigenvalues
 *  -kweights: Real vector containing weights of k-points
 *  -BZ: k-points of irreducable reciprocal cell
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 *  -mu: Chemical potential
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */
{	
	const int num_kpoints_BZ = BZ.size();  
	
	int count = 0;                                                      
	double N_tot;	
	double mu_old;    
	double deviation = 1.0;
	mu = mu_init;
	double num_kpoints_BZ_full = 0.;
	dvec EVALS(num_kpoints_BZ*NATOM);                                   
	
	for(int k=0; k<kweights.size(); k++)
		num_kpoints_BZ_full += kweights[k][0];							
	
	for(int k=myrank; k<num_kpoints_BZ; k+=numprocs)
	{		
		set_Hk0(lambda, BZ[k], Hk, lvec, UNIT_CELL);
		diagonalize_eig(Hk, evals);                                 	
		for(int i=0; i<NATOM; i++)
		{			
			EVALS[fq(k,i,NATOM)] = evals[i];	
		}
	}
	
	while(deviation > dev_mu)
	{
		count++;
					
		mu_old = mu;	
	    
		N_tot = 0.;
		for(int k=myrank; k<num_kpoints_BZ; k+=numprocs)
		{		
			for(int i=0; i<NATOM; i++)
			{			
				N_tot +=  fermi(EVALS[fq(k,i,NATOM)], mu)*kweights[k][0];	
			}
		}

#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &N_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
		mu += -DELTA*(N_tot-(0.5*double(NATOM)+double(FILLING))*num_kpoints_BZ_full);	                       
		
		deviation = abs(mu-mu_old);
		if(myrank==0){
			cout << "loop #" << count << ": deviation = " << deviation << endl;
			cout << "N_tot " << N_tot/num_kpoints_BZ_full << endl;
			cout << "chemical potential mu = " << mu << endl;
		}	
	}
	if(myrank==0)
	{
		cout << "Particle number per unit cell: " << N_tot/num_kpoints_BZ_full << endl;
		ofstream myfile ("DATA/mu.dat");
		if (myfile.is_open())
		{
			myfile << mu;
			myfile.close();
		}	
	else cout << "Unable to open file" << endl;	
	}
}


void set_Hk_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time)
/**
  *	Set downfolded td Hamiltonian 
  * -dim_new: integer value of reduced leading order of Hamiltonian
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -ASD: Gauge field of source-drain field
  * -Pol: double to set chirality
  * -time: tiome variable
  */
{
	double AX = Ax_t(time);
	double AY = Ay_t(time);
	
	for(int i=0; i<dim_new*dim_new; ++i){
		Hk_DOWN[i] = (*Hk_DOWN_LIST[0])[i] + FIRST*((*Hk_DOWN_LIST[1])[i]*AX + (*Hk_DOWN_LIST[2])[i]*AY) + SECOND*1./2.*((*Hk_DOWN_LIST[4])[i]*AX*AX + 2.*(*Hk_DOWN_LIST[5])[i]*AX*AY + (*Hk_DOWN_LIST[6])[i]*AY*AY);	
	}	
}	


void set_dHkdAx_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time) 
/**
  * Set downfolded t.-d. derivative by Ax of Hamiltonian (original band basis)
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -time: tiome variable
  */
{
	double AX = Ax_t(time);
	double AY = Ay_t(time);
	
	for(int i=0; i<dim_new*dim_new; ++i){
		Hk_DOWN[i] =  (*Hk_DOWN_LIST[1])[i] + FIRST*(*Hk_DOWN_LIST[4])[i]*AX + (*Hk_DOWN_LIST[5])[i]*AY; 	
	}		
}	


void set_dHkdAy_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time) 
/**
  * Set downfolded t.-d. derivative by Ay of Hamiltonian (original band basis)
  * -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -time: tiome variable
  */
{
	double AX = Ax_t(time);
	double AY = Ay_t(time);

	for(int i=0; i<dim_new*dim_new; ++i){
		Hk_DOWN[i] = (*Hk_DOWN_LIST[2])[i] + FIRST*(*Hk_DOWN_LIST[5])[i]*AX + (*Hk_DOWN_LIST[6])[i]*AY;  		
	}		
}	


void Hk_bands_Floquet_DOWN(dvec &BANDS_FLOQUET, dvec &OVERLAP_FLOQUET, cvec &Hk_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &K_PATH, vector<cvec*> Hk_DOWN_LIST, int &numprocs, int &myrank)
/**
 *	Calculate Floquet bands by truncated expansion in Floquet eigenfunctions
 *  -BANDS_FLOQUET: Real vector to store Floquet eigenvalues of all k-points 
 *  -OVERLAP_FLOQUET: Real vector[num_kpoints_PATHxNATOMx(2*n_max+1)] to store overlap ov Flquet bands with equilibrium bands
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store Flqoeut Hamiltonian matrix
 *  -evals_FLOQUET: Real vector[(M_max+1) x NATOM] to store Floquet eigenvalues
 *  -K_PATH: vector of high-symmetry path vectors
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
{
	const int count = COUNT;
	const int DIM_FLOQ = (2*n_max+1)*dim_new;
	
	const int num_kpoints_path = K_PATH.size();
	const double T = 2.*M_PI/w_peierls;
	const double dt = T/double(timesteps_F-1);
	
	cvec *TEMP1 = new cvec(dim_new*dim_new);
	cvec *TEMP2 = new cvec(dim_new*dim_new); 
	double temp; 
	
	for(int k=myrank; k<num_kpoints_path; k+=numprocs)
	{
#ifndef NO_OMP    	
	#pragma omp parallel for collapse(4)                                
#endif						                
		for(int m=-m_max; m<m_max+1; m++)
		{
			for(int n=-n_max; n<n_max+1; n++)
			{					
				for(int i=0; i<dim_new; i++)
				{
					for(int j=0; j<dim_new; j++)
					{
						Hk_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] = 0.0;
					}
				}
			}
		}												
		if(myrank==0) cout << endl; 
		if(myrank==0) cout << "k = " << k << endl;	
		
		ReadInMAT(Hk_DOWN_LIST,"HK_TAYLOR/HK_TAYLOR_PATH_"+to_string(k)+".dat");

		// Perform integration over one period T
		for(double t=0; t<T-dt/2.; t+=dt)
		{	
			if(myrank==0) cout << "time step: " << t/dt <<  endl;
			set_Hk_DOWN(TEMP1[0], Hk_DOWN_LIST, t);
			set_Hk_DOWN(TEMP2[0], Hk_DOWN_LIST, t+dt);							
			for(int m=-m_max; m<m_max+1; m++)
			{
				for(int n=-n_max; n<n_max+1; n++)
				{		
#ifndef NO_OMP    	
			#pragma omp parallel for                   
#endif								
					for(int i=0; i<dim_new; i++)
					{
						for(int j=0; j<dim_new; j++)
						{
							Hk_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP1)[fq(i,j,dim_new)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP2)[fq(i,j,dim_new)])*dt + double(m)*w_peierls*delta(i,j)*delta(m,n)/double(timesteps_F-1);
						}
					}				
				}
			}
		}

		// Diagonalize Floquet Hamiltonian in order to get eigenvalues and eigenvectors		
		diagonalize_F_DOWN(Hk_FLOQUET, evals_FLOQUET);  
	
		for(int jj=0; jj<DIM_FLOQ; jj++)
		{
			BANDS_FLOQUET[fq(k,jj,DIM_FLOQ)] = evals_FLOQUET[jj];
		}	
		// Calculate squared overlap of Floquet eigenstates with eigenstates of eq. Hamiltonian
		for(int i=0; i<DIM_FLOQ; ++i) // select Floquet band
		{
			temp = 0.;
			for(int j=0; j<dim_new; ++j) // sum over vector of length dim_new
			{
				temp += real(Hk_FLOQUET[fq(i,dim_new*m_max+j,DIM_FLOQ)]*conj(Hk_FLOQUET[fq(i,dim_new*m_max+j,DIM_FLOQ)]));
			}		
			OVERLAP_FLOQUET[fq(k,i,DIM_FLOQ)] = temp; 
		}			
	}
	delete TEMP1, TEMP2;	

	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &BANDS_FLOQUET[0], DIM_FLOQ*num_kpoints_path, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &OVERLAP_FLOQUET[0], DIM_FLOQ*num_kpoints_path, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	// Store data
	if(myrank==0)
	{
		ofstream myfile ("DATA/bands_floquet_DOWN_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<DIM_FLOQ; m++)
				{
					myfile << BANDS_FLOQUET[fq(k,m,DIM_FLOQ)] << " " ;
				}
			myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
	if(myrank==0)
	{
		ofstream myfile ("DATA/overlap_floquet_DOWN_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<DIM_FLOQ; m++)
				{
					myfile << OVERLAP_FLOQUET[fq(k,m,DIM_FLOQ)] << " " ;
				}
			myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
}


void Set_Hk_Floquet_DOWN(cvec &Hk_FLOQUET, cvec &dHkdx_FLOQUET, cvec &dHkdy_FLOQUET, vector<cvec*> Hk_DOWN_LIST)
/**
 *	Set Floquet Hamiltonian in k-orbital basis for use in FLOQUET_BC_LOOP()
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store Floquet Hamiltonian matrix
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store kx derivative of Floquet Hamiltonian matrix
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store ky derivative of Floquet Hamiltonian matrix
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 */
{
	// PROBLME: Calculate Hk at loop points!!!
	const double T = 2.*M_PI/w_peierls;
	const double dt = T/double(timesteps_F-1);
	
	cvec *TEMP1 = new cvec(dim_new*dim_new);
	cvec *TEMP2 = new cvec(dim_new*dim_new); 
	cvec *TEMP3 = new cvec(dim_new*dim_new);
	cvec *TEMP4 = new cvec(dim_new*dim_new); 
	cvec *TEMP5 = new cvec(dim_new*dim_new);
	cvec *TEMP6 = new cvec(dim_new*dim_new); 
	double temp; 
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
#ifndef NO_OMP    	
	#pragma omp parallel for collapse(4)                                // PERFEKTLY nested loops are collapsed into one loop
#endif						                
		for(int m=-m_max; m<m_max+1; m++)
		{
			for(int n=-n_max; n<n_max+1; n++)
			{					
				for(int i=0; i<dim_new; i++)
				{
					for(int j=0; j<dim_new; j++)
					{
						Hk_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] = 0.0;
					}
				}
			}
		}												
		for(double t=0; t<T-dt/2.; t+=dt)
		{	
			if(rank==0) cout << "time step: " << t/dt <<  endl;
			set_Hk_DOWN(TEMP1[0], Hk_DOWN_LIST, t);
			set_Hk_DOWN(TEMP2[0], Hk_DOWN_LIST, t+dt);	
			set_dHkdAx_DOWN(TEMP3[0], Hk_DOWN_LIST, t);
			set_dHkdAx_DOWN(TEMP4[0], Hk_DOWN_LIST, t+dt);
			set_dHkdAy_DOWN(TEMP5[0], Hk_DOWN_LIST, t);
			set_dHkdAy_DOWN(TEMP6[0], Hk_DOWN_LIST, t+dt);								
			for(int m=-m_max; m<m_max+1; m++)
			{
				for(int n=-n_max; n<n_max+1; n++)
				{		
#ifndef NO_OMP    	
			#pragma omp parallel for                   
#endif								
					for(int i=0; i<dim_new; i++)
					{
						for(int j=0; j<dim_new; j++)
						{
							Hk_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP1)[fq(i,j,dim_new)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP2)[fq(i,j,dim_new)])*dt + double(m)*w_peierls*delta(i,j)*delta(m,n)/double(timesteps_F-1);
							dHkdx_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP3)[fq(i,j,dim_new)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP4)[fq(i,j,dim_new)])*dt;
							dHkdy_FLOQUET[f_FL_DOWN(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP5)[fq(i,j,dim_new)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP6)[fq(i,j,dim_new)])*dt;
						}
					}				
				}
			}
			
		}	 				
	delete TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6;	
}


void FLOQUET_BC_DOWN(int k, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<cvec*> Hk_DOWN_LIST, const string &filename)
/** 
 * 	Calculate quantum geoemtric tensor of expanded Floquet Hamiltonian at k (unit is Angstroem^2)
 *  -k:integer to pick k-vec from array
 *  -evals_FLOQUET: Vector[(2*m_max+1)xdim_new] to store Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vector array to stroe Quantum geoemtric tensor
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *	-filename: String to save data
 */
{
	int DIM_FLOQ = (2*n_max+1)*dim_new;
	double temp;	
	int count = COUNT;
		
	dvec evalsx(DIM_FLOQ);
	dvec evalsy(DIM_FLOQ);
	dvec OVERLAP_00(DIM_FLOQ);
	
	cvec *TEMP0 = new cvec(DIM_FLOQ*DIM_FLOQ);
	cvec *TEMP1 = new cvec(DIM_FLOQ*DIM_FLOQ); 
	cvec *TEMP2 = new cvec(DIM_FLOQ*DIM_FLOQ);
	cvec *TEMP_PROD = new cvec(DIM_FLOQ*DIM_FLOQ);
	
	// Calculate eigenvectors of gridpoints along loop

	ReadInMAT(Hk_DOWN_LIST,"HK_TAYLOR/HK_TAYLOR_PATH_"+to_string(k)+".dat");
	Set_Hk_Floquet_DOWN(TEMP0[0], TEMP1[0], TEMP2[0], Hk_DOWN_LIST);
	diagonalize_F_DOWN(TEMP0[0], evals_FLOQUET);	
	
	// Calculate spectral weight in zero photon sector
	for(int i=0; i<DIM_FLOQ; ++i) // select Floquet band
	{
		temp = 0.;
		for(int j=0; j<dim_new; ++j) // sum over vector of length dim_new
		{
			temp += real((*TEMP0)[fq(i,dim_new*m_max+j,DIM_FLOQ)]*conj((*TEMP0)[fq(i,dim_new*m_max+j,DIM_FLOQ)]));
		}		
		OVERLAP_00[i] = temp; 
	}
	
	// Transform Floquet-Current matrices to Floquet energy eigenbasis
	// Ax
	times_nd(TEMP1[0], TEMP0[0], TEMP_PROD[0]);	
	times(TEMP0[0], TEMP_PROD[0], TEMP1[0]);
	// Ay
	times_nd(TEMP2[0], TEMP0[0], TEMP_PROD[0]);	
	times(TEMP0[0], TEMP_PROD[0], TEMP2[0]);
	
	
	// Calculate Phase around loop
	for(int m=0; m<DIM_FLOQ; m++)	
	{
		for(int j=0; j<8; j++){
			bands_BCs_FLOQUET[j][m] = 0.;
		}	
		for(int n=0; n<DIM_FLOQ; n++)
		{
			if(m==n)
			{
				continue;
			}	
			bands_BCs_FLOQUET[0][m] += real((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[1][m] += real((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[2][m] += real((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[3][m] += real((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[4][m] += -2.*imag((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[5][m] += -2.*imag((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[6][m] += -2.*imag((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[7][m] += -2.*imag((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
		}			
	}	
		
	if(filename!="no_file")
	{
		ofstream myfile1 (filename);
		if (myfile1.is_open())
		{
			for(int n=0; n<DIM_FLOQ; ++n) 
			{
				//  Berry curvature equal to phase diveded by area of loop
				for(int j=0; j<8; j++){
					myfile1 << bands_BCs_FLOQUET[j][n] << " ";
				}	
				myfile1 << endl;
			}	
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	delete TEMP0;
	delete TEMP1;
	delete TEMP2; 
	delete TEMP_PROD;
}	


void FLOQUET_BC_PATH_DOWN(vector<dvec> &K_PATH, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<cvec*> Hk_DOWN_LIST, int &numprocs, int &myrank)
/** 
 * 	Calculate quantum geoemtric tensor of expanded Floquet Hamiltonian for k-array (unit is Angstroem^2)
 *  -K_PATH: vector of high-symmetry path vectors
 *  -evals_FLOQUET: Vector[(2*m_max+1)xdim_new] to store Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vector array to stroe Quantum geoemtric tensor
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
{
	int DIM_FLOQ = (2*n_max+1)*dim_new;
	int num_kpoints = K_PATH.size();
	dvec GE_ARRAY(4*num_kpoints*DIM_FLOQ );       
	dvec BC_ARRAY(4*num_kpoints*DIM_FLOQ );                                    
	
	for(int k=myrank; k<num_kpoints; k+=numprocs)
	{
		FLOQUET_BC_DOWN(k, evals_FLOQUET, bands_BCs_FLOQUET, Hk_DOWN_LIST, "no_file");
		for(int n=0; n<DIM_FLOQ ; ++n) 
		{
			for(int j=0; j<4; j++)
			{
				GE_ARRAY[fq(k,n+j*DIM_FLOQ,4*DIM_FLOQ)] = bands_BCs_FLOQUET[j][n];
				BC_ARRAY[fq(k,n+j*DIM_FLOQ,4*DIM_FLOQ)] = bands_BCs_FLOQUET[j+4][n];
			}	
		}	
	}	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &GE_ARRAY[0], 4*num_kpoints*DIM_FLOQ , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &BC_ARRAY[0], 4*num_kpoints*DIM_FLOQ , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	int count = COUNT;
	if(myrank==0)
	{	
		ofstream myfile1 ("DATA/FLOQUET_GE_PATH_DOWN_"+to_string(count)+".dat");
		if (myfile1.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*DIM_FLOQ ; ++n) 
				{
					myfile1 << GE_ARRAY[fq(k,n,4*DIM_FLOQ)] << " ";
				}	
				myfile1 << endl;
			}
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile1 ("DATA/FLOQUET_BC_PATH_DOWN_"+to_string(count)+".dat");
		if (myfile1.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*DIM_FLOQ; ++n) 
				{
					myfile1 << BC_ARRAY[fq(k,n,4*DIM_FLOQ)] << " ";
				}	
				myfile1 << endl;
			}
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
}	


void Hk_bands_Floquet(const double &lambda, dvec &BANDS_FLOQUET, dvec &OVERLAP_FLOQUET, cvec &Hk_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &K_PATH, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, int &numprocs, int &myrank)
/**
 *	Calculate Floquet bands with full single-particle Hmailtonian
 *  -lambda: renormalization constant
 *  -BANDS_FLOQUET: Real vector to store Floquet eigenvalues of all k-points 
 *  -OVERLAP_FLOQUET: Real vector[num_kpoints_PATHxNATOMx(2*n_max+1)] to store overlap ov Flquet bands with equilibrium bands
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store Flqoeut Hamiltonian matrix
 *  -evals_FLOQUET: Real vector[(M_max+1) x NATOM] to store Floquet eigenvalues
 *  -K_PATH: vector of high-symmetry path vectors
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
{
	const int DIM_FLOQ = (2*n_max+1)*NATOM;
	const int num_kpoints_path = K_PATH.size();
	const double T = 2.*M_PI/w_peierls;
	const double dt = T/double(timesteps_F-1);
	
	cvec *TEMP1 = new cvec(NATOM*NATOM);
	cvec *TEMP2 = new cvec(NATOM*NATOM); 
	double temp; 
	
	for(int k=myrank; k<num_kpoints_path; k+=numprocs)
	{
#ifndef NO_OMP    	
	#pragma omp parallel for collapse(4)                                
#endif						                
		for(int m=-m_max; m<m_max+1; m++)
		{
			for(int n=-n_max; n<n_max+1; n++)
			{					
				for(int i=0; i<NATOM; i++)
				{
					for(int j=0; j<NATOM; j++)
					{
						Hk_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] = 0.0;
					}
				}
			}
		}												
		if(myrank==0) cout << endl; 
		if(myrank==0) cout << "k = " << k << endl;	
		// Perform integration over one period T		
		for(double t=0; t<T-dt/2.; t+=dt)
		{	
			if(myrank==0) cout << "time step: " << t/dt <<  endl;
			
			set_Hk(lambda, K_PATH[k], TEMP1[0], lvec, UNIT_CELL, t);
			set_Hk(lambda, K_PATH[k], TEMP2[0], lvec, UNIT_CELL, t+dt);								
			for(int m=-m_max; m<m_max+1; m++)
			{
				for(int n=-n_max; n<n_max+1; n++)
				{		
#ifndef NO_OMP    	
			#pragma omp parallel for                   
#endif								
					for(int i=0; i<NATOM; i++)
					{
						for(int j=0; j<NATOM; j++)
						{
							Hk_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(II*w_peierls*double(m-n)*t)*(*TEMP1)[fq(i,j,NATOM)] + exp(II*w_peierls*double(m-n)*(t+dt))*(*TEMP2)[fq(i,j,NATOM)])*dt + double(m)*w_peierls*delta(i,j)*delta(m,n)/double(timesteps_F-1);
						}
					}				
				}
			}
		}
		// Diagonalize Floquet Hamiltonian in order to get eigenvalues and eigenvectors		
		diagonalize_F(Hk_FLOQUET, evals_FLOQUET);  		
		for(int jj=0; jj<NATOM*(2*n_max+1); jj++)
		{
			BANDS_FLOQUET[fq(k,jj,NATOM*(2*n_max+1))] = evals_FLOQUET[jj];
		}	
		// Calculate squared overlap of Floquet eigenstates with eigenstates of eq. Hamiltonian
		for(int i=0; i<DIM_FLOQ; ++i) // select Floquet band
		{
			temp = 0.;
			for(int j=0; j<NATOM; ++j) // sum over vector of length dim_new
			{
				temp += real(Hk_FLOQUET[fq(i,NATOM*m_max+j,DIM_FLOQ)]*conj(Hk_FLOQUET[fq(i,NATOM*m_max+j,DIM_FLOQ)]));
			}		
			OVERLAP_FLOQUET[fq(k,i,DIM_FLOQ)] = temp; 
		}	  				
	}
	delete TEMP1, TEMP2;	
	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &BANDS_FLOQUET[0], NATOM*(2*n_max+1)*num_kpoints_path, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &OVERLAP_FLOQUET[0], NATOM*(2*n_max+1)*num_kpoints_path, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	// Store data
	int count = COUNT;
	if(myrank==0)
	{
		ofstream myfile ("DATA/bands_floquet_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<NATOM*(2*n_max+1); m++)
				{
					myfile << BANDS_FLOQUET[fq(k,m,NATOM*(2*n_max+1))] << " " ;
				}
			myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
	if(myrank==0)
	{
		ofstream myfile ("DATA/overlap_floquet_"+to_string(count)+".dat");
		if (myfile.is_open())
		{
			for(int k=0; k<num_kpoints_path; k++)
			{
				for(int m=0; m<NATOM*(2*n_max+1); m++)
				{
					myfile << OVERLAP_FLOQUET[fq(k,m,NATOM*(2*n_max+1))] << " " ;
				}
			myfile  << endl;
			}
		myfile.close();
		}
		else cout << "Unable to open file" << endl;
	}
}


void Set_Hk_Floquet(dvec &kvec, const double &lambda, cvec &Hk_FLOQUET, cvec &dHkdx_FLOQUET, cvec &dHkdy_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec)
/**
 *	Set Floquet Hamiltonian in k-orbital basis for use in FLOQUET_BC_LOOP()
 * 	-kvec: Real vector of the reciprocal space
 *  -lambda: renormalisation constant
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store Floquet Hamiltonian matrix
 *  -dHkdx_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store kx-derivative of Floquet Hamiltonian matrix
 *  -dHkdy_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store ky-derivative ofF loquet Hamiltonian matrix
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 */
{
	// PROBLME: Calculate Hk at loop points!!!
	const double T = 2.*M_PI/w_peierls;
	const double dt = T/double(timesteps_F-1);

	cvec *TEMP1 = new cvec(NATOM*NATOM);
	cvec *TEMP2 = new cvec(NATOM*NATOM); 
	cvec *TEMP3 = new cvec(NATOM*NATOM);
	cvec *TEMP4 = new cvec(NATOM*NATOM); 
	cvec *TEMP5 = new cvec(NATOM*NATOM);
	cvec *TEMP6 = new cvec(NATOM*NATOM); 
	double temp; 
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
#ifndef NO_OMP    	
	#pragma omp parallel for collapse(4)                                // PERFEKTLY nested loops are collapsed into one loop
#endif						                
		for(int m=-m_max; m<m_max+1; m++)
		{
			for(int n=-n_max; n<n_max+1; n++)
			{					
				for(int i=0; i<NATOM; i++)
				{
					for(int j=0; j<NATOM; j++)
					{
						Hk_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] = 0.0;
					}
				}
			}
		}												
		for(double t=0; t<T-dt/2.; t+=dt)
		{	
			if(rank==0) cout << "time step: " << t/dt <<  endl;
			set_Hk(lambda, kvec, TEMP1[0], lvec, UNIT_CELL, t);
			set_Hk(lambda, kvec, TEMP2[0], lvec, UNIT_CELL, t+dt);	
			set_dHkdAx(lambda, kvec, TEMP3[0], lvec, UNIT_CELL, t);
			set_dHkdAx(lambda, kvec, TEMP4[0], lvec, UNIT_CELL, t+dt);
			set_dHkdAy(lambda, kvec, TEMP5[0], lvec, UNIT_CELL, t);
			set_dHkdAy(lambda, kvec, TEMP6[0], lvec, UNIT_CELL, t+dt);								
			for(int m=-m_max; m<m_max+1; m++)
			{
				for(int n=-n_max; n<n_max+1; n++)
				{		
#ifndef NO_OMP    	
			#pragma omp parallel for                   
#endif								
					for(int i=0; i<NATOM; i++)
					{
						for(int j=0; j<NATOM; j++)
						{
							Hk_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP1)[fq(i,j,NATOM)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP2)[fq(i,j,NATOM)])*dt + double(m)*w_peierls*delta(i,j)*delta(m,n)/double(timesteps_F-1);
							dHkdx_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP3)[fq(i,j,NATOM)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP4)[fq(i,j,NATOM)])*dt;
							dHkdy_FLOQUET[f_FL(m+m_max, n+n_max, i, j)] += 0.5/T*(exp(-II*w_peierls*double(m-n)*t)*(*TEMP5)[fq(i,j,NATOM)] + exp(-II*w_peierls*double(m-n)*(t+dt))*(*TEMP6)[fq(i,j,NATOM)])*dt;
						}
					}				
				}
			}
			
		}	 				
	delete TEMP1, TEMP2, TEMP3, TEMP4, TEMP5, TEMP6;	
}


void FLOQUET_BC(dvec &kvec, const double &lambda, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, const string &filename)
/** 
 * 	Calculate quantum geometric tensor of full Floquet Hamiltonian at kvec
 * 	-kvec: Real vector of the reciprocal space
 *  -lambda: renormalisation constant
 *  -evals_FLOQUET: vector to store of Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vectro of vectros to store quantum geometric tensor 
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 */
{
	int DIM_FLOQ = (2*n_max+1)*NATOM;
	double temp;	
		
	dvec evalsx(DIM_FLOQ);
	dvec evalsy(DIM_FLOQ);
	dvec OVERLAP_00(DIM_FLOQ);
	
	cvec *TEMP0 = new cvec(DIM_FLOQ*DIM_FLOQ);
	cvec *TEMP1 = new cvec(DIM_FLOQ*DIM_FLOQ); 
	cvec *TEMP2 = new cvec(DIM_FLOQ*DIM_FLOQ);
	cvec *TEMP_PROD = new cvec(DIM_FLOQ*DIM_FLOQ);
	
	// Calculate eigenvectors of gridpoints along loop
	Set_Hk_Floquet(kvec, lambda, TEMP0[0], TEMP1[0], TEMP2[0], evals_FLOQUET, UNIT_CELL, lvec);
	diagonalize_F(TEMP0[0], evals_FLOQUET);	
	
	for(int i=0; i<DIM_FLOQ; ++i) // select Floquet band
	{
		temp = 0.;
		for(int j=0; j<NATOM; ++j) // sum over vector of length dim_new
		{
			temp += real((*TEMP0)[fq(i,NATOM*m_max+j,DIM_FLOQ)]*conj((*TEMP0)[fq(i,NATOM*m_max+j,DIM_FLOQ)]));
		}		
		OVERLAP_00[i] = temp; 
	}	
	
	// Transform Floquet-Current matrices to Floquet energy eigenbasis
	// Ax
	times_nd(TEMP1[0], TEMP0[0], TEMP_PROD[0]);	
	times(TEMP0[0], TEMP_PROD[0], TEMP1[0]);
	// Ay
	times_nd(TEMP2[0], TEMP0[0], TEMP_PROD[0]);	
	times(TEMP0[0], TEMP_PROD[0], TEMP2[0]);
	
	// Calculate Phase around loop
	for(int m=0; m<DIM_FLOQ; m++)	
	{
		for(int j=0; j<8; j++){
			bands_BCs_FLOQUET[j][m] = 0.;
		}	
		for(int n=0; n<DIM_FLOQ; n++)
		{
			if(m==n)
			{
				continue;
			}	
			bands_BCs_FLOQUET[0][m] += real((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[1][m] += real((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[2][m] += real((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[3][m] += real((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[4][m] += -2.*imag((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[5][m] += -2.*imag((*TEMP1)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[6][m] += -2.*imag((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP1)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
			bands_BCs_FLOQUET[7][m] += -2.*imag((*TEMP2)[fq(m,n,DIM_FLOQ)]*(*TEMP2)[fq(n,m,DIM_FLOQ)]/((evals_FLOQUET[m]-evals_FLOQUET[n])*(evals_FLOQUET[m]-evals_FLOQUET[n])))*OVERLAP_00[n];
		}			
	}	
		
	if(filename!="no_file")
	{
		ofstream myfile1 (filename);
		if (myfile1.is_open())
		{
			for(int n=0; n<DIM_FLOQ; ++n) 
			{
				//  Berry curvature equal to phase diveded by area of loop
				for(int j=0; j<8; j++){
					myfile1 << bands_BCs_FLOQUET[j][n] << " ";
				}	
				myfile1 << endl;
			}	
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	delete TEMP0;
	delete TEMP1;
	delete TEMP2; 
	delete TEMP_PROD;
}	


void FLOQUET_BC_PATH(const double &lambda, vector<dvec> &K_PATH, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, int &numprocs, int &myrank)
/** 
 * 	Calculate quantum geometric tensor of full Floquet Hamiltonian for momentum array
 *  -lambda: renormalisation constant
 *  -K_PATH: Vector of high-symmetry path vectors
 *  -evals_FLOQUET: vector to store of Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vectro of vectros to store quantum geometric tensor 
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 *  -mprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
{
	int num_kpoints = K_PATH.size();
	int DIM_FLOQ = (2*n_max+1)*NATOM;
	dvec GE_ARRAY(4*num_kpoints*DIM_FLOQ);       
	dvec BC_ARRAY(4*num_kpoints*DIM_FLOQ);                                    
	dvec kvec;
	
	for(int k=myrank; k<num_kpoints; k+=numprocs)
	{
		kvec = K_PATH[k];
		FLOQUET_BC(kvec, lambda, evals_FLOQUET, bands_BCs_FLOQUET, UNIT_CELL, lvec, "no_file");
		for(int n=0; n<DIM_FLOQ; ++n) 
		{
			for(int j=0; j<4; j++){
				GE_ARRAY[fq(k,n+j*DIM_FLOQ,4*DIM_FLOQ)] = bands_BCs_FLOQUET[j][n];
				BC_ARRAY[fq(k,n+j*DIM_FLOQ,4*DIM_FLOQ)] = bands_BCs_FLOQUET[j+4][n];
			}	
		}	
	}	
#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &GE_ARRAY[0], 4*num_kpoints*DIM_FLOQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &BC_ARRAY[0], 4*num_kpoints*DIM_FLOQ, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	int count = COUNT;
	if(myrank==0)
	{	
		ofstream myfile1 ("DATA/FLOQUET_GE_PATH_"+to_string(count)+".dat");
		if (myfile1.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*DIM_FLOQ; ++n) 
				{
					myfile1 << GE_ARRAY[fq(k,n,4*DIM_FLOQ)] << " ";
				}	
				myfile1 << endl;
			}
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	if(myrank==0)
	{	
		ofstream myfile1 ("DATA/FLOQUET_BC_PATH_"+to_string(count)+".dat");
		if (myfile1.is_open())
		{
			for(int k=0; k<num_kpoints; ++k)
			{
				for(int n=0; n<4*DIM_FLOQ; ++n) 
				{
					myfile1 << BC_ARRAY[fq(k,n,4*DIM_FLOQ)] << " ";
				}	
				myfile1 << endl;
			}
			myfile1.close();
		}
		else cout << "Unable to open file" << endl;	
	}
}	

