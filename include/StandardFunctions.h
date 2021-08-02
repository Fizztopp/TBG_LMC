#ifndef TBG_LMC_StandardFunctions_H
#define TBG_LMC_StandardFunctions_H

#include "Constants.h"

//INLINE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline int fq(int i, int j, int N)
/**
 *  MAT[i,j] = Vec[fq(i,j,N)] with row index i and column index j
 */
{
    return i*N+j;
}



inline int f_FL(int m, int n, int i, int j)
/**
 *	Wrapper for Floquet matrix MAT[m, n, i, j], (2*m_max+1) x 2*n_max+1) x NATOM x NATOM block matrix element where i,j in {0,..,NATOM-1}, m in {-m,...,0,...+m}, n in {-n,...,0,...+n}
 */
{
	return (2*n_max+1)*NATOM*NATOM*m + NATOM*n + (2*n_max+1)*NATOM*i + j;
}



inline int f_FL_DOWN(int m, int n, int i, int j)
/**
 *	Wrapper for truncated Floquet matrix MAT[m, n, i, j], (2*m_max+1) x (2*n_max+1)  dim_new x dim_new block matrix element where i,j in {0,..,NATOM-1}, m in {-m,...,0,...+m}, n in {-n,...,0,...+n}
 */
{
	return (2*n_max+1)*dim_new*dim_new*m + dim_new*n + (2*n_max+1)*dim_new*i + j;
}



inline double delta(int a, int b)
/**
 *  Delta function
 */
{
	if (a==b)
		return 1.;
	else
		return 0.;
}



template <class Vec>
inline void print(Vec vec)
/**
 *	Print out vector
 */
{
	for(int i=0; i<vec.size(); i++)
		{
	    	cout << vec[i] << " ";
	    }	
	cout << endl;
}



inline double fermi(double energy, double mu)
{
/**
 *	Fermi distribution:
 *	-energy: Energy eigenvalue
 *	-mu: Chemical potential
 */
    return 1./(exp((energy-mu)/TT) + 1.);
}



inline double fermi_zero(double energy, double mu)
{
/**
 *	Fermi distribution:
 *	-energy: Energy eigenvalue
 *	-mu: Chemical potential
 */
	if(energy<mu){
		return 1.0;
	}
	else if(energy==mu){
		return 0.5;
	}
	else{ 
		return 0.0;
	}
}



inline double Ax_t(double time)
{
/**
 *	Peierls field for electrons in x-direction:
 *  -time: Real time coordinate
 */
    return Ax_peierls*sin(w_peierls*time);
}



inline double Ay_t(double time)
{
/**
 *	Peierls field for electrons in y-direction:
 *  -time: Real time coordinate
 */
    return Ay_peierls*cos(w_peierls*time);
}



inline double Az_t(double time)
{
/**
 *	Peierls field for electrons in z-direction:
 *  -time: real time coordinate
 */
    return Az_peierls*sin(w_peierls*time);
}



// Template calss functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template <class Vec>
void times(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product of quadratic matrices: $C = A \cdot B$
 */
{
    int dim = sqrt(A.size());
	Vec TEMP(dim*dim);
    // Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
	    for(int j=0; j<dim; j++) {
		    TEMP[fq(j,i,dim)] = B[fq(i,j,dim)];
		   }
    }
#ifndef NO_OMP 		
	#pragma omp parallel for
#endif 		    
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += A[fq(i,k,dim)]*TEMP[fq(j,k,dim)]; 
			}
		}
	}	
}



template <class Vec>
void times_dn(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of first factor: $C = A^\dagger \cdot B$
 */
{
	int dim = sqrt(A.size());
	Vec TEMP1(dim*dim);
	Vec TEMP2(dim*dim);
	// Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++) {
			TEMP1[fq(j,i,dim)] = A[fq(i,j,dim)];
			TEMP2[fq(j,i,dim)] = B[fq(i,j,dim)];
		}
	}	
#ifndef NO_OMP 		
	#pragma omp parallel for
#endif 			
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += conj(TEMP1[fq(i,k,dim)])*TEMP2[fq(j,k,dim)];
			}
		}
	}		
}


template <class Vec>
void times_nd(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of second factor: $C = A \cdot B^\dagger$
 */
{
	int dim = sqrt(A.size());
#ifndef NO_OMP 		
	#pragma omp parallel for
#endif 				
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
					C[fq(i,j,dim)] += A[fq(i,k,dim)]*conj(B[fq(j,k,dim)]);
			}
		}
	}	
}


//VOID FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ReadIn(vector<dvec> &MAT, const string& filename);
/**
 *	Read in real valued matrix
 */



void diagonalize(cvec &Hk, dvec &evals);
/**
 *  Diagonalization of matrix Hk. Stores eigenvalues in real vector evals and eigenvectors in complex vector Hk
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */
 

        
void diagonalize_eig(cvec &Hk, dvec &evals);
/**
 *  Diagonalization of matrix Hk. Stores ONLY eigenvalues in real vector evals
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */



void diagonalize_BdG(cvec &Hk_BdG, dvec &evals_BdG);
/**
 *  Diagonalization of Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM]
 *  -Hk: Complex vector[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[(2*m_max+1) x NATOM] to store eigenvalues
 */
               
                    

void diagonalize_F(cvec &Hk_FLOQUET, dvec &evals_FLOQUET);
/**
 *  Diagonalization of Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM]
 *  -Hk: Complex vector[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[(2*m_max+1) x NATOM] to store eigenvalues
 */



void diagonalize_DOWN(cvec &Hk_DOWN, dvec &evals_DOWN);
/**
 *  Diagonalization of truncated matrix Hk_DOWN. Writes eiegenvalues to vector evals and eigenvectors (normalized!) to matrix Hk_DOWN
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian --> transformation matrices
 * 	-evals_DOWN: Real vector[dim_new] to store eigenvalues
 *  -dim_new: integer value of reduced leading order of Hamiltonian
 */



void diagonalize_F_DOWN(cvec &Hk_FLOQUET, dvec &evals_FLOQUET);
/**
 *  Diagonalization of truncated Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x dim_new x dim_new]
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1) x (2*m_max+1) x dim_new x dim_new] to store Hamiltonian --> transformation matrices
 * 	-evals_FLOQUET: Real vector[(2*m_max+1) x dim_new] to store eigenvalues
 */


#endif //TBG_LMC_STANDARDFUNCTIONS_H
