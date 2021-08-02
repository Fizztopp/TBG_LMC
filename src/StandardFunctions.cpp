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


void ReadIn(vector<dvec> &MAT, const string& filename)
{
/**
 *	Read in real valued matrix
 */
	ifstream in(filename);
	string record;
	if(in.fail()){
		cout << "file" << filename << "could not be found!" << endl;
	}
	while (getline(in, record))
	{
		istringstream is( record );
		dvec row((istream_iterator<double>(is)),
		istream_iterator<double>());
		MAT.push_back(row);
	}
	in.close();
}



extern "C" {
/** 
 *  Computes the eigenvalues and, optionally, the eigenvectors for a Hermitian matrices H
 */
    void zheev_(char* jobz, char* uplo, int* N, cdouble* H, int* LDA, double* W, cdouble* work, int* lwork, double* rwork, int *info);
}

//'N','V':  Compute eigenvalues only, and eigenvectors
char    jobz = 'V';       
//'U','L':  Upper, Lower triangle of H is stored 
char    uplo = 'U';  
// The order of the matrix H.  NATOM >= 0
int     matsize = NATOM;    
// The leading dimension of the array H.  lda >= max(1, NATOM)
int     lda = NATOM;             
// The length of the array work.  lwork  >= max(1,2* NATOM-1)
int     lwork = 2*NATOM-1;    
// dimension (max(1, 3* NATOM-2))
double  rwork[3*NATOM-2];  
// dimension (MAX(1,LWORK))
cdouble work[2*NATOM-1];  
// Info
int	    info;



void diagonalize(cvec &Hk, dvec &evals)
{
/**
 *  Diagonalization of matrix Hk. Stores eigenvalues in real vector evals and eigenvectors in complex vector Hk
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */
    zheev_(&jobz, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}



char    jobz_eig = 'N';         
void diagonalize_eig(cvec &Hk, dvec &evals)
{
/**
 *  Diagonalization of matrix Hk. Stores ONLY eigenvalues in real vector evals
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */
    zheev_(&jobz_eig, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}   



int     matsize_BdG = NATOM*2;      									// The order of the matrix A.  N >= 0
int     lda_BdG = NATOM*2;            									// The leading dimension of the array A.  LDA >= max(1,N)
int     lwork_BdG = 2*NATOM*2-1;      									// The length of the array WORK.  LWORK >= max(1,2*N-1)
double  rwork_BdG[3*NATOM*2-2];       									// dimension (max(1, 3*N-2))
cdouble work_BdG[2*NATOM*2-1]; 



void diagonalize_BdG(cvec &Hk_BdG, dvec &evals_BdG)
{
/**
 *  Diagonalization of Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM]
 *  -Hk: Complex vector[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[(2*m_max+1) x NATOM] to store eigenvalues
 */
	zheev_(&jobz, &uplo, &matsize_BdG, &Hk_BdG[0], &lda_BdG, &evals_BdG[0], &work_BdG[0], &lwork_BdG, &rwork_BdG[0], &info);
	assert(!info);
}
                      
                      
int     matsize_F = NATOM*(2*n_max+1);      							// The order of the matrix A.  N >= 0
int     lda_F = NATOM*(2*n_max+1);            							// The leading dimension of the array A.  LDA >= max(1,N)
int     lwork_F = 2*NATOM*(2*n_max+1)-1;      							// The length of the array WORK.  LWORK >= max(1,2*N-1)
double  rwork_F[3*NATOM*(2*n_max+1)-2];       							// dimension (max(1, 3*N-2))
cdouble work_F[2*NATOM*(2*n_max+1)-1]; 



void diagonalize_F(cvec &Hk_FLOQUET, dvec &evals_FLOQUET)
{
/**
 *  Diagonalization of Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM]
 *  -Hk: Complex vector[(2*m_max+1) x (2*m_max+1) x NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[(2*m_max+1) x NATOM] to store eigenvalues
 */
	zheev_(&jobz, &uplo, &matsize_F, &Hk_FLOQUET[0], &lda_F, &evals_FLOQUET[0], &work_F[0], &lwork_F, &rwork_F[0], &info);
	assert(!info);
}



extern "C" {
/**
 *  Computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE (non Hermitian) matrix --> used in FloquetEVs
 */ void zgeev_(char* JOBVL, char* JOBVR, int* N, cdouble* A, int* LDA, cdouble* W, cdouble* VL, int* LDVL, cdouble* VR, int* LDVR, cdouble* work, int* lwork, double* rwork, int *info);
}

char JOBVL = 'N';                                                       // don't compute left eigenvectors!  
char JOBVR = 'V';                                                       // compute right eigenvectors! 
int LDVL = NATOM;
int LDVR = NATOM; 
int     lwork_GE = 2*NATOM;      										
double  rwork_GE[2*NATOM]; 



void diagonalize_DOWN(cvec &Hk_DOWN, dvec &evals_DOWN)
{
/**
 *  Diagonalization of truncated matrix Hk_DOWN. Writes eiegenvalues to vector evals and eigenvectors (normalized!) to matrix Hk_DOWN
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian --> transformation matrices
 * 	-evals_DOWN: Real vector[dim_new] to store eigenvalues
 *  -dim_new: integer value of reduced leading order of Hamiltonian
 */
	char    jobz_D = 'V';             									//'N','V':  Compute eigenvalues only/+eigenvectors
	char    uplo_D = 'U';              									//'U','L':  Upper/Lower triangle of H is stored
	int     matsize_D = dim_new;      									// The order of the matrix A.  N >= 0
	int     lda_D = dim_new;            							    // The leading dimension of the array A.  LDA >= max(1,N)
	int     lwork_D = 2*dim_new-1;      								// The length of the array WORK.  LWORK >= max(1,2*N-1)
	double  rwork_D[3*dim_new-2];       							    // dimension (max(1, 3*N-2))
	cdouble work_D[2*dim_new-1];        								// dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK
	int	    info_D;
	zheev_(&jobz_D, &uplo_D, &matsize_D, &Hk_DOWN[0], &lda_D, &evals_DOWN[0], &work_D[0], &lwork_D, &rwork_D[0], &info_D);
	assert(!info_D);
}	



void diagonalize_F_DOWN(cvec &Hk_FLOQUET, dvec &evals_FLOQUET)
{
/**
 *  Diagonalization of truncated Floquet matrix H_FLOQUET[(2*m_max+1) x (2*m_max+1) x dim_new x dim_new]
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1) x (2*m_max+1) x dim_new x dim_new] to store Hamiltonian --> transformation matrices
 * 	-evals_FLOQUET: Real vector[(2*m_max+1) x dim_new] to store eigenvalues
 */
 	char    jobz_F_DOWN = 'V';             								//'N','V':  Compute eigenvalues only/+eigenvectors
	char    uplo_F_DOWN = 'U';              							//'U','L':  Upper/Lower triangle of H is stored
	int     matsize_F_DOWN = dim_new*(2*n_max+1);      					// The order of the matrix A.  N >= 0
	int     lda_F_DOWN = dim_new*(2*n_max+1);            				// The leading dimension of the array A.  LDA >= max(1,N)
	int     lwork_F_DOWN = 2*dim_new*(2*n_max+1)-1;      				// The length of the array WORK.  LWORK >= max(1,2*N-1)
	double  rwork_F_DOWN[3*dim_new*(2*n_max+1)-2];       				// dimension (max(1, 3*N-2))
	cdouble work_F_DOWN[2*dim_new*(2*n_max+1)-1]; 
	int	    info_F_DOWN;
	zheev_(&jobz_F_DOWN, &uplo_F_DOWN, &matsize_F_DOWN, &Hk_FLOQUET[0], &lda_F_DOWN, &evals_FLOQUET[0], &work_F_DOWN[0], &lwork_F_DOWN, &rwork_F_DOWN[0], &info_F_DOWN);
	assert(!info);
}

