#ifndef TBG_LMC_Functions_H
#define TBG_LMC_Functions_H



//Void FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void set_Hk0(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL);
/**
 * 	Set eq. single-particle Hamiltonian (without external field)
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 */
 


void set_Hk(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time);
/**
 * 	Set time-dependent single-particle Hamiltonian with external field
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */



void set_dHkdAx(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time);
/**
 * 	Set kx derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */



void set_dHkdAy(const double &lambda, dvec &kvec, cvec &Hk, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, double time);
/**
 * 	Set ky derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian
 *  -Real basis vector
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -time: real time parameter
 */


void set_Hk_Taylor(const double &lambda, dvec &kvec, vector<cvec*> Hk_Taylor, vector<dvec> &lvec, vector<dvec> &UNIT_CELL);
/**
 *	Taylor expansion of single-particle Hamiltonian for small fields
 * 	Set kx derivative of time-dependent single-particle Hamiltonian with external field in
 *  -lambda: renormalisation constant
 *  -kvec: Real vector of the reciprocal space
 *  -Hk_Taylor: Vector of complex matrices[10][NATOM*NATOM] to store Taylor matrices
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 */



void set_Hk_DOWN_LIST(const double &lambda, dvec &kvec, dvec &evals, vector<cvec*> Hk_Taylor, vector<cvec*> Hk_DOWN_LIST, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, int &myrank);
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



void Calc_List(const double &lambda, dvec &evals, cvec &Hk, vector<dvec> &K_POINTS, vector<dvec> &lvec, vector<dvec> &UNIT_CELL, vector<cvec*> Hk_Taylor, vector<cvec*> Hk_DOWN_LIST, const string& switcher, int &numprocs, int &myrank);
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



void ReadInMAT(vector<cvec*> Hk_DOWN_LIST, const string& filename);
/**
  * Read in taylore matrices from disc
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  *	-filename: String to define file
  */



void Hk_bands_DOWN(vector<dvec> &K_PATH, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank);
/**
 *	Calculate bands of downfolded Hamiltonian path K_PATH through BZ 
 *  -K_PATH: vector[NPATH] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */



void FUNC_DOS(cdouble &omega, double &mu, double &DOS, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename);
/**
 *	Calulates density-of-states for one frequency
 *  -omega: frequency
 *  -mu: chemical potential
 *  -DOS: denisty-of-states
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 */



void DOS_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, double &mu, double &DOS, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank);
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



void LIN_COND(cdouble &omega, double &mu, cvec &COND_O1, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename);
/**
 *	Calulates linear conductivity for one frequency
 *  -omega: frequency
 *  -mu: chemical potential
 *  -COND_O1: linear conductivity
 *  -BZ_FULL: vector[NBZ] of real vectors[3] to store k-points 
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 * 	-filename: String to store data
 */



void LIN_COND_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, double &mu, cvec &COND_O1, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank);
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



void O2_COND(cdouble &omega1, cdouble &omega2, double &mu, cvec &COND_O2, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename);
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



void O2_COND_OMEGA(cdouble domega, cdouble omega_min, cdouble omega_max, cdouble omega2, double &mu, cvec &COND_O2, vector<dvec> &BZ_FULL, vector<cvec*> Hk_DOWN_LIST, const string& filename, int &numprocs, int &myrank);
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



void EQ_BC(int k, vector<dvec> &OFF_BC, vector<dvec> &bands_BCs, vector<cvec*> Hk_DOWN_LIST, const string& filename, const string& switcher);
/**
 * 	Calculate quantum geometric tensor (QGT) in equlibrium for single k-point by derivatives of Hamiltonian
 *  -k: integer to pick k-vector from array
 *  -OFF_BC: matric to store interband terms
 *  -bands_BCs: matrix to store intraband terms
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -filename: string to name file
 *  -switcher: string to switch between k-path and BZ
 */



void EQ_BC_PATH(vector<dvec> &K_POINTS, vector<dvec> &bands_BCs, vector<cvec*> Hk_DOWN_LIST, const string& switcher, int &numprocs, int &myrank);
/**
 * 	Calculate quantum geometric tensor (QGT) in equlibrium for momuntum array
 *  -K_POINTS: vector[NBZ] of real vectors[3] to store k-points 
 *  -bands_BCs: matrix to store intraband terms
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -switcher: string to switch between k-path and BZ
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */



void HF_EFF(vector<cvec*> Hk_DOWN_LIST, vector<dvec> &K_PATH, const string& filename, const string& switcher, int &numprocs, int &myrank);
/**
 * 	Calculate the effective Floquethamiltonian in the high frequency expansion
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *  -K_PATH: vector[NBZ] of real vectors[3] to store k-points 
 *  -filename: string to name file
 *  -switcher: string to switch between k-path and BZ
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI)
 */



void Hk_bands_SP(const double &lambda, dvec &BANDS, cvec &Hk, dvec &evals, vector<dvec> &K_PATH, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, const string& filename, int &numprocs, int &myrank);
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



void mu_SP(const double &lambda, cvec &Hk, dvec &evals, vector<dvec> &kweights, vector<dvec> &BZ, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, double &mu, int &numprocs, int &myrank);
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



void set_Hk_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time);
/**
  *	Set downfolded td Hamiltonian 
  * -dim_new: integer value of reduced leading order of Hamiltonian
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -ASD: Gauge field of source-drain field
  * -Pol: double to set chirality
  * -time: tiome variable
  */



void set_dHkdAx_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time);
/**
  * Set downfolded t.-d. derivative by Ax of Hamiltonian (original band basis)
 *  -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -time: tiome variable
  */



void set_dHkdAy_DOWN(cvec &Hk_DOWN, vector<cvec*> Hk_DOWN_LIST, double time);
/**
  * Set downfolded t.-d. derivative by Ay of Hamiltonian (original band basis)
  * -Hk_DOWN: Complex vector[dim_new x dim_new] to store Hamiltonian matrix
  * -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
  * -time: tiome variable
  */
	


void Hk_bands_Floquet_DOWN(dvec &BANDS_FLOQUET, dvec &OVERLAP_FLOQUET, cvec &Hk_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &K_PATH, vector<cvec*> Hk_DOWN_LIST, int &numprocs, int &myrank);
/**
 *	Set Floquet Hamiltonian in k-orbital basis for use in FLOQUET_BC_LOOP()
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store Floquet Hamiltonian matrix
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store kx derivative of Floquet Hamiltonian matrix
 *  -Hk_FLOQUET: Complex vector[(2*m_max+1)x(2*n_max+1)xNATOMxNATOM] to store ky derivative of Floquet Hamiltonian matrix
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 */



void FLOQUET_BC_DOWN(int k, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<cvec*> Hk_DOWN_LIST, const string &filename);
/** 
 * 	Calculate quantum geoemtric tensor of expanded Floquet Hamiltonian at k (unit is Angstroem^2)
 *  -k:integer to pick k-vec from array
 *  -evals_FLOQUET: Vector[(2*m_max+1)xdim_new] to store Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vector array to stroe Quantum geoemtric tensor
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *	-filename: String to save data
 */



void FLOQUET_BC_PATH_DOWN(vector<dvec> &K_PATH, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<cvec*> Hk_DOWN_LIST, int &numprocs, int &myrank);
/** 
 * 	Calculate quantum geoemtric tensor of expanded Floquet Hamiltonian for k-array (unit is Angstroem^2)
 *  -K_PATH: vector of high-symmetry path vectors
 *  -evals_FLOQUET: Vector[(2*m_max+1)xdim_new] to store Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vector array to stroe Quantum geoemtric tensor
 *  -Hk_DOWN_LIST: Vector of complex matrices[10][dim_new x dim_new] to store truncated Taylor matrices in initial band basis
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
	


void Hk_bands_Floquet(const double &lambda, dvec &BANDS_FLOQUET, dvec &OVERLAP_FLOQUET, cvec &Hk_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &K_PATH, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, int &numprocs, int &myrank);
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



void Set_Hk_Floquet(dvec &kvec, const double &lambda, cvec &Hk_FLOQUET, cvec &dHkdx_FLOQUET, cvec &dHkdy_FLOQUET, dvec &evals_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec);
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



void FLOQUET_BC(dvec &kvec, const double &lambda, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, const string &filename);
/** 
 * 	Calculate quantum geometric tensor of full Floquet Hamiltonian at kvec
 * 	-kvec: Real vector of the reciprocal space
 *  -lambda: renormalisation constant
 *  -evals_FLOQUET: vector to store of Floquet eigenvalues
 *  -bands_BCs_FLOQUET: vectro of vectros to store quantum geometric tensor 
 *  -UNIT_CELL: Vector[NATOM] of real vectors[4] containing atomic positions and sublattice info
 *  -lvec: Real vector[4] of superlattice bravis translational vectors (in lconst*Angstroem)
 */



void FLOQUET_BC_PATH(const double &lambda, vector<dvec> &K_PATH, dvec &evals_FLOQUET, vector<dvec> &bands_BCs_FLOQUET, vector<dvec> &UNIT_CELL, vector<dvec> &lvec, int &numprocs, int &myrank);
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
	


#endif //TBG_LMC_FUNCTIONS_H
