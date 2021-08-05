# TBG_LMC
This project aims to calculate several observables for ultrafast spectroscopic experiments on twisted bilayer graphene. An efficient calculation for small angles with huge supercells is ensured by the possibility to truncate the Hamiltonian matrix in the eigen-energy basis systematically. Depending on the initial dimension of the Hamiltonian (set by the twist angle), hybrid parallel computations with OpenMP and MPI is possible. 

Author:
    [Gabriel E. Topp](fizztopp@gmail.com)

This code was used to generate the results presented in [Gabriel E Topp, Christian J Eckhardt, Dante M Kennes, Michael A Sentef, Päivi Törmä, arXiv:2103.04967](https://arxiv.org/abs/2103.04967).

The main code is written in [C++](https://isocpp.org/). The construction of the atomic lattice / k-vectors and plotting is done with [Python](https://www.python.org/).

Manual to use the code:

1) CONSTRUCT_MODEL_GEOMETRY.py has to be executed. For a choice of commensurate twist angle (defined via parameters m,n), the program computes the real space Bravais vectors (L_VECS.dat), the atomic lattice (Unit_Cell.dat) as well as the k-vectors of the 1st BZ (k_BZ.dat).

2) The main code conists three source files: TBG_LMC.cpp (main), StandardFunctions.cpp(standard mathematical operations), and Functions.cpp(project specific functions). Constants.h sets all model and simulation parameters.
   The truncation of the original Hamiltonian matrix is done by functions set_Hk_DOWN_LIST() and alc_List() from  StandardFunctions.cpp. Importantly, the systematic truncation is controlled via two parameters in Constants.h. dim_MAT defines the leading order of the truncated and stored matrices. dim_new can introduce a second truncation for reading matrices from disk to memory, which is convenient for convergence checks. 
   After truncation, different observables like density-of-states, linear optical conductivity, non-linear optical conductivity, quantum geometric tensor, Floquet bands can be calculated.
   
3) The C++ programm can be compiled via the CMake file CMakeLists.txt.

4) job.sge is an example batch file to run the programm on the [TRITON cluster](https://scicomp.aalto.fi/triton/).

4) File PLOTTING.py is an example plotting file that was used to create the output shown in Fig.5 of https://arxiv.org/abs/2103.04967.
