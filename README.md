# Genetic Algorithm
A Genetic Algorithm code for global optimization of nano clusters
This code is used for global optimization of the geometry of a atom cluster, given the number of atoms, with respect to the total energy 
*Parallel library MPI was used for parallel computing
*Either Lammps or DFTB package can be used for energy evaluation
*Efficiency and robustness of this algorithm has been benchmarked by flowing papers:

Shi, H., Auerbach, S. M, et al. "First-Principles Predictions of Structure Function Relationships of Graphene-Supported Platinum Nanoclusters". The Journal of Physical Chemistry C. 2016, 120, 11899-11909.


Input file: GA_input
Source file: cal.c, Lammps.c, read.c

How to compile.
Make sure you have LAMMPS or DFTB package installed so related libraries can be linked  
./compile.sh in Linux environment
