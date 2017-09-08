#!/bin/sh
#SBATCH --job-name DSM-ker
#SBATCH --nodes 8 
#SBATCH --ntasks-per-node=16
#SBATCH --partition cpunormal
#SBATCH --exclusive
module purge
#module load intel/16.0.1 bullxmpi
module load intel/compiler intel/mpi intel/mkl slurm




ulimit -s unlimited
ulimit -a

srun    ../../bin/mpiKernelMaker < ./kernel.inf
