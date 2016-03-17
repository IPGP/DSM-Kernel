#!/bin/sh
#SBATCH --job-name DSM-db
#SBATCH --nodes 16 
#SBATCH --ntasks-per-node=16
#SBATCH --partition cpunormal
#SBATCH --exclusive

module purge
module load intel/compiler intel-mpi/64 intel/mkl slurm
mpirun -tmi -np 256  ../../bin/mpiSGTpsv < ./database.inf 
mpirun -tmi -np 256  ../../bin/mpiSGTsh < ./database.inf 
