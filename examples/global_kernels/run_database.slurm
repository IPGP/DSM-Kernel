#!/bin/sh
#SBATCH --job-name DSM-db
#SBATCH --nodes=8 
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --threads-per-core=1
#SBATCH --time=08:20:00
#SBATCH --partition=cpunormal,data,datanew
##SBATCH --mem=120GB

module purge
#module load intel/16.0.1 bullxmpi
module load intel/compiler intel/mpi intel/mkl slurm




ulimit -s unlimited
ulimit -a

srun    ../../bin/mpiSGTpsv < ./database.inf 
