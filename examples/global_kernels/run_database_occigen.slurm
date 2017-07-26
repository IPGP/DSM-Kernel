#!/bin/sh
#SBATCH --job-name DSM-db
#SBATCH --nodes=8 
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --threads-per-core=1
#SBATCH --time=08:20:00
##SBATCH --mem=120GB
#SBATCH --constraint="HSW24"

module purge
#module load intel/16.0.1 bullxmpi
#module load intel/compiler intel-mpi/64 intel/mkl slurm
module load intel intelmpi slurm


ulimit -s unlimited
ulimit -a

mpirun -tmi -np 128  ../../bin/mpiSGTpsv < ./database.inf 
mpirun -tmi -np 128  ../../bin/mpiSGTsh < ./database.inf
