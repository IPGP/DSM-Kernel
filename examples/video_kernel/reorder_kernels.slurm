#!/bin/sh
#SBATCH --job-name kernelconverter
#SBATCH --nodes 1
#SBATCH --partition cpunormal
#SBATCH --exclusive

module purge
module load intel/compiler intel-mpi/64 intel/mkl slurm
/home/meschede/gpfs/DSM_KERNELS/RUN/combine_kernels
