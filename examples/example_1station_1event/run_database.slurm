#!/bin/sh
#SBATCH --job-name deepEarth 
#SBATCH --nodes 16 
#SBATCH --ntasks-per-node=16
#SBATCH --partition cpunormal
#SBATCH --exclusive

module purge
module load intel/compiler intel-mpi/64 intel/mkl slurm
#cd $PBS_O_WORKDIR
#. $MODULESHOME/init/sh
mpirun -tmi -np 256  ../../bin/mpiSGTpsv < ./database.inf 
