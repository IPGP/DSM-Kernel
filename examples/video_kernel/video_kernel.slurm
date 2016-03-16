#!/bin/sh
#SBATCH --job-name dsmkernels
#SBATCH --nodes 16 
#SBATCH --ntasks-per-node=16
#SBATCH --partition cpunormal
#SBATCH --exclusive

module purge
module load intel/compiler intel-mpi/64 intel/mkl slurm
#cd $PBS_O_WORKDIR # on se place dans le repertoire courant
#. $MODULESHOME/init/sh
#module load default-intel-mpi
#cat $PBS_NODEFILE | uniq > mpd.hosts10
#cat $PBS_NODEFILE > nodes10
#mpdboot --rsh=ssh -v -n `cat mpd.hosts10 | wc -l` -f mpd.hosts10
mpirun -tmi -np 256   /gpfs/users/meschede/DSM_KERNELS/DSM-KS2014gamma/bin/mpiKernelMaker < ./matthias-video.inf 
