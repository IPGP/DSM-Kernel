## single station and event example

### Overview
This example demonstrates the two central steps of this code:

1. the computation of a database of Green functions for a specified region of the Earth.
2. the computation of a kernel using the previously computed database

### Instructions
#### The Green function database
To compute the database you need the two programs SGTpsv and SGTsh. That
compute the databases necessary for PSV and SH waves, respectively. The sources
of these codes are stored in `DSM-Kernel/src/SGTpsv/` and `DSM-Kernel/src/SGTsh/`

##### Compilation
1. adjust the Makefiles (e.g. `DSM-Kernel/src/SGTpsv/Makefile`). You should for
   example chose you fortran compiler, e.g. `FC=mpif90` or `FC=mpiifort`.
2. type: `make` in the SGTsh and in the SGTpsv directories.
3. check that there were no errors and if you can find the binaries:
   `bin/mpiSGTpsv` and `DSM-Kernel/bin/mpiSGTsh`

##### Configuration and Running
1. copy the file `database.inf.example` to `database.inf`. In this file:
 * adjust the paths (replace [xy] with the path of your DSM-Kernel installation)
 * adjust the output directory path [mydatabase], where the database is stored.
   The database needs some disk space. Make sure to have xy Gb free.
 * generate the following subdirectories in your database folder: `[mydatabase]/log`,
   `[mydatabase]/RSGT`and `[mydatabase]/TSGT`
 * you can adjust the remaining parameters that control e.g. frequency range,
   domain size and sampling of the database. In particular chose the source depth
   for which you want to compute the database!
2. submit the command `mpirun -tmi -np 256  ../../bin/mpiSGTpsv < ./database.inf`
   to compute the psv database. Select the number of processors that you want to use
   with the -np parameter (256 in this example). `run_database.slurm` shows an
   example script that submits this process to the IPGP cluster malbec.

#### Kernel Calculation
1. copy the file `kernel.inf.example` to `kernel.inf`. In this file:
 * adjust the directory `[xy]` to point to your DSM-Kernels installation
 * check out the other parameters that are used to control the source, filters
   and other parameters. The time window controls for which phase you compute
   the sensitivity kernel. The source depth has to correspond to the one specified
   in the database configuration (`database.inf`).

#### Visualization
