program DSM3D
  ! 
  ! 
  !   Yes, this will be the first 3D DSM code for 3D Earth with 1D TI normal mode cards and
  !   3D perturbation to them written in vtk file
  ! 
  ! 
  !                                     Nobuaki Fuji 2017 
  !                                     Institut de Physique du Globe de Paris
  !


  use mpi
  use parameters
  use inputFiles
  implicit none

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  
  if(my_rank.eq.0) then
     call pinput
  
