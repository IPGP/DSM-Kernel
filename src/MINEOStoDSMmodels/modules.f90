module parameters

  implicit none
  
  ! MINEOS
  
  character(120) :: modelid,model1d
  integer :: ifanis,ifdeck
  double precision :: tref
  integer :: nrmod,nicb,ncmb,nmspl

  double precision, allocatable, dimension(:) :: rmod,dnm,vpv,vsv
  double precision, allocatable, dimension(:) :: qk,qm,vph,vsh,eta
  

  ! DSM

  character(120) :: psvmodel
  integer :: nzone
  double precision, allocatable, dimension (:) :: vrmin, vrmax,qmuD,qkappaD
  integer, allocatable, dimension (:) :: irmin, irmax,irlength
  double precision, allocatable, dimension (:,:) :: rrhoD,vpvD,vphD,vsvD,vshD,etaD


end module parameters
  
