program modelConverter
  ! Nobuaki Fuji for DSM Kernel
  !  July 2016, Institut de Physique du Globe de Paris
  
  implicit none
  integer :: i
  
  call getarg(1,model1d)
  model1d = trim(model1d)
  call getarg(2,psvfile)
  psvfile = trim(psvfile)

  call readMINEOScard

  nzone = 1

  do i = 1,nrmod-1
     if(rmod(i).eq.rmod(i+1)) nzone = nzone +1
  enddo

  allocate(vrmin(1:nzone),vrmax(1:nzone),qmuD(1:nzone),qkappaD(1:nzone))
  allocate(rrhoD(1:4,1:nzone),vpvD(1:4,1:nzone),vphD(1:4,1:nzone),vsvD(1:4,1:nzone),vshD(1:4,1:nzone),etaD(1:4,1:nzone))

  do i = 1,nzone
     
     
  
  
  

  

  
