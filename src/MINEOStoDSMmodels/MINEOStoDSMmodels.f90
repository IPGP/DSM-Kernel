program modelConverter
  ! Nobuaki Fuji for DSM Kernel
  !  July 2016, Institut de Physique du Globe de Paris
  
  use parameters
  implicit none
  integer :: i,j,k,kmax
  double precision, allocatable, dimension (:,:) :: a,at,ata,atainv,tmparray
  double precision :: rnormalised,coef
  
  call getarg(1,model1d)
  model1d = trim(model1d)
  call getarg(2,psvfile)
  psvfile = trim(psvfile)

  call readMINEOScard

  nzone = 1

  do i = 1,nrmod-1
     if(rmod(i).eq.rmod(i+1)) nzone = nzone +1
  enddo

  allocate(vrmin(1:nzone),vrmax(1:nzone),irmin(1:nzone),irmax(1:nzone),irlength(1:nzone),qmuD(1:nzone),qkappaD(1:nzone))
  allocate(rrhoD(1:4,1:nzone),vpvD(1:4,1:nzone),vphD(1:4,1:nzone),vsvD(1:4,1:nzone),vshD(1:4,1:nzone),etaD(1:4,1:nzone))

  
  j = 1
  vrmin(j) = rmod(1)
  irmin(j) = 1

  do i = 1,nrmod-1
     if(rmod(i).eq.rmod(i+1)) then
        irmax(j) = i
        vrmax(j) = rmod(i)
        ! NF in DSM, qmu and qkappa are constant through one zone
        qmuD(j) = qm(i)
        qkappaD(j) = qk(i)
        ! NF can of course change this easily though ...
        j = j+1
        irmin(j) = i+1
        vrmin(j) = rmod(i+1)
     endif
  enddo
  j = nzone
  vrmax(j) = rmod(nrmod)
  irmax(j) = nrmod

  do j = 1,nzone
     irlength(j) = irmax(j)-irmin(j)+1
  enddo
  
  rrhoD = 0.d0
  vpvD = 0.d0
  vphD = 0.d0
  vsvD = 0.d0
  vshD = 0.d0
  etaD = 0.d0

  
  do j = 1,nzone
     kmax = min(4,irlength(j))
     
     allocate(a(1:irlength(j),1:kmax))
     allocate(at(1:kmax,1:irlength(j)))
     allocate(ata(1:kmax,1:kmax))
     allocate(atainv(1:kmax,1:kmax))
     allocate(tmparray(1:irlength))
     a = 0.d0
     at = 0.d0
     ata = 0.d0
     atainv = 0.d0
     
     do i = irmin(j),irmax(j)
        coef = 1.d0
        rnormalised = rmod(i)/vrmax(nzone)
        do k = 1,kmax
           a(i,k) = coef
           coef = coef * rnormalised 
        enddo

     enddo
     
     at = transpose(a)

     ata = matmul(at,a)
     
     call inverseLU(kmax,ata,atainv)

     tmparray = 0.d0
     tmparray(1:irlength(j)) = dnm(irmin(j):irmax(j))
     tmparray = matmul(at,tmparray)
     rrhoD(1:kmax,j) = matmul(atainv,tmparray)
     
     
     
     deallocate(a,at,ata,atainv,tmparray)
  enddo

  call writepsvmodel
        
end program modelConverter
  
  
  

  

  
