program modelConverter
  ! Nobuaki Fuji for DSM Kernel
  !  July 2016, Institut de Physique du Globe de Paris
  
  use parameters
  implicit none
  integer :: i,j,k,kmax,iirlength
  double precision, allocatable, dimension (:,:) :: a,at,ata,atainv
  double precision, allocatable, dimension (:) :: tmparray,tmparray2
  double precision :: rnormalised,coef
  
  call getarg(1,model1d)
  model1d = trim(model1d)
  call getarg(2,psvmodel)
  psvmodel = trim(psvmodel)

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

  do i = 2,nrmod-1
     if(rmod(i).eq.rmod(i+1)) then
        irmax(j) = i
        vrmax(j) = rmod(i)
        ! NF in DSM, qmu and qkappa are constant through one zone
        qmuD(j) = qm(i)
        qkappaD(j) = qk(i)
        etaD(1,j) = 1.d0
        ! NF can of course change this easily though ...
        j = j+1
        irmin(j) = i+1
        vrmin(j) = rmod(i+1)
     endif
  enddo
  j = nzone
  vrmax(j) = rmod(nrmod)
  irmax(j) = nrmod
  qmuD(j)=qm(nrmod)
  qkappaD(j)=qk(nrmod)
  etaD(1,j)=1.d0

  do j = 1,nzone
     irlength(j) = irmax(j)-irmin(j)+1
  enddo
  
  rrhoD = 0.d0
  vpvD = 0.d0
  vphD = 0.d0
  vsvD = 0.d0
  vshD = 0.d0
  !etaD = 0.d0
  
 
  do j = 1,nzone

     !print *,j,irmin(j),irmax(j),vrmin(j),vrmax(j),irlength(j)
     iirlength = irlength(j)
     kmax = min(4,iirlength)
     allocate(tmparray(1:iirlength))
     allocate(tmparray2(1:kmax))
     allocate(a(1:iirlength,1:kmax))
     allocate(at(1:kmax,1:iirlength))
     allocate(ata(1:kmax,1:kmax))
     allocate(atainv(1:kmax,1:kmax))

     a = 0.d0
     at = 0.d0
     ata = 0.d0
     atainv = 0.d0
     
     do i = irmin(j),irmax(j)
        coef = 1.d0
        rnormalised = rmod(i)/vrmax(nzone)
        do k = 1,kmax
           a(i-irmin(j)+1,k) = coef
           coef = coef * rnormalised 
        enddo

     enddo
     
     at = transpose(a)

     ata = matmul(at,a)
     

     call inverseLU(kmax,ata,atainv)



     ! rho
     tmparray = 0.d0
     tmparray2 = 0.d0
     tmparray(1:iirlength) = dnm(irmin(j):irmax(j))
     tmparray2 = matmul(at,tmparray)
     rrhoD(1:kmax,j) = matmul(atainv,tmparray2)
    
     ! vpv
     tmparray = 0.d0
     tmparray2 =0.d0
     tmparray(1:iirlength) = vpv(irmin(j):irmax(j))
     tmparray2 = matmul(at,tmparray)
     vpvD(1:kmax,j) = matmul(atainv,tmparray2)

     ! vsv
     tmparray = 0.d0
     tmparray2 =0.d0
     tmparray(1:iirlength) = vsv(irmin(j):irmax(j))
     tmparray2 = matmul(at,tmparray)
     vsvD(1:kmax,j) = matmul(atainv,tmparray2)
     
     ! vph
     tmparray = 0.d0
     tmparray2 =0.d0
     tmparray(1:iirlength) = vph(irmin(j):irmax(j))
     tmparray2 = matmul(at,tmparray)
     vphD(1:kmax,j) = matmul(atainv,tmparray2)
    
     ! vsh
     tmparray = 0.d0
     tmparray2 =0.d0
     tmparray(1:iirlength) = vsh(irmin(j):irmax(j))
     tmparray2 = matmul(at,tmparray)
     vshD(1:kmax,j) = matmul(atainv,tmparray2)
     

     
     deallocate(a,at,ata,atainv,tmparray,tmparray2)
  enddo

  call writepsvmodel
        
end program modelConverter
  
  
  

  

  
