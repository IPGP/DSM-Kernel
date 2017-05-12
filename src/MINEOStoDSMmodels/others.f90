subroutine readMINEOScard

  use parameters
  implicit none
  integer :: i
  
  open(1,file=model1d,status='old')
  read(1,'(a)') modelid
  read(1,*) ifanis,tref,ifdeck
  read(1,*) nrmod,nicb,ncmb,nmspl
  
  allocate(rmod(nrmod),dnm(nrmod),vpv(nrmod),vsv(nrmod))
  allocate(qk(nrmod),qm(nrmod),vph(nrmod),vsh(nrmod),eta(nrmod))
  
  do i=1,nrmod
     read(1,*) rmod(i),dnm(i),vpv(i),vsv(i), &
          qk(i),qm(i),vph(i),vsh(i),eta(i)
     rmod(i)=rmod(i)*1.d-3
     dnm(i)=dnm(i)*1.d-3
     vpv(i)=vpv(i)*1.d-3
     vsv(i)=vsv(i)*1.d-3
     vph(i)=vph(i)*1.d-3
     vsh(i)=vsh(i)*1.d-3
  enddo
  close(1)


end subroutine readMINEOScard



  
subroutine writepsvmodel
  use parameters
  implicit none
  integer :: i

  open(20, file = psvmodel, status = 'unknown', action='write')
  write(20,*) nzone
  do i = 1, nzone
     write(20,123) vrmin(i),vrmax(i), rrhoD(1,i), rrhoD(2,i), rrhoD(3,i), rrhoD(4,i)
     write(20,124) vpvD(1,i), vpvD(2,i), vpvD(3,i), vpvD(4,i)
     write(20,124) vphD(1,i), vphD(2,i), vphD(3,i), vphD(4,i)
     write(20,124) vsvD(1,i), vsvD(2,i), vsvD(3,i), vsvD(4,i)
     write(20,124) vshD(1,i), vshD(2,i), vshD(3,i), vshD(4,i)
     write(20,123) etaD(1,i), etaD(2,i), etaD(3,i), etaD(4,i), qmuD(i), qkappaD(i)
     !write (20, *) vrmin(i), vrmax(i), rrhoD(1,i), rrhoD(2,i), rrhoD(3,i), rrhoD(4,i), vpvD(1,i), vpvD(2,i), vpvD(3,i), vpvD(4,i), vphD(1,i), vphD(2,i), vphD(3,i), vphD(4,i), vsvD(1,i), vsvD(2,i), vsvD(3,i), vsvD(4,i), vshD(1,i), vshD(2,i), vshD(3,i), vshD(4,i), etaD(1,i), etaD(2,i), etaD(3,i), etaD(4,i), qmuD(i), qkappaD(i)
  enddo
  close(20)
123 format(6(1x,F10.4))
124 format(4(1x,F10.4))

end subroutine writepsvmodel


