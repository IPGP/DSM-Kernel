program ModelDrawer
  
  use parameters
  implicit none

  real(kind(0d0)) :: tmparray(1:9),radius
  integer :: iradius,nradius,i
  real(kind(0d0)), parameter :: pas = 1.d0 !  (km)
  character(120) :: dummy


  psvmodel = 'tmpworkingfile_for_psvmodel' 

  
  open(unit=1, file=psvmodel,status='unknown')
100 continue
  read(5,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  
  psvmodel = 'tmpworkingfile_for_psvmodel'
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  close(20)

  allocate(vrminD(1:nzone))
  allocate(vrmaxD(1:nzone))
  allocate(rrhoD(1:4,1:nzone))
  allocate(vpvD(1:4,1:nzone))
  allocate(vphD(1:4,1:nzone))
  allocate(vsvD(1:4,1:nzone))
  allocate(vshD(1:4,1:nzone))
  allocate(etaD(1:4,1:nzone))
  allocate(qmuD(1:nzone))
  allocate(qkappaD(1:nzone))


  psvmodel = 'tmpworkingfile_for_psvmodel' 
   
  open(20, file = psvmodel, status = 'old', action='read', position='rewind')
  read(20,*) nzone
  do i = 1, nzone
     read (20, *) vrminD(i), vrmaxD(i), rrhoD(1,i), rrhoD(2,i), rrhoD(3,i), rrhoD(4,i), vpvD(1,i), vpvD(2,i), vpvD(3,i), vpvD(4,i), vphD(1,i), vphD(2,i), vphD(3,i), vphD(4,i), vsvD(1,i), vsvD(2,i), vsvD(3,i), vsvD(4,i), vshD(1,i), vshD(2,i), vshD(3,i), vshD(4,i), etaD(1,i), etaD(2,i), etaD(3,i), etaD(4,i), qmuD(i), qkappaD(i)
  enddo
  close(20)
  
  nradius = int((vrmaxD(nzone)- vrminD(1))/pas)
 
 

  do iradius = 0, nradius
     radius = vrminD(1)+pas*dble(iradius)
     call calstg_for_card(radius,nzone,vrminD,vrmaxD,rrhoD,vpvD,vphD,vsvD,vshD,etaD,qmuD,qkappaD,tmparray(1:9))
     !print *, 6371.d0-radius,",",tmparray(2)*1.d-3,",",tmparray(3)*1.d-3,",",tmparray(4)*1.d-3
     !print *, 6371.d0-radius,",",tmparray(4)*1.d-3
     write(*,130)radius,tmparray(2)*1.d-3,tmparray(3)*1.d-3,tmparray(4)*1.d-3 
130  format (4(f8.2))
  enddo
end program ModelDrawer

  
