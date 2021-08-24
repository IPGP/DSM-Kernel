subroutine pinputKernel

  use parameters
  use kernels
  implicit none
  character(120) :: tmpfile
  character(120) :: dummy
  integer :: iitmp,maxlmax
  character(200) :: commandline
  integer, external :: getpid


  
  write(tmpfile,"(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_SynViewer'//tmpfile
  
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(5,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)

  open(unit=1,file=tmpfile,status='unknown')
  ! 0a
  read(1,110) SGTinfo
  SGTinfo = trim(SGTinfo)
  ! 0b
  read(1,110) parentDir
  parentDir = trim(parentDir)
  ! 1a
  read(1,110) eventName
  eventName = trim(eventName)
  ! 1b
  read(1,*) evla, evlo, evdepth
  ! 1c
  read(1,*) mt(1), mt(2), mt(3), mt(4), mt(5), mt(6)
  ! 2a
  read(1,110) stationName
  stationName = trim(stationName)
  ! 2b 
  read(1,*) stla, stlo
  ! 3
  read(1,110) phase
  phase = trim(phase)
  ! 4
  read(1,110) compo
  compo = trim(compo)
  ! 5
  read(1,110) paramWRT
  paramWRT = trim(paramWRT)
  timeincrementV=0.d0
  if((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'betaV').or.&
     (trim(paramWRT).eq.'allV').or.(trim(paramWRT).eq.'video').or.&
     (trim(paramWRT).eq.'vRSGT').or.(trim(paramWRT).eq.'vTSGT')) then
     ! alphaV, betaV, allV are for waveform partials, video for forward modelling (Green's function)
    ! 5-bis
      read(1,*) timeincrementV ! in second 
  endif  
  ! 6a
  read(1,*) ibwfilt
  ! 6b
  read(1,110) dummy
  freqid(0) = trim(dummy)
  if(ibwfilt.eq.1) then
     ! 6c
     read(1,*) fclp(0), fchp(0), npButterworth
     ! in the near future we will calculate multiple frequencies in the same time
  endif
  ! 7
  twin(1:4) = 0.d0
  read(1,*) twin(1),twin(2),twin(3),twin(4)
  ! 8
  read(1,*) itranslat
  ! Aa
  read(1,*) ipdistance
  ! Ab
  read(1,*) c_red_reci
  ! Ba
  read(1,*) ifastFFT
  if(ifastFFT.eq.1) then
     ! Bb
     read(1,*) fmin, fmax
  endif
  ! Ca
  read(1,*) dph, ph1
  ! Cb
  read(1,*) dth, thw
  ! Cd
  read(1,*) rmin,rmax,rdelta
  ! Da
  read(1,*) start, end
  ! Db
  read(1,*) samplingHz
  ! Ea
  read(1,*) calculrapide
  if(calculrapide.ne.0.d0) then
     ! Eb
     read(1,*)  nntype
  else
     nntype = 1
  endif
  allocate(minici(1:nntype,0:nfilter))
  allocate(idecidetype(1:nntype))
  
  if(calculrapide.ne.0.d0) then
     do iitmp = 1, nntype
        ! Ec 
        read(1,*) idecidetype(iitmp)
     enddo
  endif
  read(1,*) iPSVSH
  close(1,status='delete')


  commandline = 'mkdir -p '//trim(parentDir)
  call system(commandline)
  commandline = 'mkdir -p '//trim(parentDir)//'/log'
  call system(commandline)
  commandline = 'mkdir -p '//trim(parentDir)//'/tmp'
  call system(commandline)
  commandline = 'mkdir -p '//trim(parentDir)//'/tmpvideo'
  call system(commandline)
  commandline = 'mkdir -p '//trim(parentDir)//'/tmpvideoparts'
  call system(commandline)      
    commandline = 'mkdir -p '//trim(parentDir)//'/seriousfrechet'
  call system(commandline)   

  call pinputDSM(DSMconfFile,PoutputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch,SGTinfo)
  call readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
 
  write(tmpfile,"(Z5.5)") getpid()
  tmpfile = 'tmpworkingfile_for_psvmodel'//tmpfile
  
  call readpsvmodel(psvmodel,tmpfile)
  INFO_TSGT = trim(parentDir)//"/INFO_TSGT.TXT"
  INFO_RSGT = trim(parentDir)//"/INFO_RSGT.TXT"
  rsampletxt = trim(parentDir)//"/rsample.txt"
  modelcard = trim(parentDir)//"/"//trim(modelname)//".card"

  synnfile = trim(parentDir)//"/"//trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"s.dat"




  if(itranslat.eq.1) then
     call translat (stla,stla)
     call translat (evla,evla)
  endif


  return
end subroutine pinputKernel
     

subroutine pinputDSM(DSMconfFile,outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch,SGTinfo)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_SGTcalcul'
  character(120) :: dummy,outputDir,psvmodel,modelname,DSMconfFile,SGTinfo
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch,tsgtswitch,synnswitch,dummyint
  integer, external :: getpid
  character(120) :: tmpfile

  write(tmpfile, "(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_SGTcalcul'//tmpfile

  open(unit=2, file=SGTinfo)
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
 
  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) DSMconfFile
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  read(1,*) tlen
  read(1,*) rmin_,rmax_,rdelta_
  read(1,*) r0min
  r0max=r0min
  r0delta=20.d0
  read(1,*) thetamin,thetamax,thetadelta
  read(1,*) imin,imax
  read(1,*) rsgtswitch,tsgtswitch,synnswitch,dummyint
  close(1,status='delete')

end subroutine pinputDSM

subroutine readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_DSMconf'
  character(120) :: dummy,DSMconfFile
  real(kind(0d0)) :: re,ratc,ratl,omegai
  integer  :: maxlmax
  integer, external :: getpid
  character(120) :: tmpfile

  write(tmpfile,"(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_DSMconf'//tmpfile


  open(unit=2, file=DSMconfFile, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
  
 
  open(unit=1,file=tmpfile,status='unknown')
  read(1,*) re
  read(1,*) ratc
  read(1,*) ratl
  read(1,*) omegai
  read(1,*) maxlmax
  close(1,status='delete')


end subroutine readDSMconf
  
subroutine readpsvmodel(psvmodel,tmpfile)
  implicit none
  character(120) :: psvmodel, tmpfile, dummy
  open(unit=2, file=psvmodel, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
end subroutine readpsvmodel



subroutine translat(geodetic,geocentric)

  implicit none
  real(kind(0d0)),parameter ::  flattening = 1.d0 / 298.25d0
  real(kind(0d0)), parameter :: pi = 3.1415926535897932d0 
  real(kind(0d0)) :: geocentric, geodetic 
  integer :: flag
  flag = 0
  if(geodetic .gt. 90.d0) then
     geodetic = 1.8d2 - geodetic
     flag = 1
  endif
  
  geodetic = geodetic / 1.8d2 * pi
  geocentric = datan((1.d0-flattening)*(1.d0-flattening)* dtan(geodetic) )
  geocentric = geocentric * 1.8d2 / pi
  
  if(flag .eq. 1) then
     geocentric = 1.8d2 - geocentric
  endif

  return
end subroutine translat


subroutine lsmoothfinder (tlen, np0, freq, lsmooth)

  implicit none
  real(kind(0d0)) :: tlen, freq
  integer :: np0, np, lsmooth, i
  np = 1
  do while (np<np0)
     np = np*2
  enddo
  lsmooth = int(0.5*tlen*freq/dble(np))
  i = 1

  do while (i<lsmooth)
     i = i*2
  enddo

  lsmooth = i

  return 

end subroutine lsmoothfinder

subroutine find_cmb(rcmb,nzone,vrmin,vrmax,vsv,vsh)
  implicit none
  double precision :: rcmb
  integer :: nzone
  double precision :: vrmin(nzone), vrmax(nzone)
  double precision :: vsv(4,nzone),vsh(4,nzone)
  integer :: izone
 
  do izone=1,nzone
     if((vsv(1,izone).eq.0.d0).and.(vsv(2,izone).eq.0.d0).and. &
          (vsv(3,izone).eq.0.d0).and.(vsv(4,izone).eq.0.d0).and. &
          (vsv(1,izone).eq.0.d0).and.(vsv(2,izone).eq.0.d0).and. &
          (vsv(3,izone).eq.0.d0).and.(vsv(4,izone).eq.0.d0)) then
        rcmb=vrmax(izone)
     endif
  enddo
  
end subroutine find_cmb
     

subroutine calstg_for_card(r,nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa,array)

  ! Computing the structure grid points.
  implicit none
  integer:: nzone
  real(kind(0d0)):: r,rrho(4,nzone),vpv(4,nzone),vph(4,nzone),vsv(4,nzone),vsh(4,nzone),eta(4,nzone)
  real(kind(0d0)):: qmu(nzone), qkappa(nzone),vrmin(nzone),vrmax(nzone)
  !real(kind(0d0)), parameter:: rmax  = 6371.d0
  real(kind(0d0)):: rmax
  real(kind(0d0)):: rho,ecKx,ecKy,ecKz
  real(kind(0d0)):: ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  complex(kind(0d0)):: coef1,coef2
  integer:: izone,j
  real(kind(0d0)):: array(1:9)
 
  rmax = vrmax(nzone)
  array = 0.d0
  do izone = 1, nzone
     if((r.gt.vrmin(izone)).and.(r.le.vrmax(izone))) then
        
        coef1 = cmplx(qmu(izone))
        coef2 = cmplx(qkappa(izone))
        trho = 0.d0
        tvpv = 0.d0
        tvph = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        teta = 0.d0
        do j=1,4
           if ( j.eq.1 ) then
              coef = 1.d0
           else
              coef = coef * (r / rmax )
           endif
           trho  = trho  + rrho(j,izone)  * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        rho = trho
        ecL  = rho * tvsv * tvsv
        ecN  = rho * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * ecL )
        !kappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
        ecKx = ecA - 4.d0 / 3.d0 * ecN
        ecKy = ecF + 2.d0 / 3.d0 * ecN
        ecKz = ( ecC + 2.d0 * ecF ) / 3.d0


        array(1) = 1.d3 * r
        array(2) = 1.d3 * rho
        array(3) = 1.d3 * tvpv
        array(4) = 1.d3 * tvsv
        array(5) = qkappa(izone)
        array(6) = qmu(izone)
        array(7) = 1.d3 * tvph
        array(8) = 1.d3 * tvsh
        array(9) = teta
     endif
  enddo
 

  
  return
end subroutine calstg_for_card


subroutine azimth(ellips,slat,slon,rlat,rlon,delta,azim,bazim)

  !   This routine uses Euler angles to find the geocentric distance,
  !   azimuth, and back azimuth for a source-reciever pair.
  !
  !   Input
  !
  !     slat  - source geographic latitude in decimal degrees
  !     slon  - source longitude in decimal degrees
  !     rlat  - reciever geographic latitude in decimal degrees
  !     rlon  - reciever longitude in decimal degrees
  !
  !   Output
  !
  !     delta - geocentric source-reciever distance in decimal degrees of arc
  !     azim  - geocentric azimuth from the source to the reciever
  !     bazim - geocentric back azimuth from the reciever to the source
  !
  !   The distance calculated here delta is always between 0 and 180 degrees. 
  !   Accordingly, the azimuth and back azimuth are defined for the minor 
  !   arc between (slat,slon) and (rlat,rlon).
  !
  !   if ellips = 0 then geocentric = geocentric
  !     because in NF version it is already taken into account so ellips should be 0 always

  
  implicit none
  real(kind(0d0)) :: dtor,e,slatra,slat,w,s,scolat,rlatra,rlat,rcolat,slonra,rlon,c2,s2,c1,s1,slatrc,x0,y0,z0,x1,y1,z1,x2,y2,z2,slon,rlonra,delta,azim,bazim,pi
  real(kind(0d0)), parameter :: flt = 298.25d0
  integer :: ellips
  
  
  dtor=4.d0*datan(1.d0)/180.d0
  pi=4.d0*datan(1.d0)
  
  if(ellips.ne.0) then
     e=1.d0/flt
  else
     e=0.d0
  endif

  
  !   Convert to geocentric coordinates and from latitude to colatitude.

  slatra=dtor*slat
  w=dsin(slatra)
  s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(slatra)
  !scolat=1.5707963d0-slatra+s
  scolat=pi*5.d-1-slatra+s  
  rlatra=dtor*rlat
  w=dsin(rlatra)
  s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(rlatra)
  !rcolat=1.5707963d0-rlatra+s
  rcolat=pi*5.d-1-rlatra+s

  slonra=slon*dtor
  rlonra=rlon*dtor
  c2=dcos(scolat)
  s2=dsin(scolat)
  c1=dcos(slonra)
  s1=dsin(slonra)
  slatrc=dsin(rcolat)



  
  !  Find the azimuth and distance by rotating the source to the north pole.
  
  x0=slatrc*dcos(rlonra)
  y0=slatrc*dsin(rlonra)
  z0=dcos(rcolat)
  x1=c1*x0+s1*y0
  
  z0=dcos(rcolat)
  x1=c1*x0+s1*y0
  y1=-s1*x0+c1*y0
  z1=z0
  x2=c2*x1-s2*z1
  y2=y1
  z2=c2*z1+s2*x1
  call angles(x2,y2,z2,delta,azim)
  azim=180.d0-azim
  
  !  Find the back azimuth by rotating the receiver to the north pole.
  
  c2=dcos(rcolat)
  s2=dsin(rcolat)
  c1=dcos(rlonra)
  s1=dsin(rlonra)
  slatrc=dsin(scolat)
  x0=slatrc*dcos(slonra)
  y0=slatrc*dsin(slonra)
  z0=dcos(scolat)
  x1=c1*x0+s1*y0
  y1=-s1*x0+c1*y0
  z1=z0
  x2=c2*x1-s2*z1
  y2=y1
  z2=c2*z1+s2*x1
  call angles(x2,y2,z2,delta,bazim)
  bazim=180.d0-bazim
  
  return
end subroutine azimth


subroutine angles(x,y,z,theta,phi)

  !   Finds the angles theta and phi of a spherical polar coordinate
  !   system from the cartesion coordinates x, y, and z.
  
  implicit none
  real(kind(0d0)) :: pi,rtod,arg1,x,y,theta,phi,z
  ! real(kind(0d0)), parameter :: eps = 1.d-14
  real(kind(0d0)), parameter :: eps = 0.d0
  

  pi=4.d0*datan(1.d0)
  
  rtod=180.d0/pi
  arg1=dsqrt(x*x+y*y)
  theta=datan2(arg1,z)
  if(dabs(x).le.eps.and.dabs(y).le.eps) then
     phi=0.d0
  else
     phi=datan2(y,x)
  endif
  phi=phi*rtod
  theta=theta*rtod
  
  return
end subroutine angles


subroutine rotmat(thetas0,phis0,thetar0,phir0,azim0)

  use rotate
  
  !   See notes for the definitions and transformations of coordinate systems.
  implicit none

  real(kind(0d0)):: thetas0,phis0,thetar0,phir0,azim0
  real(kind(0d0)) :: thetas,phis,thetar,phir,azim,alpha,pi
  real(kind(0d0)) :: r1(3,3),r2(3,3),r3(3,3)
  real(kind(0d0)) :: rr(3,3),rt(3,3),rtmp(3,3)
  integer :: i,j,k
  
  !common/rotate/rr,rt
  
  rr(1:3,1:3)=cc(1:3,1:3)
  rt(1:3,1:3)=ct(1:3,1:3)

  pi=4.d0*datan(1.d0)
  
  thetas=thetas0
  phis=phis0
  thetar=thetar0
  phir=phir0
  azim=azim0*pi/180.d0
  
  alpha=pi/2.d0-azim
  if(alpha.le.(-pi)) alpha=2.d0*pi+alpha
  
  !   Initializing all matrices.
  
  do i=1,3
     do j=1,3
        r1(i,j)=0.d0
        r2(i,j)=0.d0
        r3(i,j)=0.d0
        rtmp(i,j)=0.d0
        rr(i,j)=0.d0
        rt(i,j)=0.d0
     enddo
  enddo
  
  !   Rotation matrix r4 from geographic (x, y, z) to path defined (x', y', z').
	
  r1(1,1)=dcos(phis)
  r1(1,2)=dsin(phis)
  r1(2,1)=-r1(1,2)
  r1(2,2)=r1(1,1)
  r1(3,3)=1.d0
  r2(1,1)=dsin(thetas)
  r2(1,3)=dcos(thetas)
  r2(2,2)=1.d0
  r2(3,1)=-r2(1,3)
  r2(3,3)=r2(1,1)
  r3(1,1)=1.d0
  r3(2,2)=dcos(alpha)
  r3(2,3)=dsin(alpha)
  r3(3,2)=-r3(2,3)
  r3(3,3)=r3(2,2)
  
  do i=1,3
     do j=1,3
        do k=1,3
           rtmp(i,j)=rtmp(i,j)+r3(i,k)*r2(k,j)
        enddo
     enddo
  enddo
  do i=1,3
     do j=1,3
        do k=1,3
           rr(i,j)=rr(i,j)+rtmp(i,k)*r1(k,j)
        enddo
     enddo
  enddo
        
  !   Setting neglegible elements in rr to zero.
  
  do i=1,3
     do j=1,3
        if(dabs(rr(i,j)).le.1.d-14) rr(i,j)=0.d0
     enddo
  enddo

  !   Getting the transpose of rr.
  
  do i=1,3
     do j=1,3
        rt(i,j)=rr(j,i)
     enddo
  enddo
	
  do i=1,3
     do j=1,3
        rtmp(i,j)=0.d0
        do k=1,3
           rtmp(i,j)=rtmp(i,j)+rr(i,k)*rt(k,j)
        enddo
     enddo
  enddo
	
  cc(1:3,1:3)=rr(1:3,1:3)
  ct(1:3,1:3)=rt(1:3,1:3)


  return
end subroutine rotmat



subroutine hgridsimple(distan,edge,dph0,width,dth0)
  use angles
  implicit none
  real(kind(0d0)) :: dph0,dth0,dph,dth,distan,edge,width 
  real(kind(0d0)) :: midi
  real(kind(0d0)), parameter :: dismin = 0.d0
  integer :: i,j,ip,ith


  

  nphi0=int((distan+2.d0*edge+0.00001)/dph0)+1
  dph=(distan+2.d0*edge)/dble(nphi0-1)
  allocate(phi00(nphi0))
  do i=1,nphi0
     phi00(i)=-edge+dble(i-1)*dph
  enddo
      
  ntheta=1
  if(width.gt.0.d0) ntheta=2*int(width/dth0)+1	
  allocate(theta0(ntheta))
  
  if(width.eq.0.d0) then
     theta0(1)=90.d0
     dth=dph
  else	
     midi=int(ntheta/2)+1
     theta0(midi)=90.d0
     dth=width/dble(midi-1)
     do i=midi-1,1,-1
        theta0(i)=theta0(i+1)-dth
     enddo
     do i=midi+1,ntheta
        theta0(i)=theta0(i-1)+dth
     enddo
  endif
  theta0(ntheta/2+1)=90.d0
      
  
      
  !     Delete phi samples within dismin from source or receiver. 
  
  nphi=0
  do i=1,nphi0
     if((dabs(phi00(i)).lt.dismin).or.(dabs(phi00(i)-distan).lt.dismin)) then
     else
        nphi=nphi+1
     endif
  end do
  allocate(phi0(nphi))
  j=0
  do i=1,nphi0
     if((dabs(phi00(i)).lt.dismin).or.(dabs(phi00(i)-distan).lt.dismin)) then
     else
        j=j+1
        phi0(j)=phi00(i)
     endif
  enddo
  
      
  allocate(phitheta(nphi,ntheta))
  allocate(thetaphi(nphi,ntheta))
  
  do ip=1,nphi
     do ith=1,ntheta
        phitheta(ip,ith)= phi0(ip)
        thetaphi(ip,ith)= theta0(ith)
     enddo
  enddo
  return
end subroutine hgridsimple


subroutine convertPath2Geo
  use rotate
  use angles
  implicit none
  real(kind(0d0)) :: xp,yp,zp,x,y,z,pi
  !real(kind(0d0)) :: cc(3,3),ct(3,3)
  integer :: ith, ip
  !common/rotate/cc,ct
  
  pi=4.d0*datan(1.d0)
  
  ! phitheta and thetaphi : Path-specific system
  ! phi and theta : Geographic system
  allocate(phi(nphi,ntheta),theta(nphi,ntheta))
  do ith=1,ntheta
     do ip=1,nphi
        xp=dsin(thetaphi(ip,ith)*pi/180.d0)* &
             dcos(phitheta(ip,ith)*pi/180.d0)
        yp=dsin(thetaphi(ip,ith)*pi/180.d0)* &
             dsin(phitheta(ip,ith)*pi/180.d0)
        zp=dcos(thetaphi(ip,ith)*pi/180.d0) 
        x=ct(1,1)*xp+ct(1,2)*yp+ct(1,3)*zp
        y=ct(2,1)*xp+ct(2,2)*yp+ct(2,3)*zp
        z=ct(3,1)*xp+ct(3,2)*yp+ct(3,3)*zp
        theta(ip,ith)=datan(dsqrt(x**2+y**2)/z)
        if(z.lt.0.d0) theta(ip,ith)=theta(ip,ith)+pi
        if((x.gt.0.d0).and.(y.eq.0.d0)) then
           phi(ip,ith)=0.d0
        elseif((x.gt.0.d0).and.(y.gt.0.d0)) then
           phi(ip,ith)=datan(y/x)
        elseif((x.eq.0.d0).and.(y.gt.0.d0)) then
           phi(ip,ith)=pi/2.d0
        elseif((x.lt.0.d0).and.(y.gt.0.d0)) then
           phi(ip,ith)=pi-datan(y/(-x))
        elseif((x.lt.0.d0).and.(y.eq.0.d0)) then
           phi(ip,ith)=pi
        elseif((x.lt.0.d0).and.(y.lt.0.d0)) then
           phi(ip,ith)=pi+datan(y/x)
        elseif((x.eq.0.d0).and.(y.lt.0.d0)) then
           phi(ip,ith)=3.d0*pi/2.d0
        else
           phi(ip,ith)=2.d0*pi-datan(-y/x)
        endif
     enddo
  enddo
  return
end subroutine convertPath2Geo


subroutine calculateSineCosine(azim)
  
  !   Compute the distances, azimuths and back-azimuths needed in calculating 
  !   the kernels. Note: these calculations are done in the path specific 
  !   coordinate system in which the source and receiver are both at zero
  !   latitude and 0 and distan longitudes, respectively.

  use angles
  implicit none
  real(kind(0d0)) :: pi,xlat,xlon,phisq,phiqs,phirq,phirs
  integer :: ith,ip
  real(kind(0d0)) :: distanx,azimr,bazimr,azims,bazims
  real(kind(0d0)) :: azim
  
  pi=4.d0*datan(1.d0) 
  do ith=1,ntheta
     do ip=1,nphi
        xlat=90.d0-theta(ip,ith)*180.d0/pi
        xlon=phi(ip,ith)*180.d0/pi
        call azimth(0,xlat,xlon,rlat,rlon,distanx,azimr,bazimr)
        deltar(ip,ith)=distanx
        call azimth(0,slat,slon,xlat,xlon,distanx,azims,bazims)
        deltas(ip,ith)=distanx
        phirq=pi-azimr*pi/180.d0
        if(phirq.lt.0.d0) phirq=phirq+2.d0*pi
        if(phirq.gt.(2.d0*pi)) phirq=phirq-2.d0*pi
        phiqs=pi-azims*pi/180.d0
        if(phiqs.lt.0.d0) phiqs=phiqs+2.d0*pi
        if(phiqs.gt.(2.d0*pi)) phiqs=phiqs-2.d0*pi
        phisq=pi-bazims*pi/180.d0
        if(phisq.lt.0.d0) phisq=phisq+2.d0*pi
        if(phisq.gt.(2.d0*pi)) phisq=phisq-2.d0*pi
        crq(ip,ith)=dcos(phirq)
        srq(ip,ith)=dsin(phirq)
        crq2(ip,ith)=dcos(2.d0*phirq)
        srq2(ip,ith)=dsin(2.d0*phirq)
        csq(ip,ith)=dcos(phisq)
        ssq(ip,ith)=dsin(phisq)
        csq2(ip,ith)=dcos(2.d0*phisq)
        ssq2(ip,ith)=dsin(2.d0*phisq)
         cqs(ip,ith)=dcos(phiqs)
        sqs(ip,ith)=dsin(phiqs)
        cqs2(ip,ith)=dcos(2.d0*phiqs)
        sqs2(ip,ith)=dsin(2.d0*phiqs)
     enddo
  enddo

  ! for computing the reference seismogram. 
  
  phirs=pi-azim*pi/180.d0
  crq(0,0)=dcos(phirs)
  srq(0,0)=dsin(phirs)
  crq2(0,0)=dcos(2.d0*phirs)
  srq2(0,0)=dsin(2.d0*phirs)

  return
end subroutine calculateSineCosine



subroutine fwinDeterminator
  use parameters
  use tmpSGTs
  implicit none
  real(kind(0d0)) :: dt,tpw1,tpw2,xfwin
  integer :: istart,iend,ift,i
  
  istart = iWindowStart
  iend = iWindowEnd

  dt = dtn
  tpw1 = twin(1)
  tpw2 = twin(4)
  
  if(twin(1).lt.t(istart)) twin(1)=t(istart)
  if(twin(4).gt.t(iend)) twin(4)=t(iend)
  itwin=int(twin/dt)

  
  do ift = 0, nfilter
     do i = istart,iend
        if((i.le.itwin(1)).or.(i.ge.itwin(4))) then
           fwin(ift,i)=0.d0
        elseif((i.gt.itwin(1)).and.(i.lt.itwin(2))) then
           xfwin=dsin(0.5d0*pi*dble(i-itwin(1))/dble(itwin(2)-itwin(1)))
           fwin(ift,i)=xfwin*xfwin
        elseif((i.ge.itwin(2)).and.(i.le.itwin(3))) then
           fwin(ift,i)=1.d0
        else
           xfwin=dsin(0.5d0*pi*dble(i-itwin(4))/dble(itwin(3)-itwin(4)))
           fwin(ift,i)=xfwin*xfwin
        endif
     enddo
     nt1(ift) = itwin(1)
     nt2(ift) = itwin(4)
  enddo
  
  return
end subroutine fwinDeterminator

subroutine coeffCalculator
  use parameters
  use tmpSGTs
  implicit none
  
  ! Here I replace those expressions for Q with Fuji et al. 2010
  !
  !    I replace coeff(8) with J^lnA
  !      2013.5.6.
  !    I define jacobianFuji(ir,ifreq) 2013.5.28.


  integer :: ift,it,ir,ifreq
  real(kind(0d0)) :: f0
  real(kind(0d0)), parameter :: gnormt = 1.d0, gnorma = 1.d0
  
  denomv=0.d0
  denomu=0.d0
  do ift=0,nfilter	
     do it=nt1(ift),nt2(ift)
        if((it.eq.nt1(ift)).or.(it.eq.nt2(ift))) then
           denomv(ift)=denomv(ift)+0.5d0*v0(ift,it)**2*dtn
           denomu(ift)=denomu(ift)+0.5d0*u0(ift,it)**2*dtn
        else
           denomv(ift)=denomv(ift)+v0(ift,it)**2*dtn
           denomu(ift)=denomu(ift)+u0(ift,it)**2*dtn
        endif
     enddo
  enddo

  
  allocate(coeff(0:nfilter,1:nktype,1:nr,iWindowStart:iWindowEnd))
  allocate(coeffV(1:nkvtype,1:nr))
  allocate(jacobianFuji(1:nr,fmin:fmax))

  coeff=0.d0  
  coeffV=0.d0

  do ir=1,nr
     coeffV(1,ir) = -2.d0*gnormt*rhom(ir)*vpm(ir)**2
     coeffV(2,ir) = -4.d0*gnormt*rhom(ir)*vsm(ir)**2        
  enddo

  do ir=1,nr
     do ift=0,nfilter
        f0=(fclp(ift)+fchp(ift))/2.d0
        do it=nt1(ift),nt2(ift)
           coeff(ift,1,ir,it)=-2.d0*gnormt*rhom(ir)*vpm(ir)**2*&
                v0(ift,it)*dtn/denomv(ift)
           
           coeff(ift,3,ir,it)= 2.d0*gnorma*rhom(ir)*vpm(ir)**2*&
                u0(ift,it)*dtn/denomu(ift)
           coeff(ift,4,ir,it)=-2.d0*gnorma*rhom(ir)* &
                (vpm(ir)**2-4.d0*vsm(ir)**2/3.d0)* &
     	        (2.d0*dlog(f0)*u0(ift,it)/pi-hu0(ift,it))*dtn/ &
     	        (qkp(ir)*denomu(ift))
           coeff(ift,5,ir,it)=-4.d0*gnormt*rhom(ir)*vsm(ir)**2* &
     	        v0(ift,it)*dtn/denomv(ift)
           coeff(ift,7,ir,it)=4.d0*gnorma*rhom(ir)*vsm(ir)**2* &
                u0(ift,it)*dtn/denomu(ift)

           ! replace the older version for Q
           !if(qmm(ir).ne.0.) coeff(ift,8,ir,it)=2.d0*gnorma*rhom(ir)* &
     	   !     vsm(ir)**2*(2.d0*dlog(f0)*u0(ift,it)/pi-hu0(ift,it))* &
     	   !     dtn/(qmm(ir)*denomu(ift))


           ! This is the J^lnAn *q
           
           coeff(ift,8,ir,it)=gnorma/qmm(ir)* &
                u0(ift,it)*dtn/denomu(ift)

           !print *, coeff(ift,7,ir,it),coeff(ift,8,ir,it)

           ! This version is isotropic

           !coeff(ift,9,ir,it)=-2.d0*gnormt*rhom(ir)*vpm(ir)**2* &
     	   !     v0(ift,it)*dtn/denomv(ift)
           !coeff(ift,10,ir,it)=2.d0*gnorma*rhom(ir)*vpm(ir)**2* &
     	   !     u0(ift,it)*dtn/denomu(ift)
           !coeff(ift,11,ir,it)=-gnormt*rhom(ir)*vpm(ir)**2* &
     	   !     v0(ift,it)*dtn/denomv(ift)
           !coeff(ift,12,ir,it)=gnorma*rhom(ir)*vpm(ir)**2* &
     	   !     u0(ift,it)*dtn/denomu(ift)
           !coeff(ift,13,ir,it)=-4.d0*gnormt*rhom(ir)*vsm(ir)**2* &
     	   !     v0(ift,it)*dtn/denomv(ift)
           !coeff(ift,14,ir,it)=4.d0*gnorma*rhom(ir)*vsm(ir)**2* &
     	   !     u0(ift,it)*dtn/denomu(ift)
           ! denomu for anisotropy in aniso version I will include this, too
           !coeff(ift,15,ir,it)=-4.d0*gnormt*rhom(ir)*vsm(ir)**2* &
           !     v0(ift,2,it)*dtn/denomv(ift,2)
           !coeff(ift,16,ir,it)=-4.d0*gnormt*rhom(ir)*vpm(ir)**2* &
     	   !     v0(ift,2,it)*dtn/denomv(ift,2)
           !coeff(ift,17,ir,it)=-2.d0*gnormt*rhom(ir)*vpm(ir)**2* &
           !     v0(ift,2,it)*dtn/denomv(ift,2)
           !coeff(ift,18,ir,it)=-8.d0*gnormt*rhom(ir)*vsm(ir)**2* &               
	   !     v0(ift,2,it)*dtn/denomv(ift,2)
           !coeff(ift,19,ir,it)=-2.d0*gnormt*rhom(ir)*vsm(ir)**2* &
     	   !     v0(ift,2,it)*dtn/denomv(ift,2)
           !coeff(ift,20,ir,it)=-4.d0*gnormt*rhom(ir)*vsm(ir)**2* &
           !     v0(ift,2,it)*dtn/denomv(ift,2)
        enddo
     enddo
  enddo

  ! jacobianFuji (ir,ifreq)

  do ir=1,nr
     if(qmm(ir).ne.0.d0) then

        do ifreq=fmin,fmax
            if(ifreq.ne.0) then
           !jacobianFuji(ir,ifreq)=(rhom(ir)*vsm(ir)**2*cmplx(2.d0*log(omega(ifreq)/pi),(1+4.d0*log(omega(ifreq)/pi)/pi/qmm(ir)))) &
            !    / (1+2.d0*log(omega(ifreq)/pi/qmm(ir)))/cmplx(1.d0,1.d0/qmm(ir))
             jacobianFuji(ir,ifreq)=(rhom(ir)*vsm(ir)**2)*cmplx(2.d0*log(omega(ifreq)/2.d0/pi)/pi,1.d0+4.d0*log(omega(ifreq)/2.d0/pi)/pi/qmm(ir))/ &
                     cmplx(1.d0+2*log(omega(ifreq)/2.d0/pi)/pi/qmm(ir),0.d0)/cmplx(1.d0,1.d0/qmm(ir))
            else
           endif
       enddo
     endif
  enddo
                       
        


  return
end subroutine coeffCalculator
