subroutine pinput(DSMconfFile,outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch)
  implicit none
  character(120), parameter :: tmpfile='tmpworkingfile_for_SGTforPinv'
  character(120) :: dummy,outputDir,psvmodel,modelname,DSMconfFile
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch,tsgtswitch,synnswitch
  character(120) :: commandline

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
  read(1,*) rsgtswitch,tsgtswitch,synnswitch
  close(1)
  
  ! making directories

  commandline = 'mkdir '//trim(outputDir)
  call system(commandline)
  commandline = 'mkdir '//trim(outputDir)//'/RSGT'
  call system(commandline)
  commandline = 'mkdir '//trim(outputDir)//'/TSGT'
  call system(commandline)
  commandline = 'mkdir '//trim(outputDir)//'/log'
  call system(commandline)  


end subroutine pinput

subroutine readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
  implicit none
  character(120), parameter :: tmpfile='tmpworkingfile_for_DSMconf'
  character(120) :: dummy,DSMconfFile
  real(kind(0d0)) :: re,ratc,ratl,omegai
  integer  :: maxlmax

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
  close(1)


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

!
subroutine udertorsgtPinv(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:2)

  ! icomp = 1
  ! This is for P inversion
  
  rsgt(icomp) = rsgt(icomp) + uder(1,1) + uder(2,2) + uder(3,3)
  
  return
end subroutine udertorsgtPinv

!

subroutine udertotsgtPinv(imt,uder,tsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:4)
  
  ! 1 <= imt <= 4

  ! This is for P inversion

  tsgt(imt) = tsgt(imt) + uder(1,1) + uder(2,2) + uder(3,3)
  return
end subroutine udertotsgtPinv


subroutine utosynn(imt,u,synn)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: u(1:3), synn(1:10)
  
  if(imt.eq.1) then
     synn(1) = synn(1) + u(1)
     synn(5) = synn(5) + u(2)
   endif

   if(imt.eq.2) then
      synn(3) = synn(3) + 5.d-1*u(1)
      synn(4) = synn(4) - 5.d-1*u(1)
      synn(6) = synn(6) - 5.d-1*u(2)
      synn(7) = synn(7) - 5.d-1*u(2)
   endif
   
   if(imt.eq.3) then
      synn(3) = synn(3) - 5.d-1*u(1)
      synn(4) = synn(4) - 5.d-1*u(1)
      synn(6) = synn(6) - 5.d-1*u(2)
      synn(7) = synn(7) + 5.d-1*u(2)
   endif

   if(imt.eq.4) then
      synn(2) = synn(2) + u(1)
      synn(9) = synn(9) + u(2)
   endif

   if(imt.eq.5) then
      synn(10)= synn(10)+ u(3)
   endif

   if(imt.eq.6) then
      synn(8) = synn(8) - u(3)
   endif
   return
 end subroutine utosynn

  



subroutine udertorsgt(icomp,uder,rsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:10)

  if(icomp.eq.1) then ! vertical component
     rsgt(1) = rsgt(1) + uder(1,1)
     rsgt(2) = rsgt(2) - 5.d-1*uder(2,1) - 5.d-1*uder(1,2)
     rsgt(3) = rsgt(3) + 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(4) = rsgt(4) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  endif

  if(icomp.eq.2) then ! radial component
     rsgt(5) = rsgt(5) + uder(1,1)
     rsgt(6) = rsgt(6) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(7) = rsgt(7) - 5.d-1*uder(2,2) + 5.d-1*uder(3,3)
     rsgt(9) = rsgt(9) + 5.d-1*uder(1,2) + 5.d-1*uder(2,1) 
  endif

  if(icomp.eq.3) then ! transverse component
     rsgt(8) = rsgt(8) - 5.d-1*uder(2,3) - 5.d-1*uder(3,2)
     rsgt(10) = rsgt(10) + 5.d-1*uder(1,3) + 5.d-1*uder(3,1)
  endif
end subroutine udertorsgt

!

subroutine udertotsgt(imt,uder,tsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:20)
  if(imt.eq.1) then ! rr source
     tsgt(1) = tsgt(1) + uder(1,1)
     tsgt(3) = tsgt(3) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     tsgt(6) = tsgt(6) - 5.d-1*uder(1,2) - 5.d-1*uder(2,1) 
     tsgt(10)= tsgt(10)+ 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  endif
  if(imt.eq.2) then ! tt source
     tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(8) = tsgt(8) + 2.5d-1*uder(1,2) + 2.5d-1*uder(2,1)
     tsgt(9) = tsgt(9) + 5.d-1*uder(1,1) 
     tsgt(11)= tsgt(11)- 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12)- 2.5d-1*uder(2,2) - 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17)- 1.25d-1*uder(1,2) - 1.25d-1*uder(2,1)
     tsgt(18)= tsgt(18)+ 1.25d-1*uder(1,2) + 1.25d-1*uder(2,1)
     tsgt(19)= tsgt(19)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
  endif
  if(imt.eq.3) then ! pp source
     tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(8) = tsgt(8) - 2.5d-1*uder(1,2) - 2.5d-1*uder(2,1)
     tsgt(9) = tsgt(9) - 5.d-1*uder(1,1)
     tsgt(11)= tsgt(11) -2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12) +2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17) -1.25d-1*uder(1,2) - 1.25d-1*uder(2,1)
     tsgt(18)= tsgt(18) +1.25d-1*uder(1,2) + 1.25d-1*uder(2,1)
     tsgt(19)= tsgt(19) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
  endif
  if(imt.eq.4) then ! rt source
     tsgt(5) = tsgt(5) + uder(1,1)
     tsgt(7) = tsgt(7) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     tsgt(13)= tsgt(13)- 2.5d-1*uder(1,2) - 2.5d-1*uder(2,1)
     tsgt(14)= tsgt(14)+ 2.5d-1*uder(1,2) + 2.5d-1*uder(2,1)
     tsgt(15)= tsgt(15)+ 2.5d-1*uder(2,2) - 2.5d-1*uder(3,3)
     tsgt(16)= tsgt(16)- 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
  endif
  if(imt.eq.5) then ! rp source
     tsgt(13)= tsgt(13)+ 2.5d-1*uder(1,3) + 2.5d-1*uder(3,1)
     tsgt(14)= tsgt(14)+ 2.5d-1*uder(1,3) + 2.5d-1*uder(3,1)
     tsgt(15)= tsgt(15)- 2.5d-1*uder(2,3) - 2.5d-1*uder(3,2)
     tsgt(16)= tsgt(16)- 2.5d-1*uder(2,3) - 2.5d-1*uder(3,2)
  endif
  if(imt.eq.6) then ! tp source
     tsgt(17)= tsgt(17)+ 2.5d-1*uder(1,3) + 2.5d-1*uder(3,1)
     tsgt(18)= tsgt(18)+ 2.5d-1*uder(1,3) + 2.5d-1*uder(3,1)
     tsgt(19)= tsgt(19)- 2.5d-1*uder(2,3) - 2.5d-1*uder(3,2)
     tsgt(20)= tsgt(20)+ 2.5d-1*uder(2,3) + 2.5d-1*uder(3,2)
  endif
end subroutine udertotsgt




!

subroutine setmt(imt,mt)
  implicit none
  ! We use the  rr, tt, pp, rt, rp, tp order here!!

  real(kind(0d0)) :: mt(3,3)
  integer :: imt
  mt = 0.d0
  if(imt.eq.1) mt(1,1) = 1.d0
  if(imt.eq.2) mt(2,2) = 1.d0
  if(imt.eq.3) mt(3,3) = 1.d0
  if(imt.eq.4) mt(1,2) = 1.d0
  if(imt.eq.5) mt(1,3) = 1.d0
  if(imt.eq.6) mt(2,3) = 1.d0
  return
end subroutine setmt

!

subroutine locallyCartesianDerivatives (u,udr,udt,udp,uder,r,theta)
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  complex(kind(0d0)):: u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)
  real(kind(0d0)) :: r,theta
  real(kind(0d0)) :: thetasin,thetacot

  if((theta.eq.0.d0).or.(theta.eq.pi)) then
     uder(1,1) = udr(1)
     uder(2,2) = u(1)/cmplx(r)
     uder(3,3) = u(1)/cmplx(r)     
     return
  endif



  thetasin = sin(theta)
  thetacot = cos(theta)/thetasin

  ! 1,2,3: r,theta,phi; , denotes the partial derivatives

  uder(1,1) = udr(1)
  uder(1,2) = (udt(1)-u(2))/cmplx(r)
  uder(1,3) = (udp(1)/cmplx(thetasin)-u(3))/cmplx(r)

  uder(2,1) = udr(2)
  uder(2,2) = (udt(2)+u(1))/cmplx(r)
  uder(2,3) = (udp(2)/cmplx(thetasin)-u(3)*cmplx(thetacot))/cmplx(r)

  uder(3,1) = udr(3)
  uder(3,2) = udt(3)/cmplx(r)
  uder(3,3) = (udp(3)/cmplx(thetasin)+u(1)+u(2)*cmplx(thetacot))/cmplx(r)

  return
end subroutine locallyCartesianDerivatives

!

subroutine calnl( nzone,vs,iphase,nsl,nll )

  ! counting of nsl and nll.
  implicit none
  integer:: nzone,iphase(*),nsl,nll
  real(kind(0d0)):: vs(4,*)
  integer:: i
  
  nsl = 0
  nll = 0
  do i=1,nzone
     if ( ( vs(1,i).eq.0.d0 ).and.( vs(2,i).eq.0.d0 ).and.( vs(3,i).eq.0.d0 ).and.( vs(4,i).eq.0.d0 ) ) then
        nll = nll + 1
        iphase(i) = 2
     else
        nsl = nsl + 1
        iphase(i) = 1
     endif
  enddo
  return
end subroutine calnl


!


subroutine calgrid( nzone,vrmin,vrmax,vp,vs,rmin,rmax, imax,lmin,tlen,vmin,gridpar,dzpar )
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0 
  integer:: nzone,imax,lmin
  real(kind(0d0)):: vrmin(*),vrmax(*),vp(4,*),vs(4,*)
  real(kind(0d0)):: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
  integer:: izone,i,j
  real(kind(0d0)):: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp

  do izone=1,nzone
     ! computing the S-velocity at each zone
     if ( vs(1,izone).eq.0.d0 ) then
        do i=1,4
           v(i) = vp(i,izone)
        enddo
     else
        do i=1,4
           v(i) = vs(i,izone)
        enddo
     endif
     vs1 = 0.d0
     vs2 = 0.d0
     do j=1,4
        if ( j.eq.1 ) then
           coef1 = 1.d0
        else
           coef1 = coef1 * ( vrmin(izone) / rmax )
        endif
        if ( j.eq.1 ) then
           coef2 = 1.d0
        else
           coef2 = coef2 * ( vrmax(izone) / rmax )
        endif
        vs1 = vs1 + v(j) * coef1
        vs2 = vs2 + v(j) * coef2
     enddo
     rh = vrmax(izone) - vrmin(izone)
     ! computing omega,amax
     omega = 2.d0 * pi * dble(imax) / tlen
     if ( vs1.ge.vs2 ) then
        vmin(izone) = vs2
     else
        vmin(izone) = vs1
     endif
     amax = vrmax(izone)
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )   - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) ) / ( amax * amax )
     if ( gtmp.gt.0.d0 ) then
        dzpar(izone)   = dsqrt( 1.d0/gtmp )
        gridpar(izone) = rh / dzpar(izone)
     else
        dzpar(izone)   = 0.d0
        gridpar(izone) = 0.d0
     endif
  enddo
  ! rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     if ( gridpar(izone).gt.0.d0 ) then
        gridpar(izone) = gridpar(izone) / gtmp
     else
        rh = vrmax(izone) - vrmin(izone)
        gridpar(izone) = rh / ( rmax - rmin ) * 0.1d0
     endif
  enddo
! re-rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     gridpar(izone) = gridpar(izone) / gtmp
  enddo
  return
end subroutine calgrid


!

subroutine calsp( maxnzone,ndc,nsl,nll,iphase,nlayer,nslay,nllay,isp,jsp,ksp,issp,ilsp,lsp,jssp,isdr,jsdr,ildr,jdr,kdr )
  ! Computing the stack points.
  implicit none
  integer:: maxnzone
  integer:: ndc,nsl,nll,iphase(*),nlayer(maxnzone)
  integer:: nslay,nllay
  integer:: isp(maxnzone),jsp(maxnzone),ksp(maxnzone)
  integer:: issp(maxnzone),ilsp(maxnzone)
  integer:: lsp(maxnzone),jssp(maxnzone)
  integer:: isdr,jsdr,ildr,jdr,kdr
  integer:: i,isl,ill

  ! Initialization of the data
  do i=1,maxnzone
     isp(i)  = 0
     jsp(i)  = 0
     ksp(i)  = 0
     issp(i) = 0
     ilsp(i) = 0
     lsp(i)  = 0
     jssp(i) = 0
  enddo
  isdr = 0
  jsdr = 0
  ildr = 0
  jdr = 0
  kdr = 0
  ! computation of isp,jsp,ksp,issp,ilsp,lsp
  isp(1)  = 1
  jsp(1)  = 1
  ksp(1)  = 1
  issp(1) = 1
  ilsp(1) = 1
  lsp(1)  = 1
  jssp(1) = 1
  isl = 0
  ill = 0
  do i=1,ndc
     isp(i+1) = isp(i) + nlayer(i)
     if ( iphase(i).eq.1 ) then
        jsp(i+1) = jsp(i) + 16 * nlayer(i)
        ksp(i+1) = ksp(i) + 2 * ( nlayer(i) + 1 )
        lsp(i+1) = lsp(i) + 4 * nlayer(i)
        isl = isl + 1
        if ( isl.ne.nsl ) then
           issp(isl+1) = issp(isl) + 4 * nlayer(i)
           jssp(isl+1) = jssp(isl) + nlayer(i) + 1
        endif
     else
        jsp(i+1) = jsp(i) + 4 * nlayer(i)
        ksp(i+1) = ksp(i) + ( nlayer(i) + 1 )
        lsp(i+1) = lsp(i) + 2 * nlayer(i)
        ill = ill + 1
        if ( ill.ne.nll ) ilsp(ill+1) = ilsp(ill) + 4 * nlayer(i)
     endif
  enddo
  isdr = 0
  jsdr = 0
  ildr = 0
  jdr = 0
  isdr = isdr + issp(nsl)-1 + 4 * nlayer(ndc+1)
  jsdr = jsdr + jssp(nsl)-1 + nlayer(ndc+1) + 1
  ildr = ildr + 4 * nllay
  jdr =  jdr  + jsp(ndc+1)-1 + 16 * nlayer(ndc+1)
  kdr =  kdr + ksp(ndc+1)-1 + 2 * ( nlayer(ndc+1)+1 )
  
  return
end subroutine calsp

!


subroutine calspo( maxnlay,maxnzone,ndc,rdc,iphase,inlayer,r0,rmin,rmax,ra,isp,spo,spn )
  ! Computing the source location.
  implicit none
  integer:: maxnlay,maxnzone,ndc,iphase(*)
  integer:: inlayer,isp(maxnzone),spn
  real(kind(0d0)):: rdc(*),r0,rmin,rmax,ra(maxnlay+maxnzone+1),spo
  integer:: itmp,idr

  ! checking the parameter
  if ( (r0.lt.rmin).or.(r0.gt.rmax) ) pause 'The source location is improper.(calspo)'
  spo = 0
  ! computing 'spo'
  if ( r0.eq.rmax ) then
     spo = dble( inlayer ) - 0.01d0
     r0 = ra(inlayer) + (spo-dble(inlayer-1)) * ( ra(inlayer+1) -ra(inlayer) )
  else
     itmp = 2
110  continue
     if ( r0.lt.ra(itmp) ) then
        continue
     else
        itmp = itmp + 1
        goto 110
     endif
     spo = dble(itmp-2)  + ( r0-ra(itmp-1) )   / ( ra(itmp)-ra(itmp-1) )
     ! temporal handling
     if ( (spo-dble(itmp-2)).lt.0.01d0 ) then
        spo = dble(itmp-2) + 0.01d0
        r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
     endif
     if ( (spo-dble(itmp-2)).gt.0.99d0 ) then
        spo = dble(itmp-2) + 0.99d0
        r0 = ra(itmp-1) + (spo-dble(itmp-2)) * ( ra(itmp)-ra(itmp-1) )
     endif
  endif
  ! computing 'spn'
  spn = 0
  itmp = 1
 130	continue
  if ( iphase(itmp).eq.1 ) then
     spn = spn + 1
     if ( r0.le.rdc(itmp) ) then
        continue
     else
        itmp = itmp + 1
        goto 130
     endif
  else
     spn = spn + 1
     if ( r0.le.rdc(itmp) ) pause 'The source is in the liquid layer.(calspo)'
     itmp = itmp + 1
     goto 130
  endif
  ! changing 'spo'
  spo = spo - dble( isp(spn) - 1 )
  
  return
end subroutine calspo


!


subroutine calstg( maxnlay,maxnzone,nzone,iphase,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,vnp,vra,rho,kappa,ecKx,ecKy,ecKz, mu,ecL,ecN,r0,spn,ecC0,ecF0,ecL0 )

  ! Computing the structure grid points.
  implicit none
  integer:: maxnlay,maxnzone,nzone,iphase(*),nnl(*),vnp,spn
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
  real(kind(0d0)):: ra(*),rmax
  real(kind(0d0)):: vra(*),rho(*),kappa(*),ecKx(*),ecKy(*),ecKz(*)
  real(kind(0d0)):: mu(*),ecL(*),ecN(*)
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: r0,ecA0,ecC0,ecF0,ecL0
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  integer:: izone,i,j,itmp,jtmp
     
  ! initializing the data
  call vecinit( maxnlay+2*maxnzone+1,vra )
  call vecinit( maxnlay+2*maxnzone+1,rho )
  call vecinit( maxnlay+2*maxnzone+1,kappa )
  call vecinit( maxnlay+2*maxnzone+1,ecKx )
  call vecinit( maxnlay+2*maxnzone+1,ecKy )
  call vecinit( maxnlay+2*maxnzone+1,ecKz )
  call vecinit( maxnlay+2*maxnzone+1,mu )
  call vecinit( maxnlay+2*maxnzone+1,ecL )
  call vecinit( maxnlay+2*maxnzone+1,ecN )
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        vra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
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
              coef = coef * ( vra(itmp) / rmax )
           endif
           trho  = trho  + rrho(j,izone)  * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        rho(itmp) = trho
        ecL(itmp)  = rho(itmp) * tvsv * tvsv
        ecN(itmp)  = rho(itmp) * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * ecL(itmp) )
        kappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
        ecKx(itmp) = ecA - 4.d0 / 3.d0 * ecN(itmp)
        ecKy(itmp) = ecF + 2.d0 / 3.d0 * ecN(itmp)
        ecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
     enddo
     jtmp = jtmp - 1
  enddo
  vnp = itmp
  
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
        coef = coef * ( r0 / rmax )
     endif
     trho  = trho  + rrho(j,spn)  * coef
     tvpv  = tvpv  + vpv(j,spn)   * coef
     tvph  = tvph  + vph(j,spn)   * coef
     tvsv  = tvsv  + vsv(j,spn)   * coef
     tvsh  = tvsh  + vsh(j,spn)   * coef
     teta  = teta  + eta(j,spn)   * coef
  enddo
  ecL0 = trho * tvsv * tvsv
  ecA0 = trho * tvph * tvph
  ecC0 = trho * tvpv * tvpv
  ecF0 = teta * ( ecA0 - 2.d0 * ecL0 )
  
  return
end subroutine calstg

!


subroutine caltstg( maxnlay,maxnzone,nzone,rrho,vpv,vph,vsv,vsh,eta,nnl,ra,rmax,tvra,tkappa,tecKx,tecKy,tecKz,tmu,tecL,tecN)

  ! Computing the structure grid points.
  
  implicit none
  integer:: maxnlay,maxnzone,nzone,nnl(*)
  real(kind(0d0)):: rrho(4,*),vpv(4,*),vph(4,*),vsv(4,*),vsh(4,*),eta(4,*)
  real(kind(0d0)):: ra(*),rmax
  real(kind(0d0)):: tvra(*),tkappa(*),tmu(*)
  real(kind(0d0)):: tecKx(*),tecKy(*),tecKz(*),tecL(*),tecN(*)
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  real(kind(0d0)):: ecA,ecC,ecF
  integer:: izone,i,j,itmp,jtmp

  call vecinit( maxnlay+2*maxnzone+1,tvra )
  call vecinit( maxnlay+2*maxnzone+1,tkappa )
  call vecinit( maxnlay+2*maxnzone+1,tecKx )
  call vecinit( maxnlay+2*maxnzone+1,tecKy )
  call vecinit( maxnlay+2*maxnzone+1,tecKz )
  call vecinit( maxnlay+2*maxnzone+1,tmu )
  call vecinit( maxnlay+2*maxnzone+1,tecL )
  call vecinit( maxnlay+2*maxnzone+1,tecN )
  ! computing the structure grid points
  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        tvra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
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
              coef = coef * ( tvra(itmp) / rmax )
           endif
           trho = trho + rrho(j,izone) * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        tecL(itmp)  = trho * tvsv * tvsv
        tecN(itmp)  = trho * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * tecL(itmp) )
        tkappa(itmp) = ( 4.d0 * ecA + ecC + 4.d0 * ecF - 4.d0 * tecN(itmp) )/ 9.d0
        tecKx(itmp) = ecA - 4.d0 / 3.d0 * tecN(itmp)
        tecKy(itmp) = ecF + 2.d0 / 3.d0 * tecN(itmp)
        tecKz(itmp) = ( ecC + 2.d0 * ecF ) / 3.d0
     enddo
     jtmp = jtmp - 1
  enddo
  
  return
end subroutine caltstg

!

subroutine calinv(vnp,rho,kappa,rhoinv,kappainv)
  ! Computing the inverse of density and elastic constant.
  implicit none
  integer:: vnp,i
  real(kind(0d0)):: rho(*),kappa(*),rhoinv(*),kappainv(*)
  
  do i=1,vnp
     rhoinv(i)   = 1.d0 / rho(i)
     kappainv(i) = 1.d0 / kappa(i)
  enddo
  
  return
end subroutine calinv

!


subroutine submat( nlayer,ha,hb,h )

  ! Subtracting matrix `hb' from matrix `ha'.
  implicit none
  integer:: nlayer
  real(kind(0d0)):: ha(*),hb(*),h(*)
  integer:: i
  
  do i=1,4*nlayer
     h(i) = ha(i) - hb(i)
  enddo
  return
end subroutine submat


!


subroutine calspdr( maxnzone,nzone,iphase,nlayer,jjdr,kkdr )

  implicit none
  integer:: maxnzone,nzone,iphase(*)
  integer:: nlayer(maxnzone),jjdr(*),kkdr(*)
  integer:: izone
     
  jjdr(1) = 1
  kkdr(1) = 1
  do izone=1,nzone-1
     if ( iphase(izone).eq.1 ) then
        jjdr(izone+1) = jjdr(izone) + 16 * nlayer(izone)
        if ( iphase(izone+1).eq.1 ) then
           kkdr(izone+1) = kkdr(izone) + 2 * nlayer(izone)
        else
           kkdr(izone+1) = kkdr(izone) + 2 * ( nlayer(izone)+1 )
        endif
     else
        jjdr(izone+1) = jjdr(izone) + 4 * nlayer(izone)
        if ( iphase(izone+1).eq.1 ) then
           kkdr(izone+1) = kkdr(izone) + ( nlayer(izone)+1 )
        else
           kkdr(izone+1) = kkdr(izone) + nlayer(izone)
        endif
     endif
  enddo

  return
end subroutine calspdr

!



subroutine calmdr( omega,l,nzone,vrmin,vrmax,vmin,dzpar,rmax,sufzone )
  implicit none
  real(kind(0d0)), parameter:: pi=3.1415926535897932d0 
  integer:: l,nzone,sufzone
  real(kind(0d0)):: omega,vrmin(*),vrmax(*),vmin(*),dzpar(*),rmax
  integer:: izone
  real(kind(0d0)):: gtmp,tdzpar
  sufzone = 0
  do izone=1,nzone
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) )   - ( (dble(l)+0.5d0) * (dble(l)+0.5d0) ) / ( vrmax(izone) * vrmax(izone) )
     if ( gtmp.gt.0.d0 ) then
        tdzpar = sqrt( 1.d0/gtmp )
     else
        if ( vrmax(izone).gt.rmax*(1-2.d0*pi/(dble(l)+0.50)) )  then
           tdzpar = 0.d0
        else
           sufzone = izone
           tdzpar = 0.d0
        endif
     endif
  enddo

  return
end subroutine calmdr

!


subroutine calu0( c0,bvec,u )
  implicit none
  complex(kind(0d0)):: c0,bvec,u

  u = u + c0 * bvec
  
  return
end subroutine calu0

!
     

subroutine calulcd0( c0,c0der,rsta,theta, bvec,bvecdt,bvecdp,ulcd )
  
  implicit none
  complex(kind(0d0)):: c0,c0der,bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
  real(kind(0d0)):: rsta,theta
  complex(kind(0d0)):: u1,uder11,uder12,uder13
  
  u1 = c0 * bvec(1)
  uder11 = c0der * bvec(1)
  uder12 = c0 * bvecdt(1)
  uder13 = c0 * bvecdp(1)
  
  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + uder12 / rsta
  ulcd(3) = ulcd(3) + uder13 / rsta / dsin(theta)
  ulcd(5) = ulcd(5) + u1 / rsta
  ulcd(9) = ulcd(9) + u1 / rsta
   
  return
end subroutine calulcd0

!
     

subroutine calu( c0,lsq,bvec,u )
  implicit none
  real(kind(0d0)):: lsq
  complex(kind(0d0)):: c0(2),bvec(3),u(3)
  
  u(1) = u(1) + c0(1) * bvec(1)
  u(2) = u(2) + c0(2) * bvec(2) / dcmplx(lsq)
  u(3) = u(3) + c0(2) * bvec(3) / dcmplx(lsq)

  return
end subroutine calu


!

subroutine calup(c1,c2,lsq,bvec,u)
  implicit none
  real(kind(0d0)) :: lsq
  complex(kind(0d0)) :: c1,c2, bvec(1:3), u(1:3)

  u(1) = u(1) + c1*bvec(1)
  u(2) = u(2) + c2*bvec(2)/dcmplx(lsq)
  u(3) = u(3) + c2*bvec(3)/dcmplx(lsq)
  
  return
end subroutine calup

!

subroutine calup0(c1,bvec,u)
  implicit none
  complex(kind(0d0)) :: c1,u(1:3),bvec(1:3)

  u(1) = u(1) + c1*bvec(1)
  
  return
end subroutine calup0





!


subroutine calulcd( c0,c0der,lsq,rsta,theta, bvec,bvecdt,bvecdp,ulcd )

  implicit none
  real(kind(0d0)):: lsq,rsta,theta
  complex(kind(0d0)):: c0(2),c0der(2)
  complex(kind(0d0)):: bvec(3),bvecdt(3),bvecdp(3),ulcd(9)
  
  complex(kind(0d0)):: u1,u2,u3
  complex(kind(0d0)):: uder11,uder12,uder13
  complex(kind(0d0)):: uder21,uder22,uder23
  complex(kind(0d0)):: uder31,uder32,uder33

  u1 = c0(1) * bvec(1)
  u2 = c0(2) * bvec(2) / dcmplx(lsq)
  u3 = c0(2) * bvec(3) / dcmplx(lsq)
  ! partial derivatives of u
  uder11 = c0der(1) * bvec(1)
  uder12 = c0(1) * bvecdt(1)
  uder13 = c0(1) * bvecdp(1)
  uder21 = c0der(2) * bvec(2) / dcmplx(lsq)
  uder22 = c0(2) * bvecdt(2) / dcmplx(lsq)
  uder23 = c0(2) * bvecdp(2) / dcmplx(lsq)
  uder31 = c0der(2) * bvec(3) / dcmplx(lsq)
  uder32 = c0(2) * bvecdt(3) / dcmplx(lsq)
  uder33 = c0(2) * bvecdp(3) / dcmplx(lsq)
  ! locally Cartesian derivatives of u
  ulcd(1) = ulcd(1) + uder11
  ulcd(2) = ulcd(2) + ( uder12 - u2 ) / rsta 
  ulcd(3) = ulcd(3) + ( uder13 / dsin(theta) - u3 ) / rsta
  ulcd(4) = ulcd(4) + uder21
  ulcd(5) = ulcd(5) + ( uder22 + u1 ) / rsta
  ulcd(6) = ulcd(6) + ( uder23 - u3 * dcos(theta) )  / rsta / dsin(theta)
  ulcd(7) = ulcd(7) + uder31
  ulcd(8) = ulcd(8) + uder32 / rsta
  ulcd(9) = ulcd(9)  + ( ( uder33 + u2 * dcos(theta) ) / dsin(theta) + u1 ) / rsta
  return
end subroutine calulcd
    
!


subroutine matinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  real(kind(0d0)):: a(n1,*)

  do j=1,n2
     do i=1,n1
        a(i,j) = 0.d0
     enddo
  enddo
  
  return
end subroutine matinit

!


subroutine cmatinit( n1,n2,a )
  implicit none
  integer:: n1,n2,i,j
  complex(kind(0d0)):: a(n1,*)
  
  do j=1,n2
     do i=1,n1
        a(i,j) = dcmplx( 0.d0 )
     enddo
  enddo
  return
end subroutine cmatinit

!
     

subroutine vecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none
  integer:: nn,i
  real(kind(0d0)):: b(*)
     
  do i=1,nn
     b(i) = 0.d0
  enddo
  return
end subroutine vecinit

!


subroutine cvecinit( nn,b )

  ! Filling zero to the vector 'g'.
  implicit none 
  integer:: nn,i
  complex(kind(0d0)):: b(*)
  
  do i=1,nn
     b(i) = dcmplx( 0.d0 )
  enddo
  return
end subroutine cvecinit

!


subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )

  implicit none
  integer:: ncomp,nderiv
  real(kind(0d0)):: rsta,rrsta(3)
  complex(kind(0d0)):: g(3*ncomp),u(ncomp)
  real(kind(0d0)):: dh(3)
  
  integer:: ip(3),ier,i,itmp,icomp
  complex(kind(0d0)):: a(3,3),b(3),wk(3)
  real(kind(0d0)):: eps
  eps = -1.d0
  
  do icomp=1,ncomp
     u(icomp) = dcmplx(0.d0)
  enddo
  do i=1,3
     dh(i) = rrsta(i) - rsta
  enddo
   
  if( (dh(2).eq.0.d0).and.(nderiv.eq.0)) then
     itmp = ncomp + 1
     do icomp=1,ncomp
        u(icomp) = g(itmp)
        itmp = itmp + 1
     enddo
     return
  endif    
  do i=1,3
     a(1,i) = dcmplx( 1.d0 )
     a(2,i) = dcmplx( dh(i) )
     a(3,i) = dcmplx( dh(i) * dh(i) / 2.d0 )
  enddo
  call fillinpb(nderiv,b)
  call glu(a,3,3,b,eps,wk,ip,ier)
    
  
  do icomp=1,ncomp
     do i=1,3
        u(icomp) = u(icomp) + b(i) * g( ncomp * (i-1) + icomp )
     enddo
  enddo
end subroutine interpolate

!


subroutine fillinpb( nderiv,b )

  implicit none
  integer:: nderiv
  complex(kind(0d0)):: b(3)
  
  if( (nderiv.ne.0).and.(nderiv.ne.1).and.(nderiv.ne.2) ) pause 'invalid argument (fillinpb)'
  if(nderiv.eq.0) then
     b(1) = dcmplx( 1.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 0.d0 )
  elseif(nderiv.eq.1) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 1.d0 )
     b(3) = dcmplx( 0.d0 )
  elseif(nderiv.eq.2) then
     b(1) = dcmplx( 0.d0 )
     b(2) = dcmplx( 0.d0 )
     b(3) = dcmplx( 1.d0 )
  endif
  
  return
end subroutine fillinpb


!



subroutine calamp( g,l,lsuf,maxamp,ismall,ratl )
  implicit none
  integer:: ismall,l,lsuf
  real(kind(0d0)):: maxamp,ratl
  complex(kind(0d0)):: g(2)
  real(kind(0d0)):: amp,ampratio
  
  ampratio = 0.d0
  amp = dsqrt( zabs( g(1) )**2 + zabs( g(2) )**2 )
  if ( amp.gt.maxamp ) maxamp = amp
  if ( (amp.ne.0.d0).and.(maxamp.ne.0.d0) ) ampratio = amp / maxamp
  if ( ( ampratio.lt.ratl ).and.( l.ge.lsuf ) ) then
     ismall = ismall + 1
  else
     ismall = 0
  endif

  return
end subroutine calamp

!


subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)
  implicit none
  integer:: nzone,lsuf
  real(kind(0d0)):: omega,vrmax(*),vsv(4,*)
  real(kind(0d0)):: tvs,coef
  integer:: i

  tvs = 0.d0
  do i=1,4
     if(i.eq.1) then
        coef = 1.d0
     else
        coef = coef 
     endif
     tvs = tvs + ( vsv(i,nzone) ) * coef
  enddo
  lsuf = int(omega * vrmax(nzone) / tvs - 0.5d0) + 1
  return
end subroutine callsuf


!

subroutine calra_psv(nlayer,inlayer,jnlayer,jnslay,jnllay,gridpar,dzpar,nzone,vrmin,vrmax,iphase,rmin,rmax,nslay,nllay,nnl,re )
  
  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0 
  integer :: nlayer,inlayer,jnlayer,jnslay,jnllay
  integer :: nzone, iphase(*), nslay, nllay, nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
  integer :: izone,itmp,i,ntmp
  real(kind(0d0)) :: rh,re

  nslay=0
  nllay=0
  inlayer=0

  nnl = 0

  jnlayer=0
  jnslay=0
  jnllay=0


  itmp=1 
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)

     if(dzpar(izone).eq.0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     ! ntmp (see Geller & Takeuchi 1995 6.2)     
     nnl(izone) = ntmp
     if ( nnl(izone).lt.5 ) nnl(izone)=5
     if(iphase(izone).eq.1) nslay = nslay+nnl(izone)
     if(iphase(izone).eq.2) nllay=nllay+nnl(izone)
  enddo

  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  nlayer = inlayer 
  jnlayer = jnlayer + inlayer
  jnslay = jnslay + nslay
  jnllay = jnllay + nllay

  return

end subroutine calra_psv

!

subroutine calra2_psv(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,ra,re,nsta, rsta, rrsta, iista,log_solid_liquid ,rs, cista,iphase,istazone,ciista)

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0 
  integer :: nlayer,inlayer
  integer :: nzone,nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,ra(nlayer+nzone+1)
  integer :: izone,itmp,i,nsta,ista, iista(1:3,1:nsta), cista
  real(kind(0d0)) :: rsta(1:nsta),rrsta(1:3,1:nsta)
  real(kind(0d0)) :: rh,re,rs,ctmp
  real(kind(0d0)) :: chikasa
  logical :: log_solid_liquid(1:nsta)
  integer :: iphase(*),istazone(1:nsta),ciista
  
  chikasa = 0.d0

  ra = 0
  ra(1) = rmin
  ciista = 0
  ctmp = 6371.d0
  
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        itmp = itmp + 1
        ra(itmp) = vrmin(izone) &
             &	+ rh * dble(i) / dble( nnl(izone) )
     enddo
  enddo


  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)
     do i=1,nnl(izone)
        do ista = 1, nsta
           if( (ra(itmp).lt.rsta(ista)) .and.(rsta(ista).le.ra(itmp+1))) then
              if(i.ne.nnl(izone)) then
                 istazone(ista) = izone
                 log_solid_liquid(ista) = .true.
                 if(iphase(istazone(ista)).eq.2) then log_solid_liquid(ista) = .false. !pause 'rsta is in liquid layer (calra2_psv)'         
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)
                     
                 iista(1,ista) = i
                 iista(2,ista) = i + 1
                 iista(3,ista) = i + 2
              else
                 istazone(ista) = izone
                 log_solid_liquid(ista) = .true.
                 if(iphase(istazone(ista)).eq.2) log_solid_liquid(ista) = .false. !pause 'rsta is in liquid layer (calra2_psv)'
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)
                 
                 iista(1,ista) = i - 1
                 iista(2,ista) = i
                 iista(3,ista) = i + 1
              endif
              
              if((abs(rs-rsta(ista)).lt.ctmp).and.(abs(rs-rsta(ista)).ge.chikasa)) then
                 cista = ista
                 ciista = itmp
                 ctmp = abs(rs-rsta(ista))
          
              endif
           endif          
        enddo
        itmp = itmp + 1     
     enddo
  enddo
 
  if(cista.eq.0) cista = 1
end subroutine calra2_psv
!
subroutine calcutd(nzone,nnlayer,nnl,tmpc,rat,nn,iphase,spo,spn, ra,kkdr,kc)
  implicit none
  integer :: nzone,nn,nnlayer,spn,kkdr(1:nzone),kc,iphase(1:nzone),nnl(1:nzone)
  complex(kind(0d0)) :: tmpc(1:nn)
  real(kind(0d0)) :: rat,spo,ra(1:nnlayer+nzone+1)
  integer :: nc
  real(kind(0d0)) :: cU(nn),cV(nn),rc
  real(kind(0d0)) :: maxamp,amp(nn)
  integer :: iz,jz,jj,i,ml(nzone),tzone

  do jj=1,nn
     cU(jj) = 0.d0
     cV(jj) = 0.d0
  enddo
  iz = 2
  jz = 1
  do jj=1,nn
     if(iz.le.nzone) then
        if(jj.eq.kkdr(iz)) then
           if(iphase(iz).ne.iphase(iz-1)) jz = jz - 1
           iz = iz + 1
        endif
     endif
     if(iphase(iz-1).eq.1) then
        if(mod((jj-kkdr(iz-1)),2).eq.1) then ! U
           cU(jz) = cdabs(tmpc(jj))
           jz = jz + 1
        else		! V
        endif
     else ! U in fluid
        cU(jz) = cdabs(tmpc(jj))
        jz = jz + 1
     endif
  enddo

  maxamp = -1.d0
  do i=1,jz-1
     amp(i) = cU(i)
     if(maxamp.lt.amp(i)) maxamp = amp(i)
  enddo
!
  maxamp = maxamp * rat ! threshold value
!
  nc = 1
  do i=1,jz-1
     if(amp(i).gt.maxamp) then
        nc = i
        cycle
     endif
  enddo
  i = 1
  do jj=1,nzone
     i = i + nnl(jj)
     ml(jj) = i
  enddo
  do jj=nzone,1,-1
     if(ml(jj).gt.nc) tzone = jj
  enddo
  rc = ra(nc)
  
  do i=1,jz-1
     if( (ra(i).le.rc).and.(rc.lt.ra(i+1)) ) then
        nc = i
        if(tzone.eq.1) then ! case(tzone is innermost zone)
           if(iphase(tzone).eq.1) kc = 1 + 2 * nc
           if(iphase(tzone).eq.2) kc = 1 + nc
        else 
           if(iphase(tzone).eq.1) then
              kc = kkdr(tzone) + 2 * (nc - ml(tzone-1))
           endif
           if(iphase(tzone).eq.2) then
              kc = kkdr(tzone) + nc - ml(tzone-1)
           endif
        endif
     endif
  enddo
  
  return
end subroutine calcutd
