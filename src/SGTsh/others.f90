subroutine pinput(DSMconfFile,outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch)
  implicit none
  !character(120), parameter :: tmpfile='tmpworkingfile_for_SGTforPinv'
  character(120) :: dummy,outputDir,psvmodel,modelname,DSMconfFile
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch,tsgtswitch,synnswitch
  character(200) :: commandline

  character(120) :: tmpfile
  integer, external :: getpid


  write(tmpfile,"(Z5.5)") getpid()
  tmpfile='tmpworkingfile_for_SGTforSinv'//tmpfile

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

  commandline = 'mkdir -p '//trim(outputDir)
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/RSGT'
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/TSGT'
  call system(commandline)
  commandline = 'mkdir -p '//trim(outputDir)//'/log'
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


subroutine utosynnSH(imt,u,ssynn)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: u(1:3), ssynn(1:5),synn(1:10)

  synn = dcmplx(0.d0)
  synn(6) = ssynn(1)
  synn(7) = ssynn(2)
  synn(8) = ssynn(3)
  synn(9) = ssynn(4)
  synn(10)= ssynn(5) 
  


  !if(imt.eq.1) then
     !synn(1) = synn(1) + u(1)
     !synn(5) = synn(5) + u(2)
   !endif

   if(imt.eq.2) then
      !synn(3) = synn(3) + 5.d-1*u(1)
      !synn(4) = synn(4) - 5.d-1*u(1)
      synn(6) = synn(6) - 5.d-1*u(2)
      synn(7) = synn(7) - 5.d-1*u(2)
   endif
   
   if(imt.eq.3) then
      !synn(3) = synn(3) - 5.d-1*u(1)
      !synn(4) = synn(4) - 5.d-1*u(1)
      synn(6) = synn(6) - 5.d-1*u(2)
      synn(7) = synn(7) + 5.d-1*u(2)
   endif

   if(imt.eq.4) then
      !synn(2) = synn(2) + u(1)
      synn(9) = synn(9) + u(2)
   endif

   if(imt.eq.5) then
      synn(10)= synn(10)+ u(3)
   endif

   if(imt.eq.6) then
      synn(8) = synn(8) - u(3)
   endif

    
   ssynn(1) = synn(6)
   ssynn(2) = synn(7)
   ssynn(3) = synn(8)
   ssynn(4) = synn(9)
   ssynn(5) = synn(10) 
   
   return
end subroutine utosynnSH
 

!






subroutine udertorsgtSH(icomp,uder,rrsgt)
  implicit none
  integer :: icomp
  complex(kind(0d0)) :: uder(1:3,1:3), rsgt(1:10),rrsgt(1:5)

  ! h6    : rsgt(1)
  ! h7    : rsgt(2)
  ! h8    : rsgt(3)
  ! h9    : rsgt(4)
  ! h10   : rsgt(5)
  
  rsgt = dcmplx(0.d0)
  rsgt(6) = rrsgt(1)
  rsgt(7) = rrsgt(2)
  rsgt(8) = rrsgt(3)
  rsgt(9) = rrsgt(4)
  rsgt(10)= rrsgt(5) 

  

  
  !if(icomp.eq.1) then ! vertical component
     !rsgt(1) = rsgt(1) + uder(1,1)
     !rsgt(2) = rsgt(2) - 5.d-1*uder(2,1) - 5.d-1*uder(1,2)
     !rsgt(3) = rsgt(3) + 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     !rsgt(4) = rsgt(4) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  !endif

  if(icomp.eq.2) then ! radial component
     !rsgt(5) = rsgt(5) + uder(1,1)
     rsgt(6) = rsgt(6) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     rsgt(7) = rsgt(7) - 5.d-1*uder(2,2) + 5.d-1*uder(3,3)
     rsgt(9) = rsgt(9) + 5.d-1*uder(1,2) + 5.d-1*uder(2,1) 
  endif

  if(icomp.eq.3) then ! transverse component
     rsgt(8) = rsgt(8) - 5.d-1*uder(2,3) - 5.d-1*uder(3,2)
     rsgt(10) = rsgt(10) + 5.d-1*uder(1,3) + 5.d-1*uder(3,1)
  endif
  

  rrsgt(1) = rsgt(6)
  rrsgt(2) = rsgt(7)
  rrsgt(3) = rsgt(8)
  rrsgt(4) = rsgt(9)
  rrsgt(5) = rsgt(10) 

  return
end subroutine udertorsgtSH


!
subroutine udertotsgtSH(imt,uder,ttsgt)
  implicit none
  integer :: imt
  complex(kind(0d0)) :: uder(1:3,1:3),tsgt(1:20),ttsgt(1:10)

  tsgt = dcmplx(0.d0)
  
  tsgt(7)  = ttsgt(1)
  tsgt(12) = ttsgt(2)
  tsgt(13) = ttsgt(3)
  tsgt(14) = ttsgt(4)
  tsgt(15) = ttsgt(5)
  tsgt(16) = ttsgt(6)
  tsgt(17) = ttsgt(7)
  tsgt(18) = ttsgt(8)
  tsgt(19) = ttsgt(9)
  tsgt(20) = ttsgt(10)



  !if(imt.eq.1) then ! rr source
     !tsgt(1) = tsgt(1) + uder(1,1)
     !tsgt(3) = tsgt(3) - 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
     !tsgt(6) = tsgt(6) - 5.d-1*uder(1,2) - 5.d-1*uder(2,1) 
     !tsgt(10)= tsgt(10)+ 5.d-1*uder(2,2) - 5.d-1*uder(3,3)
  !endif
  if(imt.eq.2) then ! tt source
     !tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     !tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     !tsgt(8) = tsgt(8) + 2.5d-1*uder(1,2) + 2.5d-1*uder(2,1)
     !tsgt(9) = tsgt(9) + 5.d-1*uder(1,1) 
     !tsgt(11)= tsgt(11)- 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12)- 2.5d-1*uder(2,2) - 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17)- 1.25d-1*uder(1,2) - 1.25d-1*uder(2,1)
     tsgt(18)= tsgt(18)+ 1.25d-1*uder(1,2) + 1.25d-1*uder(2,1)
     tsgt(19)= tsgt(19)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20)+ 1.25d-1*uder(2,2) - 1.25d-1*uder(3,3)
  endif
  if(imt.eq.3) then ! pp source
     !tsgt(2) = tsgt(2) - 5.d-1*uder(1,1)
     !tsgt(4) = tsgt(4) + 2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     !tsgt(8) = tsgt(8) - 2.5d-1*uder(1,2) - 2.5d-1*uder(2,1)
     !tsgt(9) = tsgt(9) - 5.d-1*uder(1,1)
     !tsgt(11)= tsgt(11) -2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(12)= tsgt(12) +2.5d-1*uder(2,2) + 2.5d-1*uder(3,3)
     tsgt(17)= tsgt(17) -1.25d-1*uder(1,2) - 1.25d-1*uder(2,1)
     tsgt(18)= tsgt(18) +1.25d-1*uder(1,2) + 1.25d-1*uder(2,1)
     tsgt(19)= tsgt(19) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
     tsgt(20)= tsgt(20) -1.25d-1*uder(2,2) + 1.25d-1*uder(3,3)
  endif
  if(imt.eq.4) then ! rt source
     !tsgt(5) = tsgt(5) + uder(1,1)
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


  ttsgt(1) = tsgt(7)
  ttsgt(2) = tsgt(12)
  ttsgt(3) = tsgt(13)
  ttsgt(4) = tsgt(14)
  ttsgt(5) = tsgt(15)
  ttsgt(6) = tsgt(16)
  ttsgt(7) = tsgt(17)
  ttsgt(8) = tsgt(18)
  ttsgt(9) = tsgt(19)
  ttsgt(10)= tsgt(20)
  return
end subroutine udertotsgtSH

!

subroutine locallyCartesianDerivatives (u,udr,udt,udp,uder,r,theta)
  implicit none
  complex(kind(0d0)):: u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)
  real(kind(0d0)) :: r,theta
  real(kind(0d0)) :: thetasin,thetacot
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0


  if((theta.eq.0.d0).or.(theta.eq.pi)) then
     uder(1,1) = udr(1)
     uder(2,2) = u(1)/dcmplx(r)
     uder(3,3) = u(1)/dcmplx(r)     
     return
  endif


  thetasin = sin(theta)
  thetacot = cos(theta)/thetasin

  ! 1,2,3: r,theta,phi; , denotes the partial derivatives

  uder(1,1) = udr(1)
  uder(1,2) = (udt(1)-u(2))/dcmplx(r)
  uder(1,3) = (udp(1)/dcmplx(thetasin)-u(3))/dcmplx(r)
  uder(2,1) = udr(2)
  uder(2,2) = (udt(2)+u(1))/dcmplx(r)
  uder(2,3) = (udp(2)/dcmplx(thetasin)-u(3)*dcmplx(thetacot))/dcmplx(r)
  uder(3,1) = udr(3)
  uder(3,2) = udt(3)/dcmplx(r)
  uder(3,3) = (udp(3)/dcmplx(thetasin)+u(1)+u(2)*dcmplx(thetacot))/dcmplx(r)
  return
end subroutine locallyCartesianDerivatives
  

subroutine normalisetoKM(u,r)
  implicit none
  complex(kind(0d0)) :: u(1:3)
  integer :: i
  real(kind(0d0)) :: r
  do i = 1,3
     u(i) = u(i) / dcmplx(r)
  enddo
  return
end subroutine normalisetoKM



subroutine calgrid( nzone,vrmin,vrmax,vs,rmin,rmax, &
     & imax,lmin,tlen,vmin,gridpar,dzpar )
  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0 
  integer :: nzone,imax,lmin
  real(kind(0d0)) :: vrmin(*),vrmax(*),vs(4,*)
  real(kind(0d0)) :: rmin,rmax,tlen,vmin(*),gridpar(*),dzpar(*)
  integer :: izone,i,j
  real(kind(0d0)) :: coef1,coef2,v(4),vs1,vs2,rh,omega,amax,gtmp
  do izone=1,nzone
!     computing the S-velocity at each zone
     do i=1,4
        v(i) = vs(i,izone)
     enddo
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
     !     computing rh
     rh = vrmax(izone) - vrmin(izone)
     !      computing omega,amax
     omega = 2.d0 * pi * dble(imax) / tlen
     if ( vs1.ge.vs2 ) then
        vmin(izone) = vs2
     else
        vmin(izone) = vs1
     endif
     amax = vrmax(izone)
     gtmp = ( omega * omega ) / ( vmin(izone) * vmin(izone) ) &
          &        - ( (dble(lmin)+0.5d0) * (dble(lmin)+0.5d0) )&
          &        / ( amax * amax )
     if ( gtmp.gt.0.d0 ) then
        dzpar(izone)   = dsqrt( 1.d0/gtmp )
        gridpar(izone) = rh / dzpar(izone)
     else
        dzpar(izone)   = 0.d0
        gridpar(izone) = 0.d0
     endif
  enddo
  !     rearangement of gridpar
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

  !     re-rearangement of gridpar
  gtmp = 0.d0
  do izone=1,nzone
     gtmp = gtmp + gridpar(izone)
  enddo
  do izone=1,nzone
     gridpar(izone) = gridpar(izone) / gtmp
  enddo

  
end subroutine calgrid

!-----------------------------------------------------------------------------

subroutine calra(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,re )

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0 
  integer :: nlayer,inlayer
  integer :: nzone, nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax
  integer :: izone,itmp,i,ntmp
  real(kind(0d0)) :: rh,re

  inlayer = 0
  itmp = 1
  do izone=1,nzone
     rh = vrmax(izone) - vrmin(izone)

     if(dzpar(izone).eq.0.d0) then
        ntmp = 1
     else
        ntmp = int( sqrt(3.3d0 / re ) * rh / dzpar(izone) &
             &                    / 2.d0 / pi  / 7.d-1 + 1 )
     endif
     ! ntmp (see Geller & Takeuchi 1995 6.2)     
     nnl(izone) = ntmp
     if ( nnl(izone).lt.5 ) nnl(izone)=5
  enddo

  inlayer = 0
  do izone=1,nzone
     inlayer = inlayer + nnl(izone)
  enddo
  nlayer = inlayer 

end subroutine calra


!-----------------------------------------------------------------------------

subroutine calra2(nlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nnl,ra,re,nsta,rsta,rrsta, iista)

  implicit none
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0 
  integer :: nlayer,inlayer
  integer :: nzone,nnl(nzone)
  real(kind(0d0)) :: gridpar(*),dzpar(*),vrmin(*),vrmax(*),rmin,rmax,ra(nlayer+nzone+1)
  integer :: izone,itmp,i,nsta,ista, iista(1:3,1:nsta), ciista
  real(kind(0d0)) :: rsta(1:nsta),rrsta(1:3,1:nsta)
  real(kind(0d0)) :: rh,re,rs,ctmp
  !real(kind(0d0)) :: chikasa

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
                 rrsta(1,ista) = ra(itmp)
                 rrsta(2,ista) = ra(itmp+1)
                 rrsta(3,ista) = ra(itmp+2)
                     
                 iista(1,ista) = itmp
                 iista(2,ista) = itmp + 1
                 iista(3,ista) = itmp + 2
              else
                 rrsta(1,ista) = ra(itmp-1)
                 rrsta(2,ista) = ra(itmp)
                 rrsta(3,ista) = ra(itmp+1)
                 
                 iista(1,ista) = itmp - 1
                 iista(2,ista) = itmp
                 iista(3,ista) = itmp + 1
              endif

              !if((abs(rs-rsta(ista)).lt.ctmp).and.(abs(rs-rsta(ista)).ge.chikasa)) then
              !   ciista = ista
              !   ctmp = abs(rs-rsta(ista))
              !endif
           endif          
        enddo
        itmp = itmp + 1     
     enddo
  enddo

  if(ciista.eq.0) ciista = 1
  
end subroutine calra2


!-----------------------------------------------------------------------------


subroutine calsp (ndc, nlayer, isp, jsp)
  implicit none
  integer :: ndc,nlayer(*)
  integer :: isp(*),jsp(*)
  integer :: i

  ! computation of isp,jsp,ksp,lsp
  isp(1) = 1
  jsp(1) = 1
  do i=1,ndc
     isp(i+1) = isp(i) + nlayer(i)
     jsp(i+1) = jsp(i) + 4 * nlayer(i)
  enddo

  return
end subroutine calsp

!-----------------------------------------------------------------------------

subroutine calspo( ndc,rdc,nlayer,r0,rmin,rmax,ra, isp,spo,spn )
  ! computation de la source
  implicit none
  integer :: ndc,nlayer,isp(*),spn
  real(kind(0d0)) :: rdc(*),r0,rmin,rmax,ra(*),spo
  integer :: itmp
  

  ! checquer des parameters
  if ( (r0.lt.rmin).or.(r0.gt.rmax) ) &
       & pause 'The source location is improper.(calspo)'
  ! computing 'spo'
  if ( r0.eq.rmax ) then
     spo = dble(nlayer) - 0.01d0
     r0 = ra(nlayer) &
          & + (spo-dble(nlayer-1)) * ( ra(nlayer+1)-ra(nlayer) )   
  else
     do itmp = 2, ndc+nlayer+2        
        if ( r0.lt.ra(itmp) ) exit
     enddo
     spo = dble(itmp-2) + ( r0-ra(itmp-1) ) / ( ra(itmp)-ra(itmp-1) )
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
  do itmp = 1, ndc+1    
     spn = spn + 1
     if ( r0.le.rdc(itmp) ) exit
  enddo
  ! changing 'spo'
  spo = spo - dble( isp(spn) - 1 )
  


end subroutine calspo


!-----------------------------------------------------------------------------


subroutine calgra( isp,ra,r0,spn,spo,gra )
  
  integer :: isp(*),spn,itmp
  real(kind(0d0)) :: ra(*),r0,spo,gra(*)
  
  itmp = isp(spn) + dint( spo )
  gra(1) = ra(itmp)
  gra(2) = r0
  gra(3) = ra(itmp+1)
  
end subroutine calgra

!-----------------------------------------------------------------------------

subroutine calstg( nzone,rrho,vsv,vsh,nlayer,nnl,ra,rmax,vnp,vra,rho,ecL,ecN )
  implicit none
  integer:: nzone,nlayer,nnl(nzone),vnp
  real(kind(0d0)) :: rrho(4,nzone),vsv(4,nzone),vsh(4,nzone),ra(nlayer+nzone+1),rmax
  real(kind(0d0)) :: vra(nlayer+2*nzone+1),rho(nlayer+2*nzone+1)
  real(kind(0d0)) :: ecL(nlayer+2*nzone+1),ecN(nlayer+2*nzone+1)
  real(kind(0d0)) :: trho,tvsv,tvsh,coef
  integer :: izone,i,j,itmp,jtmp

  vra = 0
  rho = 0
  ecL = 0
  ecN = 0

  itmp = 0
  jtmp = 0
  do izone=1,nzone
     do i=1,nnl(izone)+1
        itmp = itmp + 1
        jtmp = jtmp + 1
        vra(itmp) = ra(jtmp)
        ! --- evaluating the density and elastic constants at this point
        trho = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        do j=1,4
           if ( j.eq.1 ) then
              coef = 1.d0
           else
              coef = coef * ( vra(itmp) / rmax )
           endif
           trho = trho + rrho(j,izone) * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
        enddo
        rho(itmp) = trho 
        ecL(itmp)  = rho(itmp) * tvsv * tvsv
        ecN(itmp)  = rho(itmp) * tvsh * tvsh

     enddo
     jtmp = jtmp - 1
  enddo
  vnp = itmp
end subroutine calstg


!-----------------------------------------------------------------------------

subroutine calgstg(nzone,nlayer,spn,rrho,vsv,vsh, ra,vra,rmax,rho,ecL,ecN,r0,mu0 )

  implicit none
  integer :: nlayer, nzone
  integer :: spn
  real(kind(0d0)) :: rrho(4,nzone),vsv(4,nzone),vsh(4,nzone)
  real(kind(0d0)) :: ra(3),rmax
  real(kind(0d0)) :: vra(3),rho(3)
  real(kind(0d0)) :: ecL(3),ecN(3)
  real(kind(0d0)) :: r0,mu0
  real(kind(0d0)) :: trho,tvsv,tvsh,coef
  integer :: i,j

  vra = 0.d0
  rho = 0.d0
  ecL = 0.d0
  ecN = 0.d0

  do i=1,3
     vra(i) = ra(i)
     trho = 0.d0
     tvsv = 0.d0
     tvsh = 0.d0
     do j=1,4
        if ( j.eq.1 ) then
           coef = 1.d0
        else
           coef = coef * ( vra(i) / rmax )
        endif
        trho = trho + rrho(j,spn) * coef
        tvsv  = tvsv  + vsv(j,spn)   * coef
        tvsh  = tvsh  + vsh(j,spn)   * coef
     enddo
     rho(i) = trho
     ecL(i)  = rho(i) * tvsv * tvsv
     ecN(i)  = rho(i) * tvsh * tvsh
    
  enddo

  mu0 = ecL(2)

  return

end subroutine calgstg



!-----------------------------------------------------------------------------

subroutine callsuf(omega,nzone,vrmax,vsv,lsuf)

  implicit none
  integer :: nzone, lsuf
  real(kind(0d0)) :: omega, vrmax(*), vsv(4,*)
  real(kind(0d0)) :: tvs, coef
  integer :: i
  
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
end subroutine callsuf


!-----------------------------------------------------------------------------


subroutine calcoef( nzone,omega,q,coef )
  implicit none
  real(kind(0d0)), parameter ::  pi = 3.1415926535897932d0   
  integer :: izone,nzone
  real(kind(0d0)) :: omega,q(*)
  complex(kind(0d0)) :: coef(*)
  real(kind(0d0)) :: aa,bb
  

  

  do izone=1,nzone
     if(q(izone).le.0.d0) then
        coef(izone) = dcmplx(1.d0)
     else
        if ( omega.eq.0.d0 ) then
           aa = 1.d0
        else
           aa = 1.d0 + dlog( omega / ( 2.d0 * pi ) ) / ( pi * q(izone) )
        endif
        bb = 1.d0 / ( 2.d0 * Q(izone) )
        coef(izone) = dcmplx( aa, bb ) * dcmplx( aa, bb )
  
     endif
  enddo
     


end subroutine calcoef


!-----------------------------------------------------------------------------

subroutine calamp(g,l,lsuf,maxamp,ismall,ratl)
  
  implicit none
  integer :: l,lsuf,ismall
  real(kind(0d0)) :: maxamp,ratl
  complex(kind(0d0)) :: g
  real(kind(0d0)) :: amp,ampratio    
  ampratio = 0.d0
  amp = abs(g)
  if( amp.gt.maxamp ) maxamp = amp
  if ( (amp.ne.0.d0).and.(maxamp.ne.0.d0) ) then
     ampratio = amp / maxamp
  endif
  if( (ampratio.lt.ratl).and.(l.gt.lsuf) ) then
     ismall = ismall + 1
  else
     ismall = 0
  endif
  
end subroutine calamp


!-----------------------------------------------------------------------------

subroutine calu(c0,lsq,bvec,u)
  implicit none
  real(kind(0d0)) :: lsq
  complex(kind(0d0)) :: c0, bvec(3), u(1:3)
  
  u(1) = dcmplx( 0.d0 )
  u(2) = u(2) + c0 * bvec(2) / dcmplx(lsq)
  u(3) = u(3) + c0 * bvec(3) / dcmplx(lsq)

end subroutine calu


!-----------------------------------------------------------------------------

subroutine interpolate( ncomp,nderiv,rsta,rrsta,g,u )

  implicit none
  integer :: ncomp,nderiv
  real(kind(0d0)) :: rsta,rrsta(3)
  complex(kind(0d0)) :: g(3*ncomp),u(ncomp)    
  real(kind(0d0)):: dh(3)     
  integer :: ip(3),ier,i,itmp,icomp
  complex(kind(0d0)) :: a(3,3),b(3),wk(3)
  real(kind(0d0)) :: eps 
  data eps / -1.d0 /

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

  return
end subroutine interpolate


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




subroutine fillinpb( nderiv,b )

  integer :: nderiv
  complex(kind(0d0)) :: b(3)
     
  if( (nderiv.ne.0).and.(nderiv.ne.1).and.(nderiv.ne.2) ) &
       &     pause 'invalid argument (fillinpb)'
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



!------------------------------------------------------------------------------


subroutine calcutd(nzone,nnl,tmpr,rat,nn,ra,kc)

  implicit none
  integer :: nzone,nn,spn,kc,nnl(*)
  complex(kind(0d0)) :: tmpr(*)
  real(kind(0d0)) :: rat,ra(*)
  integer  :: nc
  real(kind(0d0)):: cU(nn),rc
  real(kind(0d0)) :: maxamp,amp(nn)
  integer :: iz,jz,jj,i,ml(nzone),tzone
  
  do jj=1,nn
     cU(jj) = 0.d0
  enddo
  
  iz = 2
  jz = 1
  do jj=1,nn
     cU(jj) = tmpr(jj)
  enddo

  maxamp = -1.d0
  do i=1,nn
     amp(i) = cU(i)
     if(maxamp.lt.amp(i)) maxamp = amp(i)
  enddo
  maxamp = maxamp * rat ! threshold value
  if(maxamp.eq.0.d0) then
     kc = 1
     return
  endif
  
  do i=1,nn
     if(amp(i).gt.maxamp) then
        nc = i
        exit
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
  kc = nc
  

  return
end subroutine calcutd

