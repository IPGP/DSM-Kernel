subroutine pinputKernel(kernelinf)

  use parameters

  implicit none
  character(120) :: tmpfile,kernelinf
  character(120) :: dummy
  integer :: iitmp,maxlmax

  tmpfile='tmpworkingfile_for_SynViewer'

  open(unit=2, file=kernelinf,status='unknown')

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
  close(1)

  call pinputDSM(DSMconfFile,PoutputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch,SGTinfo)
  call readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
  tmpfile = 'tmpworkingfile_for_psvmodel'
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
  character(120), parameter :: tmpfile='tmpworkingfile_for_SGTcalcul'
  character(120) :: dummy,outputDir,psvmodel,modelname,DSMconfFile,SGTinfo
  real(kind(0d0)) :: tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta
  real(kind(0d0)) :: thetamin,thetamax,thetadelta
  integer :: imin,imax,rsgtswitch,tsgtswitch,synnswitch

  
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
  read(1,*) rsgtswitch,tsgtswitch,synnswitch
  close(1)

end subroutine pinputDSM

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
