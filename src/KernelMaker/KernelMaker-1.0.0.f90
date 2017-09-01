program KernelMaker


  !
  !   
  !  KernelMaker
  !                                      by FUJI Nobuaki 2012
  !
  !       many subroutines from Li Zhao's code for normal mode
  !   
  !     0.3.1    video mode             Nobuaki Fuji 2014.3.
  !
  !     1.0.0    for new file format (compatible with SGTpsv/sh >1.0.0)
  !              


  use parameters
  use tmpSGTs
  use angles
  use kernels
  use rotate
  
  implicit none
  include 'mpif.h'
  include '../../etc/config_calcul/constants.h'

  ! file name

  character(200) :: infofile, gridfile, kerfile, kertotalfile
  character(120) :: tmpchar
  integer :: i,icheck,jt,it,k,j
  character(40) :: datex,timex
  character(120) :: list
  real(kind(0d0)) :: tmparray(1:9) ! for model parameters

  real(kind(0d0)) :: thetasgcs,phisgcs
  real(kind(0d0)) :: thetargcs,phirgcs
  real(kind(0d0)) :: distan,azim,bazim

  integer, external :: getpid  

  ! ignoring scheme
  integer :: iitype,itype
  integer(2), allocatable :: iflagForRapidity(:),iflagForRapidityOld(:), iflagForRapidityNext(:)
  integer(2), allocatable :: iflagForRapidity0(:),modiflag(:)
  integer :: searchiteration,iiphi
  
  ! filter
  integer :: ift
  
  ! the output
  real(kind(0e0)), allocatable :: totalker(:,:,:,:,:)
  integer :: kc,idum
  real(kind(0e0)) :: fdum


  ! number of components ! attention! these definitions are also necessary in calculateKenel !!
  integer, parameter :: num_tsgtPSV = 20
  integer, parameter :: num_rsgtPSV = 10
  integer, parameter :: num_synnPSV = 10
  integer, parameter :: num_tsgtSH = 10
  integer, parameter :: num_rsgtSH = 5
  integer, parameter :: num_synnSH = 5
  integer, parameter :: num_h3 = 6
  integer, parameter :: num_h4 = 6


  ! for MPI
  integer :: nproc, my_rank, ierr
  !    
  !-----------------------------------------------------------------------

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
  
  if(my_rank.eq.0) then
     ! Inputting parameters
     call pinputKernel
     ! DSM constants calculation
     omegai = - dlog(omegai) / tlen
     if(ifastFFT.ne.1) then
        fmin = imin
        fmax = imax
     endif     
  endif


  call MPI_BCAST(SGTinfo,    120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(parentDir,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(eventName,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stationName,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(phase,      120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(compo,      120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(paramWRT,   120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(synnfile,   120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(INFO_TSGT,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(INFO_RSGT,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rsampletxt, 120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(modelcard,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Poutputdir, 120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(psvmodel,   120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(modelname,  120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(DSMconfFile,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(evla,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(evlo,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(evdepth,   1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stla,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(stlo,      1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mt,        6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(twin,      4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ipdistance,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(c_red_reci,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ibwfilt,      1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(npButterworth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ifastFFT,     1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fmin,         1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fmax,         1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fclp,1+nfilter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fchp,1+nfilter,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(freqid,120*(1+nfilter),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dph,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ph1,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dth,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(thw,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(start,     1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(end,       1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iWindowStart, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iWindowEnd,   1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(samplingHz,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(calculrapide,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nntype,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iPSVSH,       1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  mtype=0
  if(compo.eq.'Z') then
     mtype=mtype+10
  elseif(compo.eq.'R') then
     mtype=mtype+20
  elseif(compo.eq.'T') then
     mtype=mtype+30
  else
     compo='X'
     mtype=99
  endif
  if(my_rank.ne.0) then
     allocate(minici(1:nntype,0:nfilter))
     allocate(idecidetype(1:nntype))
  endif
  call MPI_BCAST(idecidetype,nntype,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  ! minici should be for each node
  
  ! exporting DSM parameters
  call MPI_BCAST(re,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ratc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ratl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omegai,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rmin_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rmax_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rdelta_,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(my_rank.eq.0) then
     r_n =  int((rmax_-rmin_)/rdelta_)+1
  endif
  call MPI_BCAST(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(r_(1:r_n))
    if(my_rank.eq.0) then
     do i = 1, r_n
        r_(i) = rmin_ + dble(i-1)*rdelta_
     enddo
  endif
  call MPI_BCAST(r_,r_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  

  if(my_rank.eq.0) then
     if((rmin.eq.0.d0).and.(rmax.eq.0.d0).and.(rdelta.eq.0.d0)) then
        rmin=rmin_
        rmax=rmax_
        rdelta=rdelta_
     endif
  endif
  call MPI_BCAST(rmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rdelta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(my_rank.eq.0) then
     nr =  int((rmax-rmin)/rdelta)+1
  endif
  call MPI_BCAST(nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(r(1:nr))
    if(my_rank.eq.0) then
     do i = 1, nr
        r(i) = rmin + dble(i-1)*rdelta
     enddo
  endif
  call MPI_BCAST(r,nr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
  ! check r(:) and r_(:) because we don't interpolate for depth 
  
  do ir=1,nr
     icheck=0
     do i=1,r_n
        if(r_(i).eq.r(ir)) then
           icheck=1
        endif
     enddo
     if(icheck.eq.0) then
        print *, r(ir), "is not in the catalogue, sorry"
        stop
     endif
  enddo
  
  
     

  call MPI_BCAST(r0min,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(r0delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(my_rank.eq.0) then
     r0_n =  int((r0max-r0min)/r0delta)+1
  endif
  call MPI_BCAST(r0_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  
  allocate(r0D(1:r0_n))
  if(my_rank.eq.0) then
     do i = 1, r0_n
        r0D(i) = r0min + dble(i-1)*r0delta
     enddo
  endif
  call MPI_BCAST(r0D,r0_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)   
  ir0 = r0_n
  
  call MPI_BCAST(thetamin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(thetamax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(thetadelta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  if(my_rank.eq.0) then
     theta_n = int((thetamax-thetamin)/thetadelta)+1
  endif
  call MPI_BCAST(theta_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(thetaD(1:theta_n))
  if(my_rank.eq.0) then   
     do i = 1,theta_n
        thetaD(i) = (thetamin + dble(i-1)*thetadelta)
     enddo
  endif  
  call MPI_BCAST(thetaD,theta_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  !-----------------------------------------------------------------------
  ! model parameters
  if(my_rank.eq.0) then
     write(psvmodel,"(Z4)") getpid()
     psvmodel = 'tmpworkingfile_for_psvmodel'//psvmodel
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     close(20)
  endif
  call MPI_BCAST(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
  if(my_rank.eq.0) then
     write(psvmodel,"(Z4)") getpid()
     psvmodel = 'tmpworkingfile_for_psvmodel'//psvmodel  
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     do i = 1, nzone
        read (20, *) vrminD(i), vrmaxD(i), rrhoD(1,i), rrhoD(2,i), rrhoD(3,i), rrhoD(4,i), vpvD(1,i), vpvD(2,i), vpvD(3,i), vpvD(4,i), vphD(1,i), vphD(2,i), vphD(3,i), vphD(4,i), vsvD(1,i), vsvD(2,i), vsvD(3,i), vsvD(4,i), vshD(1,i), vshD(2,i), vshD(3,i), vshD(4,i), etaD(1,i), etaD(2,i), etaD(3,i), etaD(4,i), qmuD(i), qkappaD(i)
     enddo
     close(20,status='delete')
  endif
  call MPI_BCAST(vrminD,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vrmaxD,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rrhoD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vpvD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vphD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vsvD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vshD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(etaD,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qmuD,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qkappaD,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  

  rEarth=vrmaxD(nzone) 

  !-----------------------------------------------------------------------

  ! Record the date and time at the beginning of the job
  if(my_rank.eq.0) then
     list = trim(parentDir)//"/log/calLog"//"."// &
          trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"."//trim(paramWRT)//".log"
     open(1,file =list, status = 'unknown', form = 'formatted')
     call date_and_time(datex,timex)
     write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
          '    Starting date and time:                     ', &
          datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
          timex(1:2),':',timex(3:4),':',timex(5:8)   
     close (1)
  endif
  
  !-----------------------------------------------------------------------   
  ! lsmoothfinder for FFT
  np0=fmax
  call lsmoothfinder(tlen,np0,samplingHz,lsmooth)
  i=1
  do while (i<lsmooth)
     i = i*2
  enddo
  lsmooth = i
  i = 0
  np1 = 1
  do while (np1<np0)
     np1 = np1*2
  enddo  
  np1 = np1*lsmooth    
  
  ! redefinition of samplingHz
  samplingHz = dble(2*np1)/tlen
  dtn = 1.d0/samplingHz
  iWindowStart = int(start*samplingHz)
  iWindowEnd   = int(end*samplingHz)

  

  ! allocate SGTs, synthetics in frequency
  allocate(omega(fmin:fmax))
  do i = fmin, fmax
     omega(i) = 2.d0*pi*dble(i)/tlen
  enddo

  allocate(tmph01(fmin:fmax))
  allocate(tmph02(fmin:fmax))
  allocate(tmph03(fmin:fmax))
  allocate(tmph04(fmin:fmax))
  allocate(tmph05(fmin:fmax))
  allocate(tmph06(fmin:fmax))
  allocate(tmph07(fmin:fmax))
  allocate(tmph08(fmin:fmax))
  allocate(tmph09(fmin:fmax))
  allocate(tmph10(fmin:fmax))

  allocate(tmph11(fmin:fmax))
  allocate(tmph12(fmin:fmax))
  allocate(tmph13(fmin:fmax))
  allocate(tmph14(fmin:fmax))
  allocate(tmph15(fmin:fmax))
  allocate(tmph16(fmin:fmax))
  allocate(tmph17(fmin:fmax))
  allocate(tmph18(fmin:fmax))
  allocate(tmph19(fmin:fmax))
  allocate(tmph20(fmin:fmax))

  allocate(tmph21(fmin:fmax))
  allocate(tmph22(fmin:fmax))
  allocate(tmph23(fmin:fmax))
  allocate(tmph24(fmin:fmax))
  allocate(tmph25(fmin:fmax))
  allocate(tmph26(fmin:fmax))


  ! allocate vectors in time domain
  allocate(t(iWindowStart:iWindowEnd))
  allocate(u(iWindowStart:iWindowEnd))
  allocate(u0(0:nfilter,iWindowStart:iWindowEnd))
  allocate(v(iWindowStart:iWindowEnd))
  allocate(v0(0:nfilter,iWindowStart:iWindowEnd))
  allocate(hu(iWindowStart:iWindowEnd))
  allocate(hu0(0:nfilter,iWindowStart:iWindowEnd))
  allocate(fwin(0:nfilter,iWindowStart:iWindowEnd))  
  allocate(nt1(0:nfilter))
  allocate(nt2(0:nfilter))
  allocate(denomv(0:nfilter))
  allocate(denomu(0:nfilter))
  
  do i = iWindowStart,iWindowEnd
     t(i) = dble(i)*dtn
  enddo


  ! extract the reference model parameters on the radial grid points
  allocate(rhom(1:nr))
  allocate(vpm(1:nr))
  allocate(vsm(1:nr))
  allocate(qmm(1:nr))
  allocate(qkp(1:nr))
  do ir = 1,nr
     call calstg_for_card(r(ir),nzone,vrminD,vrmaxD,rrhoD,vpvD,vphD,vsvD,vshD,etaD,qmuD,qkappaD,tmparray)
     rhom(ir) = tmparray(2)
     vpm(ir)  = 5.d-1*tmparray(3)+5.d-1*tmparray(7)
     vsm(ir)  = 5.d-1*tmparray(4)+5.d-1*tmparray(8)
     qmm(ir)  = tmparray(6)
     qkp(ir)  = tmparray(5)
     if(qmm(ir).le.0.d0) qmm(ir)  = 1.d5
     if(qkp(ir).le.0.d0) qkp(ir)  = 1.d5
  enddo
  
  
  ! source check
  
  icheck=0
  do ir0 = 1,r0_n
     if(rEarth-evdepth.eq.r0D(ir0)) then
        icheck =1
     endif
  enddo
  if(icheck.eq.0) then
     print *, "depth",evdepth, "is not in the catalogue",Poutputdir,",sorry"
     stop
  endif
  
  sdep=evdepth
  rs=rEarth-evdepth
  rr=rEarth-rdep
  
  ! Convert the event latitude (0 to +/- pi/2) to co-latitude (0 to pi)
  !  and psuedo-longitude (0 to +/- pi) to longitude (0 to 2*pi).
  
  slat = evla
  slon = evlo
  rlat  = stla
  rlon  = stlo
  
  !thetasgcs0=(90.d0-slat)*pi/180.d0
  !phisgcs0=slon*pi/180.d0
  !if(slon.lt.0.d0) phisgcs0=(360.d0+slon)*pi/180.d0
  
  thetargcs=(90.d0-rlat)*pi/180.d0
  phirgcs=rlon*pi/180.d0
  if(rlon.lt.0.d0) phirgcs=(360.d0+rlon)*pi/180.d0
  
  thetasgcs=(90.d0-slat)*pi/180.d0
  phisgcs=slon*pi/180.d0
  if(slon.lt.0.d0) phisgcs=(360.d0+slon)*pi/180.d0


  slat=90.d0-thetasgcs*180.d0/pi
  slon=phisgcs*180.d0/pi

  ! Calculate the station epicentral distance and azimuth
  call azimth(0,slat,slon,rlat,rlon,distan,azim,bazim)
  call rotmat(thetasgcs,phisgcs,thetargcs,phirgcs,azim)

  
  if(my_rank.eq.0) then
     list = trim(parentDir)//"/log/calLog"//"."// &
          trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"."//trim(paramWRT)//".log"
     open(1,file =list, status = 'old', access='append',form = 'formatted')
     
     write(1,'(a,4(1x,f8.4))') '    Source and distance: ', &
     	  slat,slon,sdep,distan
     write(1,'(a,4(1x,f8.4))') '    Receiver location:            ', &
     	  rlat,rlon
     write(1,'(a,2(1x,f8.4))') '    azimuth and back-azimuth:     ', &
     	  azim,bazim
     close(1)
  endif
  
  ! Of course MinLengthFor0 can be chosen but nMinLengthFor0=2 works well so I fixed this value
  MinLengthFor0 = dph*2.d0
  nMinLengthFor0=int(MinLengthFor0/dph)
  ! see the explication above

  ! Creating horizontal grid 
  ! and converting horizontal grid from path-specific system to geographic system
  call hgridsimple(distan,ph1,dph,thw,dth)
  allocate(iflagForRapidity(1:nphi))
  allocate(iflagForRapidityOld(1:nphi))
  allocate(iflagForRapidityNext(1:nphi))
  allocate(iflagForRapidity0(1:nphi))
  allocate(modiflag(1:nphi-nMinLengthFor0+1))
  call convertPath2Geo
  ! Calculate the cosine and sine functions in SGT expressions for scattering points (ZC10a)
  allocate(crq(0:nphi,0:ntheta),crq2(0:nphi,0:ntheta))
  allocate(srq(0:nphi,0:ntheta),srq2(0:nphi,0:ntheta))
  allocate(csq(0:nphi,0:ntheta),csq2(0:nphi,0:ntheta))
  allocate(ssq(0:nphi,0:ntheta),ssq2(0:nphi,0:ntheta))
  allocate(cqs(0:nphi,0:ntheta),cqs2(0:nphi,0:ntheta))
  allocate(sqs(0:nphi,0:ntheta),sqs2(0:nphi,0:ntheta))
  allocate(deltar(nphi,ntheta),deltas(nphi,ntheta))
  call calculateSineCosine(azim)


  ! allocation for catalogues
  allocate(tsgtomega(1:num_tsgtPSV,fmin:fmax,1:theta_n))
  allocate(rsgtomega(1:num_rsgtPSV,fmin:fmax,1:theta_n))
  allocate(synnomega(1:num_synnPSV,fmin:fmax,1:theta_n))
  
  ! allocation for sgtF (interpolated functions)
  allocate(tsgtF(1:num_tsgtPSV,fmin:fmax))
  allocate(rsgtF(1:num_rsgtPSV,fmin:fmax))
  allocate(synnF(1:num_synnPSV,fmin:fmax))
  allocate(h3(1:num_h3,fmin:fmax))
  allocate(h4(1:num_h4,fmin:fmax))
  allocate(u_freq(fmin:fmax))
  
  ! Now computing reference synthetic (in my_rank = 0)

  synnomega=cmplx(0.d0)
  tsgtomega=cmplx(0.d0)
  rsgtomega=cmplx(0.d0)
  h3=cmplx(0.d0)
  h4=cmplx(0.d0)
  u=0.d0
  u0=0.d0
  if(my_rank.eq.0) then
     
     if(iPSVSH.ne.1) call rdsgtomega(rs,0.d0,num_synnPSV,num_synnPSV,2)
     if(iPSVSH.ne.2) call rdsgtomega(rs,0.d0,num_synnSH,num_synnPSV,1)
     call clsgt(distan,num_synnPSV,synnF(1:num_synnPSV,fmin:fmax),synnomega(1:num_synnPSV,fmin:fmax,1:theta_n))
     
     call synn2h3freq(0,0)

     u_freq(fmin:fmax)=h3(1,fmin:fmax)*cmplx(mt(1))+ h3(2,fmin:fmax)*cmplx(mt(2)) &
          +h3(3,fmin:fmax)*cmplx(mt(3))+ 2.d0*(h3(4,fmin:fmax)*cmplx(mt(4)) &
          +h3(5,fmin:fmax)*cmplx(mt(5))+ h3(6,fmin:fmax)*cmplx(mt(6)))
     call vectorFFT_double(fmin,fmax,np1,u_freq(fmin:fmax),u(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)
   
     !  Numerically differentiate displacement to obtain velocity response.

     v=0.d0
     do it=iWindowStart+1,iWindowEnd
        v(it)=(u(it)-u(it-1))/dtn
     enddo
     
     ! Calculate the Hilbert transform of displacement

     hu=0.d0
     call hilbert(iWindowEnd-iWindowStart,dtn,t,u(iWindowStart:iWindowEnd),hu(iWindowStart:iWindowEnd))     

     if(ibwfilt) then
        do ift = 0,nfilter
           call bwfilt(u(iWindowStart:iWindowEnd),u0(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
           call bwfilt(v(iWindowStart:iWindowEnd),v0(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
           call bwfilt(hu(iWindowStart:iWindowEnd),hu0(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
        enddo
     endif

     do ift = 0,nfilter
        open(1,file=trim(synnfile)//"."//trim(freqid(ift)),status='unknown',form='formatted')
        do jt=iWindowStart,iWindowEnd
           write(1,*) t(jt), u(jt), u0(ift,jt)
        enddo
        close(1)
     enddo
     
  endif
  
  if(trim(paramWRT).eq.'test') then
     if(my_rank.eq.0) then
        do jt=iWindowStart,iWindowEnd
           write(13,*) t(jt), u(jt), u0(0,jt)
        enddo
     endif
     call MPI_FINALIZE(ierr)
     stop
  endif
  
  
  call MPI_BCAST(u0,(iWindowEnd-iWindowStart+1)*(nfilter+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(u,(iWindowEnd-iWindowStart+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(v0,(iWindowEnd-iWindowStart+1)*(nfilter+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(v,(iWindowEnd-iWindowStart+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  call MPI_BCAST(hu0,(iWindowEnd-iWindowStart+1)*(nfilter+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(hu,(iWindowEnd-iWindowStart+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
  call fwinDeterminator

  ! extracting "the phase"

  u0(0:nfilter,iWindowStart:iWindowEnd)=u0(0:nfilter,iWindowStart:iWindowEnd)*fwin(0:nfilter,iWindowStart:iWindowEnd)
  v0(0:nfilter,iWindowStart:iWindowEnd)=v0(0:nfilter,iWindowStart:iWindowEnd)*fwin(0:nfilter,iWindowStart:iWindowEnd)
  hu0(0:nfilter,iWindowStart:iWindowEnd)=hu0(0:nfilter,iWindowStart:iWindowEnd)*fwin(0:nfilter,iWindowStart:iWindowEnd)
 
  ! Calculate the denominator in the kernel expressions
  
  call coeffCalculator
  ! for radial and azimuthal anisotropy
  if((sym.ge.0.d0).and.(sym.le.360.d0)) then
     csym=dcos(sym*pi/180.d0)
     ssym=dsin(sym*pi/180.d0)
  endif
 

  if(my_rank.eq.0) then   
     list = trim(parentDir)//"/log/calLog"//"."// &
           trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"."//trim(paramWRT)//".log"         
     open(1,file =list, status = 'old',access='append', form = 'formatted')
     call date_and_time(datex,timex)
     write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
          '    kernel calculation started:                     ', &
          datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
          timex(1:2),':',timex(3:4),':',timex(5:8)   
     write(1,*) 'start kernel calculations Mtype= ',mtype
     close (1)
     
  endif

  ! for video mode, number of snapshots will be decided here

  number_of_snapshots = 1
  if((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
     jtstep_timeincrementV=1
     if(timeincrementV.ne.0.d0) jtstep_timeincrementV=int(timeincrementV*samplingHz) 
     if(timeincrementV.eq.0.d0) jtstep_timeincrementV=1
     number_of_snapshots = (iWindowEnd-iWindowStart+1)/jtstep_timeincrementV
  endif
  
  ! Now loop over all grid ponts to compute the kernels
  ! for the parallelisation, I devide nr into nproc

  allocate(tmpker(0:nktype,0:nfilter))

  if((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
     allocate(tmpvideoker(0:nkvtype,0:nfilter,1:number_of_snapshots))
     allocate(videoker(1:nphi,0:nkvtype,0:nfilter,1:number_of_snapshots)) ! be careful of the difference between ker and videoker (along theta)
  endif
     
  if(trim(paramWRT).eq.'vTSGT') then
     allocate(tmpvideoker(1:num_h4,0:nfilter,1:number_of_snapshots))
     allocate(videoker(1:nphi,1:num_h4,0:nfilter,1:number_of_snapshots))
  endif

  if(trim(paramWRT).eq.'vRSGT') then
     allocate(tmpvideoker(1:num_h3,0:nfilter,1:number_of_snapshots))
     allocate(videoker(1:nphi,1:num_h3,0:nfilter,1:number_of_snapshots))
  endif
  
  allocate(ker(1:nphi,1:ntheta,0:nktype,0:nfilter))

  tmpker = 0.d0
  ker = 0.e0
  allocate(du(iWindowStart:iWindowEnd))
  allocate(duf(0:nfilter,iWindowStart:iWindowEnd)) 
  allocate(duq(iWindowStart:iWindowEnd))
  allocate(duqf(0:nfilter,iWindowStart:iWindowEnd))
  ntot=nphi*ntheta*nr
  k=0  



  do ir=1,nr
     ! if-line for parallelisation
     if((ir.ne.0).and.((mod(nr-my_rank-ir,2*nproc).eq.0).or.(mod(nr+my_rank+1-ir,2*nproc).eq.0))) then
        rx=r(ir)
        
        synnomega=cmplx(0.d0)
        tsgtomega=cmplx(0.d0)
        rsgtomega=cmplx(0.d0)
        h3=cmplx(0.d0)
        h4=cmplx(0.d0)
        
        ! SSGT reading
        if(iPSVSH.ne.1) call rdsgtomega(rx,0.d0,num_rsgtPSV,num_rsgtPSV,20)
        if((iPSVSH.ne.2).and.(vsm(ir).ne.0.d0)) call rdsgtomega(rx,0.d0,num_rsgtSH,num_rsgtPSV,10)
        ! TSGT reading
        if(iPSVSH.ne.1) call rdsgtomega(rs,rx,num_tsgtPSV,num_tsgtPSV,200)
        if((iPSVSH.ne.2).and.(vsm(ir).ne.0.d0)) call rdsgtomega(rs,rx,num_tsgtSH,num_tsgtPSV,100)
       
        

 
        


        ! Video mode will calculate for every point

        if((trim(paramWRT).eq.'allV').or.(trim(paramWRT).eq.'alphaV').or. &
           (trim(paramWRT).eq.'betaV')) then
           calculrapide=0.d0
        endif 

        if((calculrapide.eq.0.d0).and.((trim(paramWRT).ne.'betaV').and.(trim(paramWRT).ne.'allV').and.(trim(paramWRT).ne.'alphaV'))) then ! calculate for every point
           do ith = 1,ntheta
              do ip=1,nphi
                 k=k+1
                 xlat=90.d0-theta(ip,ith)
                 xlon=phi(ip,ith)
                 call calculateKernel
                 ker(ip,ith,:,:)=tmpker(:,:)
                 
              enddo
              
           enddo
        elseif((trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV').or.(trim(paramWRT).eq.'alphaV')) then
           
            do ith = 1,ntheta
               videoker=0.e0
              do ip=1,nphi
                 k=k+1
                 xlat=90.d0-theta(ip,ith)
                 xlon=phi(ip,ith)
                 tmpvideoker=0.e0
                 call calculateKernel
                 
                 videoker(ip,:,:,:)=tmpvideoker(:,:,:)
              enddo
            
              !write(tmpchar,'(I7,".",I7)') int(rx*1.d3), int(xlat*1.d3)
              write(tmpchar,'(I7,".",I7)') ir,ith
              do j=1,15
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                   //"."//"timemarching"
              open(1,file=kerfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*nphi*(1+nkvtype)*(1+nfilter)*number_of_snapshots)
              write(1,rec=1) videoker(1:nphi,0:nkvtype,0:nfilter,1:number_of_snapshots)
              close(1)            
           enddo


        elseif(trim(paramWRT).eq.'vRSGT') then


           do ith=1,ntheta
              videoker=0.e0
              do ip=1,nphi
                 k=k+1
                 xlat=90.d0-theta(ip,ith)
                 xlon=phi(ip,ith)
                 tmpvideoker=0.e0
                 call calculateRSGT
                 videoker(ip,:,:,:)=tmpvideoker(:,:,:)
              enddo


              write(tmpchar,'(I7,".",I7)') ir,ith
              do j=1,15
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                   //"."//"timemarching"
              open(1,file=kerfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*nphi*(num_h3)*(1+nfilter)*number_of_snapshots)
              write(1,rec=1) videoker(1:nphi,1:num_h3,0:nfilter,1:number_of_snapshots)
              close(1)           

           enddo



        elseif(trim(paramWRT).eq.'vTSGT') then

           
           
           do ith=1,ntheta
              videoker=0.e0
              do ip=1,nphi
                 k=k+1
                 xlat=90.d0-theta(ip,ith)
                 xlon=phi(ip,ith)
                 tmpvideoker=0.e0
                 call calculateTSGT
                 videoker(ip,:,:,:)=tmpvideoker(:,:,:)
              enddo


              write(tmpchar,'(I7,".",I7)') ir,ith
              do j=1,15
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                   //"."//"timemarching"
              open(1,file=kerfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*nphi*(num_h4)*(1+nfilter)*number_of_snapshots)
              write(1,rec=1) videoker(1:nphi,1:num_h4,0:nfilter,1:number_of_snapshots)
              close(1)           

           enddo


        else ! we ignore the value smaller than calculrapide*max

           ! First obtain the values on the great circle
           ith=ntheta/2+1
           do ip=1,nphi
              k=k+1
              xlat=90.d0-theta(ip,ith)
              xlon=phi(ip,ith)
              call calculateKernel
              ker(ip,ith,:,:)=tmpker(:,:)
           enddo
           
           ith=ntheta/2+1 ! for making criteria for grid searching
           
           minici = 0.d0
           
           ! step 1
           do ift=0,nfilter
              do iitype = 1,nntype
                 itype = idecidetype(iitype)
                 minici(iitype,ift)=calculrapide*maxval(abs(ker(1:nphi,ith,itype,ift)))
              enddo
           enddo
           
           
           ! step 2 and 3
           
           ith=ntheta/2+1 
           iflagForRapidity(1:nphi) = 1
           do ift=0,nfilter 
              do ip = 1,nphi
                 do iitype = 1,nntype
                    itype = idecidetype(iitype)
                    if(abs(ker(ip,ith,itype,ift)).gt.minici(iitype,ift)) then
                       iflagForRapidity(ip) = iflagForRapidity(ip)*0
                    endif
                 enddo
              enddo
           enddo
           
           
           
           ! step 4
           
           modiflag(1:nphi-nMinLengthFor0+1)=1
           do ip=1,nphi-nMinLengthFor0+1
              do iiphi=ip,ip+nMinLengthFor0-1
                 modiflag(ip) = modiflag(ip)*iflagForRapidity(iiphi)
              enddo
              ! iflagForRapidity(ip:ip+nMinLengthFor0-1)=tmpiflagForRapidity
           enddo
           
           iflagForRapidity(1:nphi)=1
           do ip=1,nphi-nMinLengthFor0+1
              if(modiflag(ip).eq.0) then
                 iflagForRapidity(ip:ip+nMinlengthFor0-1)=modiflag(ip)    
              endif
           enddo
           
           
           iflagForRapidity0(1:nphi)=iflagForRapidity(1:nphi)
           
           ! step 5 for +side 
           do ith=ntheta/2+2,ntheta
              iflagForRapidityOld(1:nphi) = iflagForRapidity(1:nphi)
              do searchiteration = 1, int(nphi/nMinLengthFor0)+1
                 
                 if(minval(iflagForRapidity(1:nphi)).eq.1) exit
                 !print *,    searchiteration, int(nphi/nMinlengthFor0)+1
                 do ip = 1,nphi
                     if(iflagForRapidity(ip)) then
                       !ker(ip,ith,2,0)=1.d0       
                       !cycle ! don't calculate if 1
                    else
                       k=k+1
                       xlat=90.d0-theta(ip,ith)
                       xlon=phi(ip,ith)
                       call calculateKernel
                       ker(ip,ith,:,:)=tmpker(:,:)
                    endif
                 enddo
                 
                 
                 
                 ! step 2 and 3 in step5
                 iflagForRapidity(1:nphi) = 1
                 do ift=0,nfilter 
                    do ip = 1,nphi
                       do iitype = 1,nntype
                          itype = idecidetype(iitype)
                          if(abs(ker(ip,ith,itype,ift)).gt.minici(iitype,ift)) then
                             iflagForRapidity(ip) = iflagForRapidity(ip)*0
                          endif
                       enddo
                    enddo
                 enddo
                 
                 ! step 4 in step5
                 
                 
                 
                 modiflag(1:nphi-nMinLengthFor0+1)=1
                 do ip=1,nphi-nMinLengthFor0+1
                    do iiphi=ip,ip+nMinLengthFor0-1
                       modiflag(ip) = modiflag(ip)*iflagForRapidity(iiphi)
                    enddo
                    ! iflagForRapidity(ip:ip+nMinLengthFor0-1)=tmpiflagForRapidity
                 enddo
                 
                 
                 
                 iflagForRapidity(1:nphi)=1
                 do ip=1,nphi-nMinLengthFor0+1
                    if(modiflag(ip).eq.0) then
                       iflagForRapidity(ip:ip+nMinlengthFor0-1)=modiflag(ip)    
                    endif
                 enddo
                 
                 
                 
                 
                 iflagForRapidityNext(1:nphi) = iflagForRapidity(1:nphi)                   
                 
                 do ip=1,nphi
                    if((iflagForRapidityOld(ip)-iflagForRapidity(ip)).eq.1) then
                       iflagForRapidity(ip) = 0
                    else
                       iflagForRapidity(ip) = 1
                    endif
                    
                    iflagForRapidityOld(ip)=iflagForRapidityOld(ip)*iflagForRapidityNext(ip)
                 enddo
                 
              enddo
              iflagForRapidity(1:nphi) = iflagForRapidityNext(1:nphi)
              
           enddo
           


           
           ! step 5 for -side 
           iflagForRapidity(1:nphi)=0
           iflagForRapidity(1:nphi)=iflagForRapidity0(1:nphi)
           
           do ith=ntheta/2,1,-1
              iflagForRapidityOld(1:nphi) = iflagForRapidity(1:nphi)
              do searchiteration = 1, int(nphi/nMinLengthFor0)+1
                 
                 if(minval(iflagForRapidity(1:nphi)).eq.1) exit
                 !print *,    searchiteration, int(nphi/nMinlengthFor0)+1
                 do ip = 1,nphi
                     if(iflagForRapidity(ip)) then
                       !ker(ip,ith,2,0)=1.d0       
                       !cycle ! don't calculate if 1
                    else
                       k=k+1
                       xlat=90.d0-theta(ip,ith)
                       xlon=phi(ip,ith)
                       call calculateKernel
                       ker(ip,ith,:,:)=tmpker(:,:)
                    endif
                 enddo
                 
                 
                 
                 ! step 2 and 3 in step5
                 iflagForRapidity(1:nphi) = 1
                 do ift=0,nfilter 
                    do ip = 1,nphi
                       do iitype = 1,nntype
                          itype = idecidetype(iitype)
                          if(abs(ker(ip,ith,itype,ift)).gt.minici(iitype,ift)) then
                             iflagForRapidity(ip) = iflagForRapidity(ip)*0
                          endif
                       enddo
                    enddo
                 enddo
                 
                 ! step 4 in step5
                 
                 
                 
                 modiflag(1:nphi-nMinLengthFor0+1)=1
                 do ip=1,nphi-nMinLengthFor0+1
                    do iiphi=ip,ip+nMinLengthFor0-1
                       modiflag(ip) = modiflag(ip)*iflagForRapidity(iiphi)
                    enddo
                    ! iflagForRapidity(ip:ip+nMinLengthFor0-1)=tmpiflagForRapidity
                 enddo
                 
                 
                 
                 iflagForRapidity(1:nphi)=1
                 do ip=1,nphi-nMinLengthFor0+1
                    if(modiflag(ip).eq.0) then
                       iflagForRapidity(ip:ip+nMinlengthFor0-1)=modiflag(ip)    
                    endif
                 enddo
                 
                 
                 
                 
                 iflagForRapidityNext(1:nphi) = iflagForRapidity(1:nphi)                   
                 
                 do ip=1,nphi
                    if((iflagForRapidityOld(ip)-iflagForRapidity(ip)).eq.1) then
                       iflagForRapidity(ip) = 0
                    else
                       iflagForRapidity(ip) = 1
                    endif
                    
                    iflagForRapidityOld(ip)=iflagForRapidityOld(ip)*iflagForRapidityNext(ip)
                 enddo
                 
              enddo
              iflagForRapidity(1:nphi) = iflagForRapidityNext(1:nphi)
              
           enddo


        endif ! for the fast kernel calculation method
     endif ! if-line for parallelisation

     ! write kernels for each depth

     write(tmpchar,'(I7)') int(rx*1.d3)
     do j=1,7
        if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
     enddo
     kerfile = trim(parentDir)//"/tmp/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)
     
     open(1,file=kerfile,status='unknown',form='unformatted', &
          access = 'direct', recl=kind(0e0)*nphi*ntheta*(nktype+1)*(nfilter+1))    
     write(1,rec=1) ker(1:nphi,1:ntheta,0:nktype,0:nfilter)
     close(1)
  enddo! ir-loop termine
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! write final reordered kernel output
  if(my_rank.eq.0) then   
     deallocate(ker)

     if(trim(paramWRT).eq.'vRSGT') then
        allocate(totalker(nr,nphi,ntheta,1:num_h3,0:nfilter))
     elseif(trim(paramWRT).eq.'vTSGT') then
        allocate(totalker(nr,nphi,ntheta,1:num_h4,0:nfilter))
     elseif((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
        allocate(totalker(nr,nphi,ntheta,0:nkvtype,0:nfilter))
     else
        allocate(totalker(nr,nphi,ntheta,0:nktype,0:nfilter))
     endif
     totalker = 0.e0
     do ir=1,nr
        write(tmpchar,'(I7)') int(r(ir)*1.d3)
        do j=1,7
           if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
        enddo

        kerfile = trim(parentDir)//"/tmp/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)
        
        open(1,file=kerfile,status='unknown',form='unformatted', &
             access = 'direct', recl=kind(0e0)*nphi*ntheta*(nktype+1)*(nfilter+1))    
        read(1,rec=1) totalker(ir,1:nphi,1:ntheta,0:nktype,0:nfilter)
        close(1) 
     enddo


     do ift = 0,nfilter
        kertotalfile = trim(parentDir)//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(freqid(ift))//trim(".kernel")
        open(1,file=kertotalfile,status='unknown',form='unformatted',access='sequential')
        if(compo.eq.'Z') then
           kc=1
       elseif(compo.eq.'R') then
	  kc=2
	elseif(compo.eq.'T') then
	  kc=3
       else
          kc=-1
       endif
       idum=0
       fdum=0.e0

       write(1) totalker(:,:,:,:,ift)
       close(1) 
     enddo




     if(trim(paramWRT).eq.'vRSGT') then
     
        do jt=1,number_of_snapshots
           totalker=0.e0
           do ir=1,nr
              do ith=1,ntheta
                       
                 
                 write(tmpchar,'(I7,".",I7)') ir,ith
                 do j=1,15
                    if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
                 enddo
                 kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                      //"."//"timemarching"
                 open(1,file=kerfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*nphi*(num_h3)*(1+nfilter)*number_of_snapshots)
                 read(1,rec=1) videoker(1:nphi,1:num_h3,0:nfilter,1:number_of_snapshots)
                 close(1) 
                 totalker(ir,1:nphi,ith,1:num_h3,0:nfilter)=videoker(1:nphi,1:num_h3,0:nfilter,jt)
              enddo
           enddo

                 
           do ift = 0,nfilter
              write(tmpchar,'(I7)') jt
              do j=1,7
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kertotalfile = trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"." &
                   //trim(freqid(ift))//"."//trim(tmpchar)//trim(".video")
              open(1,file=kertotalfile,status='unknown',form='unformatted',access='sequential')
              if(compo.eq.'Z') then
                 kc=1
              elseif(compo.eq.'R') then
                 kc=2
              elseif(compo.eq.'T') then
                 kc=3
              else
                 kc=-1
              endif
              idum=0
              fdum=0.e0
              
              write(1) totalker(1:nr,1:nphi,1:nth,1:num_h3,ift)
              close(1) 
              
           enddo
           

        enddo
     endif


      if(trim(paramWRT).eq.'vTSGT') then
     
        do jt=1,number_of_snapshots
           totalker=0.e0
           do ir=1,nr
              do ith=1,ntheta
                       
                 
                 write(tmpchar,'(I7,".",I7)') ir,ith
                 do j=1,15
                    if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
                 enddo
                 kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                      //"."//"timemarching"
                 open(1,file=kerfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*nphi*(num_h3)*(1+nfilter)*number_of_snapshots)
                 read(1,rec=1) videoker(1:nphi,1:num_h4,0:nfilter,1:number_of_snapshots)
                 close(1) 
                 totalker(ir,1:nphi,ith,1:num_h4,0:nfilter)=videoker(1:nphi,1:num_h4,0:nfilter,jt)
              enddo
           enddo

                 
           do ift = 0,nfilter
              write(tmpchar,'(I7)') jt
              do j=1,7
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kertotalfile = trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"." &
                   //trim(freqid(ift))//"."//trim(tmpchar)//trim(".video")
              open(1,file=kertotalfile,status='unknown',form='unformatted',access='sequential')
              if(compo.eq.'Z') then
                 kc=1
              elseif(compo.eq.'R') then
                 kc=2
              elseif(compo.eq.'T') then
                 kc=3
              else
                 kc=-1
              endif
              idum=0
              fdum=0.e0
              
              write(1) totalker(1:nr,1:nphi,1:nth,1:num_h4,ift)
              close(1) 
              
           enddo
           

        enddo
     endif
     
     if((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
        
      
        
        do jt=1,number_of_snapshots
           totalker=0.e0
           do ir=1,nr
              do ith=1,ntheta
                 write(tmpchar,'(I7,".",I7)') ir,ith
                 do j=1,15
                    if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
                 enddo
                 kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
                      //"."//"timemarching"
                 open(1,file=kerfile,status='old',form='unformatted',access='direct',recl=kind(0e0)*nphi*(1+nkvtype)*(1+nfilter)*number_of_snapshots)
                 read(1,rec=1) videoker(1:nphi,0:nkvtype,0:nfilter,1:number_of_snapshots)
                 close(1)
                 totalker(ir,1:nphi,ith,0:nkvtype,0:nfilter)=videoker(1:nphi,0:nkvtype,0:nfilter,jt)

                 ! something that we can try :
                 ! open(1,file=kerfile,status='old',form='unformatted',access='direct',recl=kind(0e0)*nphi*(1+nkvtype)*(1+nfilter))
                 ! read(1,rec=jt) totalker(ir,1:nphi,ith,0:nkvtype,0:nfilter)
                
              enddo
           enddo
           
           
           do ift = 0,nfilter
              write(tmpchar,'(I7)') jt
              do j=1,7
                 if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
              enddo
              kertotalfile = trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"." &
                   //trim(freqid(ift))//"."//trim(tmpchar)//trim(".video")
              open(1,file=kertotalfile,status='unknown',form='unformatted',access='sequential')
              if(compo.eq.'Z') then
                 kc=1
              elseif(compo.eq.'R') then
                 kc=2
              elseif(compo.eq.'T') then
                 kc=3
              else
                 kc=-1
              endif
              idum=0
              fdum=0.e0
              
              write(1) totalker(1:nr,1:nphi,1:nth,0:nkvtype,ift)
              close(1) 
              
           enddo
           

        enddo
                 
     endif

     
     
     ift=nfilter

     infofile = trim(parentDir)//trim(stationName)//"."//trim(eventName)//"."//&
                trim(phase)//"."//trim(compo)//trim(".info")
     open(1,file=infofile,status='unknown',form='unformatted',access='sequential')
     write(1) nr,nphi,ntheta,nktype,nfilter,iWindowEnd-iWindowStart+1,number_of_snapshots
     write(1) real(t(nt1(ift))),real(t(nt2(ift)))
     write(1) (real(mt(i)),i=1,6)
     write(1) real(r0D)
     write(1) real(u)
     close(1)
     
     print *, nr, nphi, ntheta,nktype

     gridfile = trim(parentDir)//trim(stationName)//"."//trim(eventName)//"."//&
                trim(phase)//"."//trim(compo)//trim(".grid")
     open(1,file=gridfile,status='unknown',form='unformatted',access='sequential')
     write(1) nr, nphi, ntheta, nktype
     write(1) real(r)
     write(1) real(phitheta)
     write(1) real(thetaphi)
     close(1)


     list = trim(parentDir)//"/log/calLog"//"."// &
          trim(stationName)//"."//trim(eventName)//"."//trim(compo)//"."//trim(paramWRT)//".log"   
     open(1,file =list, status = 'old',access='append', form = 'formatted')
     call date_and_time(datex,timex)
     write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
          '    Finishing date and time:                     ', &
          datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
          timex(1:2),':',timex(3:4),':',timex(5:8)   
     close (1)
     
  endif

  call MPI_FINALIZE(ierr)
  stop

end program KernelMaker
  
