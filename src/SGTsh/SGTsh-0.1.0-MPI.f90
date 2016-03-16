program  SGTsh


  !-----------------------------------------------------------------------
  !      SGTsh
  !
  !  
  !
  !
  !  calculation de la fonction de Green pour SH
  !       
  !                                               2002.10.KAWAI Kenji
  !                                               2009.6. FUJI Nobuaki
  !                                               2010.9. FUJI Nobuaki
  !                                               2012.3. FUJI Nobuaki
  !                                               
  !
  !
  !                 
  !
  !-----------------------------------------------------------------------
  use mpi
  implicit none
  
  !-------------------------<< input matrix >>----------------------------------
  
  character(120) :: outputDir, psvmodel,modelname,DSMconfFile
  character(120) :: list,list1,tmpfile
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)) ::re,ratc,ratl
  integer :: maxlmax
  
  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: rmin_, rmax_, rdelta_ 
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: thetamin, thetamax, thetadelta
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ciista, ir_,ir0,imt,icomp,idepth,itheta, theta_n
  
  character(120) :: coutfile
  integer :: imin, imax
  integer :: rsgtswitch, tsgtswitch, synnswitch

  
  !---------------------------<< variables >>---------------------------
  ! variable for the trial function
  
  integer :: nzone
  integer :: iimax,ii
  integer :: i ,j, ier,jj
  real(kind(0d0)) :: dummy
  real(kind(0d0)), allocatable :: vrmin(:), vrmax(:)
  real(kind(0d0)), allocatable :: rrho(:,:), vsv(:,:), vsh(:,:), qmu(:)
  real(kind(0d0)), allocatable :: vra (:), rho(:), ecL(:), ecN(:)
  real(kind(0d0)), allocatable :: gvra(:,:), grho(:,:), gecL(:,:), gecN(:,:),gra(:,:)
  complex(kind(0d0)), allocatable :: coef(:),cwork(:)
  real(kind(0d0)) :: rmin, rmax
  real(kind(0d0)), allocatable :: vmin(:), gridpar(:), dzpar(:)
  complex(kind(0d0)), allocatable :: tmpc(:)
  real(kind(0d0)) :: maxamp
  real(kind(0d0)) :: omegai
  real(kind(0d0)), allocatable :: ra(:)
  integer :: nnlayer, vnp,nn
  integer, allocatable :: nlayer(:), iphase(:)
  integer :: ioutercore
  
  ! variables pour des points stackes 
  integer, allocatable :: isp(:),jsp(:),ins(:)
  
  ! variables pour la source
  
  integer, allocatable :: spn(:),ns(:)
  real(kind(0d0)) :: mt(3,3),lsq
  real(kind(0d0)), allocatable :: mu0(:),spo(:)
  
  !-----------------------------------------------------------------------
  ! variables pour des elements de matrice 
  complex(kind(0d0)), allocatable :: a0(:,:), a2(:,:), a(:,:),dr(:),z(:)
  real(kind(0d0)), allocatable :: t(:), h1(:), h2(:), h3(:), h4(:), work(:)
  real(kind(0d0)), allocatable :: gt(:,:),gh1(:,:),gh2(:,:),gh3(:,:),gh4(:,:)
  complex(kind(0d0)),allocatable :: aa(:,:), ga(:,:),ga2(:,:,:),gdr(:,:)
  complex(kind(0d0)), allocatable :: g0(:)
  complex(kind(0d0)) :: g0tmp, g0dertmp
  ! la frequence
  real(kind(0d0)) :: omega
  integer :: lsuf
  
  ! des autres 
  integer :: lda 
  integer :: kc, ismall, llog,m, l,ig2
  real(kind(0d0)) :: eps 
  
  !-----------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: bvec(:,:,:),bvecdt(:,:,:),bvecdp(:,:,:)
  complex(kind(0d0)), allocatable :: tsgt(:,:,:,:),rsgt(:,:,:),synn(:,:)
  complex(kind(0e0)), allocatable :: tsgtsngl(:,:), rsgtsngl(:,:),synnsngl(:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  complex(kind(0d0)) :: rvec(1:3,-2:2)
  complex(kind(0d0))::u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)
  
  
  data lda / 2 /
  data eps / -1.d0 /
  
  integer, parameter :: num_tsgt = 10
  integer, parameter :: num_rsgt = 5
  integer, parameter :: num_synn = 5
  
  
  !-----------------------------------------------------------------------
  ! for MPI
  
  integer :: nproc, my_rank, ierr  
  
  !-----------------------------------------------------------------------


  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  
  if(my_rank.eq.0) then
     call pinput(DSMconfFile,outputDir,psvmodel,modelname,tlen,rmin_,rmax_,rdelta_,r0min,r0max,r0delta,thetamin,thetamax,thetadelta,imin,imax,rsgtswitch,tsgtswitch,synnswitch)
     call readDSMconf(DSMconfFile,re,ratc,ratl,omegai,maxlmax)
     tmpfile = 'tmpworkingfile_for_psvmodel'
     call readpsvmodel(psvmodel,tmpfile)
     psvmodel = 'tmpworkingfile_for_psvmodel'
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     close(20)
  endif
  
  psvmodel = 'tmpworkingfile_for_psvmodel'

  ! exporting DSM parameters
  call MPI_BCAST(re,  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ratc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ratl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(omegai,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(maxlmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  call MPI_BCAST(outputDir,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(psvmodel,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(modelname,120,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(imin,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(imax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rsgtswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(tsgtswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(synnswitch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  



  if(my_rank.eq.0) then
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     allocate(vrmin(1:nzone))
     allocate(vrmax(1:nzone))
     allocate(rrho(1:4,1:nzone))
     allocate(vsv(1:4,1:nzone))
     allocate(vsh(1:4,1:nzone))
     allocate(qmu(1:nzone))
     allocate(vmin(1:nzone))
     allocate(gridpar(1:nzone))
     allocate(dzpar(1:nzone))
     allocate(nlayer(1:nzone))
     allocate(iphase(1:nzone))
     allocate(isp(1:nzone))
     allocate(jsp(1:nzone))
     allocate(coef(1:nzone))
     do i = 1, nzone
        read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), dummy, dummy, dummy, dummy, qmu(i), dummy
        if ((vsv(1,i).eq.0.d0).and.(vsv(2,i).eq.0.d0).and.(vsv(3,i).eq.0.d0).and.(vsv(4,i).eq.0.d0)) then
           iphase(i) = 2
           ioutercore = i
        else
           iphase(i) = 1
        endif
     enddo
     close(20)

     
     
     ! CAUTION: this program can only calculate for solid media (SH) for this moment
     
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     nzone = nzone - ioutercore
     deallocate(vrmin,vrmax,rrho,vsv,vsh,qmu,vmin,gridpar,dzpar,nlayer,iphase,isp,jsp,coef)
     close(20)
  endif

  
  call MPI_BCAST(nzone,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  allocate(vrmin(1:nzone))
  allocate(vrmax(1:nzone))
  allocate(rrho(1:4,1:nzone))
  allocate(vsv(1:4,1:nzone))
  allocate(vsh(1:4,1:nzone))
  allocate(qmu(1:nzone))
  allocate(vmin(1:nzone))
  allocate(gridpar(1:nzone))
  allocate(dzpar(1:nzone))
  allocate(nlayer(1:nzone))
  allocate(iphase(1:nzone))
  allocate(isp(1:nzone))
  allocate(jsp(1:nzone))
  allocate(coef(1:nzone))


  if(my_rank.eq.0) then
     
     open(20, file = psvmodel, status = 'old', action='read', position='rewind')
     read(20,*) nzone
     nzone = nzone - ioutercore
     
     do i = 1,ioutercore 
        read (20, *) dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy
     enddo
     do i = 1, nzone
        read (20, *) vrmin(i), vrmax(i), rrho(1,i), rrho(2,i), rrho(3,i), rrho(4,i), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, vsv(1,i), vsv(2,i), vsv(3,i), vsv(4,i), vsh(1,i), vsh(2,i), vsh(3,i), vsh(4,i), dummy, dummy, dummy, dummy, qmu(i), dummy
     enddo
     close(20)
     
  endif
  
  call MPI_BCAST(vrmin,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vrmax,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(rrho,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vsv,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(vsh,4*nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(qmu,nzone,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
  rmin = vrmin(1)
  rmax = vrmax(nzone)
  omegai = - dlog(omegai) / tlen
  
  ! depths for stocking the Green function
  
  if(my_rank.eq.0) then
     r_n =  int((rmax_-rmin_)/rdelta_)+1
  endif

  call MPI_BCAST(r_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)



  allocate(r_(1:r_n))
  allocate(rrsta(1:3,1:r_n))
  allocate(iista(1:3,1:r_n))
  

  if(my_rank.eq.0) then
     do i = 1, r_n
        r_(i) = rmin_ + dble(i-1)*rdelta_
     enddo
  endif
  
 
  call MPI_BCAST(r_,r_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


  ! GCARCs

  if(my_rank.eq.0) then
     theta_n = int((thetamax-thetamin)/thetadelta)+1
  endif

  call MPI_BCAST(theta_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(theta(1:theta_n))
  allocate(bvec(1:3,-2:2,1:theta_n))
  allocate(bvecdt(1:3,-2:2,1:theta_n))
  allocate(bvecdp(1:3,-2:2,1:theta_n))
  allocate(plm(1:3,0:3,1:theta_n))

  if(my_rank.eq.0) then   
     do i = 1,theta_n
        theta(i) = (thetamin + dble(i-1)*thetadelta)
     enddo
  endif
  
  call MPI_BCAST(theta,theta_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  

  ! source depths

  if(my_rank.eq.0) then
     r0_n =  int((r0max-r0min)/r0delta)+1
  endif


  call MPI_BCAST(r0_n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  allocate(r0(1:r0_n))
  allocate(spo(1:r0_n))
  allocate(spn(1:r0_n))
  allocate(ns(1:r0_n))
  allocate(mu0(1:r0_n))
  allocate(ins(1:r0_n))
  allocate(gra(1:3,1:r0_n)) 
  allocate(gvra(1:3,1:r0_n)) 
  allocate(grho(1:3,1:r0_n))
  allocate(gecL(1:3,1:r0_n))
  allocate(gecN(1:3,1:r0_n)) 
  allocate(gt(1:8,1:r0_n))
  allocate(gh1(1:8,1:r0_n))
  allocate(gh2(1:8,1:r0_n))
  allocate(gh3(1:8,1:r0_n))
  allocate(gh4(1:8,1:r0_n)) 
  allocate(aa(1:4,1:r0_n))
  allocate(ga(1:8,1:r0_n))
  allocate(ga2(1:2,1:3,1:r0_n))
  allocate(gdr(1:3,r0_n))

  


  if(my_rank.eq.0) then
     do i = 1, r0_n
        r0(i) = r0min + dble(i-1)*r0delta
     enddo
     
  endif
  
  call MPI_BCAST(r0,r0_n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ir0 = r0_n
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)


  ! SSGT and TSGT allocation

 
  allocate(tsgt(1:num_tsgt,1:r_n,1:theta_n,1:r0_n))
  allocate(rsgt(1:num_rsgt,1:r_n,1:theta_n))
  allocate(synn(1:num_synn,1:theta_n))

  allocate(tsgtsngl(1:num_tsgt,1:theta_n))
  allocate(rsgtsngl(1:num_rsgt,1:theta_n))
  allocate(synnsngl(1:num_synn,1:theta_n))



  ! computation de nombre et la location des points de grid
  
  iimax =imax

  call calgrid( nzone,vrmin,vrmax,vsv,rmin,rmax,iimax,1,tlen,vmin,gridpar,dzpar )
  call calra(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,re )
  allocate(ra(1:nnlayer+nzone+1))
  call calra2(nnlayer,gridpar,dzpar,nzone,vrmin,vrmax,rmin,rmax,nlayer,ra,re,r_n,r_,rrsta, iista)
  ! computation de points stackes et la location de la source
  call calsp( nzone-1,nlayer,isp,jsp )
  do ir0 = 1,r0_n
     call calspo( nzone-1,vrmax,nnlayer,r0(ir0),rmin,rmax,ra,isp,spo(ir0),spn(ir0) )
     call calgra( isp,ra,r0(ir0),spn(ir0),spo(ir0),gra(1:3,ir0))
  enddo

  ! computation des elements de matrice

  allocate(vra(1:nnlayer+2*nzone+1))
  allocate(rho(1:nnlayer+2*nzone+1))
  allocate(ecL(1:nnlayer+2*nzone+1))
  allocate(ecN(1:nnlayer+2*nzone+1))
  allocate(a0(1:2,1:nnlayer+1))
  allocate(a2(1:2,1:nnlayer+1))
  allocate(a(1:2,1:nnlayer+1))
  allocate(t(1:4*nnlayer))
  allocate(cwork(1:4*nnlayer))
  allocate(h1(1:4*nnlayer))
  allocate(h2(1:4*nnlayer))
  allocate(h3(1:4*nnlayer))
  allocate(h4(1:4*nnlayer))
  allocate(work(1:4*nnlayer))
  allocate(tmpc(nnlayer+1))
  allocate(g0(1:nnlayer+1))
  allocate(dr(1:nnlayer+1))
  allocate(z(1:nnlayer+1))  
  call calstg( nzone,rrho,vsv,vsh,nnlayer,nlayer,ra,rmax,vnp,vra,rho,ecL,ecN)
  do ir0 = 1, r0_n
     call calgstg(nzone,nnlayer,spn(ir0),rrho,vsv,vsh,gra(1:3,ir0),gvra(1:3,ir0),rmax,grho(1:3,ir0),gecL(1:3,ir0),gecN(1:3,ir0),r0(ir0),mu0(ir0)) 
  enddo

  do i= 1, nzone
     call calmatc(nlayer(i),vnp,vra,rho,2,0,0,ra(isp(i)),t(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,2,1,1,ra(isp(i)),h1(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,1,1,0,ra(isp(i)),h2(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecL,0,0,0,ra(isp(i)),h3(jsp(i)),work(jsp(i)))
     call calmatc(nlayer(i),vnp,vra,ecN,0,0,0,ra(isp(i)),h4(jsp(i)),work(jsp(i)))
     call caltl(nlayer(i),vnp,vra,rho,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),t(jsp(i)),work(jsp(i)),t(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecL,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h3(jsp(i)),work(jsp(i)),h3(jsp(i)))
     call calhl(nlayer(i),vnp,vra,ecN,ra(isp(i)),work(jsp(i)))
     call calt(nlayer(i),h4(jsp(i)),work(jsp(i)),h4(jsp(i)))
  enddo
  do ir0 = 1, r0_n
     call calmatc( 2,3,gvra(1:3,ir0),grho(1:3,ir0),2,0,0,gra(1:3,ir0),gt(1:8,ir0), work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),2,1,1,gra(1:3,ir0),gh1(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),1,1,0,gra(1:3,ir0),gh2(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecL(1:3,ir0),0,0,0,gra(1:3,ir0),gh3(1:8,ir0),work )
     call calmatc( 2,3,gvra(1:3,ir0),gecN(1:3,ir0),0,0,0,gra(1:3,ir0),gh4(1:8,ir0),work )
     call caltl(2,3,gvra(1:3,ir0),grho(1:3,ir0),gra(1:3,ir0),work )
     call calt(2,gt(1:8,ir0),work,gt(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecL(1:3,ir0),gra(1:3,ir0),work)
     call calt( 2,gh3(1:8,ir0),work,gh3(1:8,ir0))
     call calhl(2,3,gvra(1:3,ir0),gecN(1:3,ir0),gra(1:3,ir0),work )
     call calt( 2,gh4(1:8,ir0),work,gh4(1:8,ir0))  
  enddo


  !computation de la dislocation
  nn = nnlayer + 1
  do ir0 = 1, r0_n
     ns(ir0) = isp(spn(ir0)) + dint(spo(ir0))
     ins(ir0) = 4 * ns(ir0) - 3
  enddo


  !******************** Computing the displacement *********************


  llog = 0
  
   write(list1, '(I7,".",I7,".",I7)') imin,imax,my_rank
  do j = 1,22
     if(list1(j:j).eq.' ')list1(j:j) = '0'
  enddo
  list1 = trim(outputDir)//"/log/listSH"//"."//trim(modelname)//"."//trim(list1)
  
  open(24, file = list1, status = 'unknown', form = 'formatted')
  write(24,*) 
  close(24)

  


  
  if(my_rank.eq.0) then
     write(list, '(I7,".",I7)') imin,imax
     do j = 1,15
        if(list(j:j).eq.' ')list(j:j) = '0'
     enddo
     list = trim(outputDir)//"/log/calLogSH"//"."//trim(modelname)//"."//trim(list)
     
     
     open(1,file =list, status = 'unknown', form = 'formatted')
     call date_and_time(datex,timex)
     write(1,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') &
          '    Starting date and time:                     ', &
          datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ', &
          timex(1:2),':',timex(3:4),':',timex(5:8)   
     close (1)
  endif
  
  
  
  
  do i = imin, imax ! each frequency        
     
     tsgt = cmplx(0.d0)
     rsgt = cmplx(0.d0)
     synn = cmplx(0.d0)
     
     tsgtsngl = cmplx(0.e0)
     rsgtsngl = cmplx(0.e0)
     synnsngl = cmplx(0.e0)

     
     omega = 2.d0 * pi * dble(i)/tlen

     
     if((i.ne.0).and.((mod(imax-my_rank-i,2*nproc).eq.0).or.(mod(imax+my_rank+1-i,2*nproc).eq.0))) then
        call callsuf(omega,nzone,vrmax,vsv,lsuf)       
        call calcoef( nzone,omega,qmu,coef)
        plm = 0.d0
        a0 = 0.d0
        a2 = 0.d0
        do j = 1, nzone
           call cala0( nlayer(j),omega,omegai,t(jsp(j)), h1(jsp(j)),&
                & h2(jsp(j)), h3(jsp(j)), h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a0( 1,isp(j) ) )
           call cala2( nlayer(j),h4(jsp(j)),coef(j), cwork(jsp(j)) )
           call overlap( nlayer(j),cwork(jsp(j)),a2( 1,isp(j) ) )
        enddo
 

        kc = 1
        ismall = 0
        maxamp = -1.d0
        llog = maxlmax
        
        do l = 0, maxlmax ! l-loop commence
           lsq = dsqrt(dble(l)*dble(l+1))
           do itheta = 1,theta_n
              call calbvecphi0(l,(theta(itheta)/180.d0*pi),plm(1,0,itheta),bvec(1,-2,itheta),bvecdt(1,-2,itheta),bvecdp(1,-2,itheta))              
           enddo
           
           rvec = cmplx(0.d0)
           call calbveczero(l,rvec(1,-2))
           
           if(ismall.gt.20) then
              if(llog.gt.l) then
                 llog =l
              endif
              cycle   !exit
           endif
           
           tmpc = 0.d0
           
           a = 0.d0
           ga2 = 0.d0
           
           call cala( nn,l,lda,a0,a2,a )
         
	   ir0=1
           call calga( 1,omega,omegai,l,t(ins(ir0)),h1(ins(ir0)),h2(ins(ir0)),h3(ins(ir0)),h4(ins(ir0)),coef(spn(ir0)),aa(1:4,ir0))
           call calga( 2,omega,omegai,l,gt(1:8,ir0),gh1(1:8,ir0),gh2(1:8,ir0),gh3(1:8,ir0),gh4(1:8,ir0),coef(spn(ir0)),ga(1:8,ir0))
           call overlap( 2,ga(1:8,ir0),ga2(1:2,1:3,ir0))
           
           
           
           do m = -2, 2 ! m-loop commence       
              if ( ( m.ne.0 ).and.( iabs(m).le.iabs(l) ) ) then
                 
                 
                 ir0=1
                 ig2 = 0
                 do imt = 2,6
                    g0 = cmplx(0.d0)
                    call setmt(imt,mt)  
                    call calg2( l,m,spo(ir0),r0(ir0),mt,mu0(ir0),coef(spn(ir0)),ga(1:8,ir0),aa(1:4,ir0),ga2(1:2,1:3,ir0),gdr(1:3,ir0),g0(isp(spn(ir0))),ig2)
                    
                    
                    if ( (m.eq.-2).or.(m.eq.-l) ) then             
                       if(ig2.eq.1) then
                          call dclisb0( a,nn,1,lda,g0,eps,dr,z,ier)
                          ig2 = ig2+1
                       else
                          call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
                       endif
                    else 
                       call dcsbsub0( a,nn,1,lda,g0,eps,dr,z,ier)
                    endif
                    
                                       
                    if((imt.eq.3).and.(ir0.eq.r0_n)) then
                       call calamp(g0(nn),l,lsuf,maxamp,ismall,ratl)
                    endif
                       

                    
                    if(synnswitch.eq.1) then
                       do itheta = 1, theta_n
                          u = cmplx(0.d0)
                          call calu(g0(nn),lsq,bvec(1:3,m,itheta),u(1:3))
                          call utosynnSH(imt,u(1:3),synn(1:num_synn,itheta))
                       enddo
                    endif
                    

                    do ir_= 1,r_n                    
                       g0tmp = 0.d0
                       g0dertmp = 0.d0
                       call interpolate(1,0,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0tmp)
                       call interpolate(1,1,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0dertmp)
                       
                       
                       
                       do itheta = 1, theta_n
                          u = cmplx(0.d0)
                          udr = cmplx(0.d0)
                          udt = cmplx(0.d0)
                          udp = cmplx(0.d0)
                          uder = cmplx(0.d0)
                          call calu(g0tmp,lsq,bvec(1:3,m,itheta),u(1:3))
                          call calu(g0dertmp,lsq,bvec(1:3,m,itheta),udr(1:3))
                          call calu(g0tmp,lsq,bvecdt(1:3,m,itheta),udt(1:3))
                          call calu(g0tmp,lsq,bvecdp(1:3,m,itheta),udp(1:3))
                          call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3),r_(ir_),theta(itheta)/180.d0*pi)
                          call udertotsgtSH(imt,uder,tsgt(1:num_tsgt,ir_,itheta,ir0))
                       enddo
                    enddo ! ir_-loop termine
                 enddo ! imt-loop termine
                 !enddo !ir0-loop termine     
                 
                 ! back-propagated wavefield
                 if(iabs(m).eq.1) then     
                    do icomp = 2,3
                       g0 = cmplx(0.d0)
                       g0(nn) = -conjg(rvec(icomp,m))/cmplx(lsq)
                       call dcsbsub0(a,nn,1,lda,g0,eps,dr,z,ier)
                       !print *, g0
                       do ir_= 1,r_n      
                          
                          g0tmp = 0.d0
                          g0dertmp = 0.d0
                          
                          call interpolate(1,0,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0tmp)
                          call interpolate(1,1,r_(ir_),rrsta(1:3,ir_),g0(iista(1:3,ir_)),g0dertmp)
                          !print *, g0tmp
                          do itheta = 1, theta_n
                             u = cmplx(0.d0)
                             udr = cmplx(0.d0)
                             udt = cmplx(0.d0)
                             udp = cmplx(0.d0)
                             uder = cmplx(0.d0)
                             call calu(g0tmp,lsq,bvec(1:3,m,itheta),u(1:3))
                             call calu(g0dertmp,lsq,bvec(1:3,m,itheta),udr(1:3))
                             call calu(g0tmp,lsq,bvecdt(1:3,m,itheta),udt(1:3))
                             call calu(g0tmp,lsq,bvecdp(1:3,m,itheta),udp(1:3))
                             
                             call locallyCartesianDerivatives(u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3),r_(ir_),theta(itheta)/180.d0*pi)
                             
                             call udertorsgtSH(icomp,uder(1:3,1:3),rsgt(1:num_rsgt,ir_,itheta))
                             
                          
                          enddo
                       enddo ! ir_-loop termine
                    enddo
                 endif
                 ! back propagated waves finishes
                 
              endif
           enddo ! m-loop termine              
        enddo !l-loop termine         
        
        
        open(24,file =list1, status = 'old',access='append', form = 'formatted')
        write(24,*) i, dble(i)/tlen, llog-1     
        close(24)
        
        
        
        
        
     
        do ir_ = 1,r_n
           
           ir0=1
           write(coutfile, '(I7,".",I7,".",I7,".TSGT_SH")') int(r0(ir0)*1.d3),int(r_(ir_)*1.d3),i
           do j = 1,29
              if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
           enddo
           
           coutfile = trim(modelname)//"."//coutfile
           coutfile = trim(outputDir)//"/TSGT/"//coutfile
           open(1,file=coutfile,status='unknown',form='unformatted', &
                access = 'direct', recl=2*num_tsgt*kind(0e0)*theta_n)
           tsgtsngl(1:num_tsgt,1:theta_n) = tsgt(1:num_tsgt,ir_,1:theta_n,ir0)
           write(1,rec=1)tsgtsngl(1:num_tsgt,1:theta_n)
           
           close(1)          
           
           if(rsgtswitch) then
              write(coutfile, '(I7,".",I7,".RSGT_SH")') int(r_(ir_)*1.d3),i
              do j = 1,21
                 if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
              enddo
              coutfile = trim(modelname)//"."//coutfile
              coutfile = trim(outputDir)//"/RSGT/"//coutfile
              open(1,file=coutfile,status='unknown',form='unformatted', &
                   access = 'direct', recl=2*num_rsgt*kind(0e0)*theta_n)
              rsgtsngl(1:num_rsgt,1:theta_n) = rsgt(1:num_rsgt,ir_,1:theta_n)
              write(1,rec=1)rsgtsngl(1:num_rsgt,1:theta_n)
              close(1)                    
           endif
        enddo
        
        ir0 = 1
        if(synnswitch) then
           write(coutfile, '(I7,".",I7,".SYNN_SH") ') int(r0(ir0)*1.d3),i
           do j = 1,21
              if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
           enddo
           coutfile = trim(modelname)//"."//coutfile
           coutfile = trim(outputDir)//"/RSGT/"//coutfile
           open(1,file=coutfile,status='unknown',form='unformatted', &
                access = 'direct', recl=2*num_synn*kind(0e0)*theta_n)
           synnsngl(1:num_synn,1:theta_n) = synn(1:num_synn,1:theta_n)
           write(1,rec=1)synnsngl(1:num_synn,1:theta_n)
           close(1)
        endif
     endif
  enddo ! omega-loop termine   
  
  close(24)
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  
  if(my_rank.eq.0) then
     
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

end program SGTsh
