module parameters



  
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
  integer :: intir0
  character(120) :: coutfile
  integer :: imin, imax
  integer :: rsgtswitch, tsgtswitch, synnswitch,psgtswitch

  
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
  complex(kind(0d0)), allocatable :: tsgt(:,:,:,:),rsgt(:,:,:),synn(:,:),psgt(:,:,:,:)
  complex(kind(0e0)), allocatable :: tsgtsngl(:,:), rsgtsngl(:,:),synnsngl(:,:),psgtsngl(:,:)
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


end module parameters
