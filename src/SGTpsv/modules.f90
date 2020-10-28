module parameters  


!-------------------------<< input matrix >>----------------------------------
  !include 'mpif.h'  
  character(120) :: outputDir, psvmodel, modelname,DSMconfFile 
  character(120) :: list,list1
  character(40) :: datex,timex
  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
  real(kind(0d0)) :: re,ratc,ratl
  integer :: maxlmax

  real(kind(0d0)) :: tlen 
  real(kind(0d0)) :: rmin_, rmax_, rdelta_ 
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: thetamin, thetamax, thetadelta
  real(kind(0d0)), allocatable :: r_(:),r0(:),theta(:)
  real(kind(0d0)), allocatable :: rrsta(:,:)
  real(kind(0d0)), allocatable :: lambda(:),qkp(:) !! only for stack points and only used for fluid regions

  logical, allocatable :: log_solid_liquid(:) ! true for solid and false for liquid
   
  integer, allocatable :: iista(:,:)
  integer :: r_n,r0_n,ciista, ir_,ir0,imt,icomp,itheta, theta_n
  integer :: intir0  
  character(120) :: coutfile
  integer :: imin, imax
  integer :: rsgtswitch, tsgtswitch, synnswitch, psgtswitch


  
  ! ---------------------------<< variables >>---------------------------
  ! variable for the trial function
  integer:: nnlayer,nlay
  integer, allocatable :: nlayer(:)
  integer:: nslay,nllay
  integer:: inlayer,jnlayer,jnslay,jnllay
  integer:: l,m
  real(kind(0d0)),allocatable:: ra(:)
  ! variable for the structure
  integer:: nzone,isl,ill,nsl,nll
  integer,allocatable:: iphase(:)
  integer::ndc,vnp
  real(kind(0d0)):: rmin,rmax
  real(kind(0d0)),allocatable:: vrmin(:),vrmax(:),rrho(:,:),vpv(:,:),vph(:,:),vsv(:,:),vsh(:,:),eta(:,:),qmu(:),qkappa(:)
  real(kind(0d0)),allocatable::vra(:),rho(:),kappa(:) 
  real(kind(0d0)),allocatable::ecKx(:) !3*Kx=3A-4N
  real(kind(0d0)),allocatable::ecKy(:) !3*Ky=3F+2N
  real(kind(0d0)),allocatable::ecKz(:) !3*Kz=2F+C
  real(kind(0d0)),allocatable::mu(:),ecL(:),ecN(:),rhoinv(:),kappainv(:)
  complex(kind(0d0)),allocatable:: coef1(:),coef2(:),coef(:)
  ! variable for the periodic range
  
  real(kind(0d0)):: omega,omegai
  ! variable for the source
  integer:: spn,ns
  real(kind(0d0)):: mt(3,3),spo
  real(kind(0d0)):: ecC0,ecF0,ecL0
  complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4)

  ! variable for the matrix elements
  complex(kind(0d0)),allocatable:: a0(:,:),a1(:,:),a2(:,:), a(:,:), c(:,:), ctmp(:,:)
  real(kind(0d0)), allocatable :: t(:)
  real(kind(0d0)), allocatable :: h1x(:), h1y(:), h1z(:), h2L(:), h2N(:), h3ax(:), h3ay(:), h3az(:), h4aL(:), h4aN(:), h5ax(:), h5ay(:), h5az(:), h6aL(:), h6aN(:), h3x(:), h3y(:), h3z(:), h4L(:), h4N(:), h5x(:), h5y(:), h5z(:), h6L(:), h6N(:), h7x(:), h7y(:), h7z(:), h8L(:), h8N(:), h3mx(:,:), h3my(:,:), h3mz(:,:), h5mx(:,:), h5my(:,:), h5mz(:,:), h4m1L(:,:), h4m1N(:,:), h4m2L(:,:), h4m2N(:,:), h6m1L(:,:), h6m1N(:,:), h6m2L(:,:), h6m2N(:,:)
  real(kind(0d0)),allocatable:: p1(:),p2(:),p3(:)
  complex(kind(0d0)),allocatable:: g0(:)
  complex(kind(0d0)),allocatable:: d0(:)
  complex(kind(0d0)):: g0tmp(2),g0dertmp(2) ! forward
  ! variable for the stack point
  integer,allocatable:: isp(:),issp(:),ilsp(:),jssp(:),jsp(:), ksp(:),lsp(:)
  integer::isdr,jsdr,ildr,cista,cksta
  ! variables for the output stack point
  integer,allocatable:: istazone(:)
  integer,allocatable:: ksta(:)   ! output stack point for g
  integer,allocatable:: jsta(:)   ! output stack point for d
  
 
  ! variables for the gridding
  integer,allocatable:: jjdr(:),kkdr(:)
  integer:: jdr,kdr
  real(kind(0d0)),allocatable:: vmin(:),gridpar(:),dzpar(:)
  ! variables for l cut off
  integer:: kc,lsuf,sufzone,ismall,llog
  real(kind(0d0)):: maxamp
  ! variables for the numerical integration
  complex(kind(0d0)):: anum(4,4,10),bnum(4,4,10)
 
  ! other variables
  integer:: i,j,nn,ier,itmp,jtmp,mtmp,kkdr0,nn0,ig2
  integer:: ll(12),lli(12),llj(12)
  real(kind(0d0)):: eps,l2,lsq
  real(kind(0d0)),allocatable:: work(:)
  complex(kind(0d0)), allocatable ::z(:), w(:),cwork(:)
  
  !-----------------------------------------------------------------------
  !complex(kind(0d0)), allocatable :: dvec(:,:,:,:),dvecdt(:,:,:,:),dvecdp(:,:,:,:)
  complex(kind(0d0)), allocatable :: dvec0(:,:,:),dvecdt0(:,:,:),dvecdp0(:,:,:)
  complex(kind(0d0)), allocatable :: tsgt(:,:,:,:),rsgt(:,:,:),synn(:,:) ,psgt(:,:,:,:)
  complex(kind(0e0)), allocatable :: tsgtsngl(:,:), rsgtsngl(:,:),synnsngl(:,:),psgtsngl(:,:)
  real(kind(0d0)), allocatable :: plm(:,:,:)
  complex(kind(0d0)) :: rdvec(1:3,-2:2)
  complex(kind(0d0))::u(1:3),udr(1:3),udt(1:3),udp(1:3),uder(1:3,1:3)



  data eps /-1.d0/

  ! for Pinv : num_tsgt = 4, num_rsgt = 1
  
  !integer, parameter :: num_tsgt = 4
  !integer, parameter :: num_rsgt = 1
  integer, parameter :: num_tsgt = 20
  integer, parameter :: num_rsgt = 10
  integer, parameter :: num_synn = 10



  
  !--------------------------------------------------------------------------
  ! for MPI

  integer :: nproc,my_rank,ierr
    
  !--------------------------------------------------------------------------
  

  
end module parameters
