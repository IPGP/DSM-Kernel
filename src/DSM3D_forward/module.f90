
!  Modules for DSMpsv3D
!                                      2016.03. FUJI Nobuaki
!
!           Copyright Institut de Physique du Globe de Paris


! module indices
! module constants
! module variables
! module receiversphericals
! module scatterersphericals
! module outputs
! module paramsmpi

! some of the features required for Frecht derivatives are ommitted (u,udr,udt,udp) but 
! one can look at SGTpsv/sh when needed to extend this forward modeling program for partial derivatives

module indices

  implicit none
  
  ! Indices fixed during some long time or throughout the whole calculation in an individual processor
  ! and some useful dummy parameters

  
  ! Strings for time calculations
  character(40) :: datex,timex


  ! dummy parameters
  integer :: idummy
  character(120) :: coutfile
  
  !! indices related to constants module

  ! indices for double couple and single force component
  integer :: imt, icomp
  
  ! indices for receivers
  integer :: ir_, i_angle, i_station

  ! indices for scatteres
  integer :: jr_, j_angle, j_scatter

  !! indices related to variables module

  integer:: inlayer,jnlayer,jnslay,jnllay
  integer:: l,m
  integer:: isl, ill
  
  ! indices for stack points
  
  integer::isdr,jsdr,ildr,cista,cksta
 
  ! indices for gridding
  
  integer:: jdr,kdr

  ! other indices

  integer:: i,j,nn,ier,itmp,jtmp,mtmp,kkdr0,nn0,ig2

end module indices

module constants
  
  implicit none

  ! Constants module for DSMpsv3D
  ! Most of them are defined in the .inf file for DSMpsv3D

  real(kind(0d0)), parameter :: pi=3.1415926535897932d0
    
  ! Directories and files to be defined in inf file
  
  character(120) :: outputDir, psvmodel, modelname,DSMconfFile 
  character(120) :: list,list1

  ! DSM constants for error control

  real(kind(0d0)) :: re,ratc,ratl
  integer :: maxlmax
  
  ! DSM calculation constants
  
  real(kind(0d0)) :: tlen ! Window length (s)
  integer :: imin, imax ! frequency min/max


  ! Source depth

  real(kind(0d0)) :: r0

  ! Source components (double couple and single force) only for sources in the solid part
  real(kind(0d0)) :: mt(3,3)
  real(kind(0d0)) :: singleforce(3)
  
  integer:: spn,ns
  real(kind(0d0)):: spo
  real(kind(0d0)):: ecC0,ecF0,ecL0
  complex(kind(0d0)):: ya(4),yb(4),yc(4),yd(4)

  
!!! RECEIVER LOCATIONS


  ! Receiver files ! NF it is better to include a preprocessing subroutine to create those files 
  character(120) :: recradiusfile ! with i_r_index(=idummy), r_(ir_index(1:r_n))
  character(120) :: recanglefile ! with i_angle_index(=idummy),theta(1:angle_n),phi(1:angle_n)
  character(120) :: reclocfile ! with idummy, i_r_index(1:station_n),i_angle_index(1:station_n)

  ! Receiver depths and locations

  integer :: r_n ! number of vertical stacking points for receivers
  integer :: angle_n ! number of lateral stacking points for receivers
  integer :: station_n ! number of receivers


  real(kind(0d0)), dimension(:), allocatable :: r_ ! with r_n elements
  real(kind(0d0)), dimension(:), allocatable :: theta, phi! with theta_n elements
  
  integer, dimension(:), allocatable :: i_r_index ! with r_n elements
  integer, dimension(:), allocatable :: i_angle_index ! with angle_n elements

  ! Interpolation matrices for vertical stacking points 
 
  real(kind(0d0)), dimension(:,:), allocatable :: rrsta
  integer, dimension(:,:), allocatable :: iista
  


 !!! SCATTERER LOCATIONS
  

  ! Scatterer files ! NF it is better to include a preprocessing subroutine to create those files 
  character(120) :: scaradiusfile ! with j_r_index(=idummy), perr_(jr_index(1:r_prn))
  character(120) :: scaanglefile ! with j_angle_index(=idummy),pertheta(1:angle_prn),perphi(1:angle_prn)
  character(120) :: scalocfile ! with idummy, j_r_index(1:scatter_n),j_angle_index(1:scatter_n),perturbationVpv,ph,sv,sh,Qp,Qs

  ! Scatterer depths and locations

  integer :: r_prn ! number of vertical stacking points for receivers
  integer :: angle_prn ! number of lateral stacking points for receivers
  integer :: scatter_n ! number of receivers


  real(kind(0d0)), dimension(:), allocatable :: perr_ ! with r_n elements
  real(kind(0d0)), dimension(:), allocatable :: pertheta, perphi! with theta_n elements
  
  integer, dimension(:), allocatable :: j_r_index ! with r_n elements
  integer, dimension(:), allocatable :: j_angle_index ! with angle_n elements

  ! Scatterer perturabations
  
  real(kind(0d0)), dimension(:), allocatable :: perturbationVpv,perturbationVph,perturbationVsv,perturbationVsh
  real(kind(0d0)), dimension(:), allocatable :: perturbationQp, perturbationQs ! with scatter_n elements


  ! Interpolation matrices for vertical stacking points 
 
  real(kind(0d0)), dimension(:,:), allocatable :: ppsta
  integer, dimension(:,:), allocatable :: jjsta


  ! Number of synthetic components
  
  integer, parameter :: num_synn = 18 ! for 6 moment tensors or 9 for 3 single forces

  ! Shallow events (to avoid evanescent waves: see Kawai et al. 2006)
  
  real(kind(0d0)), parameter :: shallowdepth = 100.d0 


end module constants

module variables
  
  implicit none

  ! variables for the trial function
  
  integer:: nnlayer,nlay
  integer, dimension(:), allocatable :: nlayer
  integer:: nslay,nllay
  real(kind(0d0)), dimension(:), allocatable:: ra


  ! variables for 1D structure
  
  integer :: nzone ! number of zones in PREM way of description
  integer :: nsl, nll ! numbers of solid and liquid zones
  integer, dimension(:), allocatable :: iphase ! whether solid or liquid with nzone elements
  integer :: ndc,vnp
  
  real(kind(0d0)):: rmin,rmax ! model radius range
  real(kind(0d0)),dimension(:,:), allocatable:: vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa
  real(kind(0d0)),dimension(:),allocatable::vra,rho,kappa
  real(kind(0d0)),dimension(:),allocatable::ecKx !3*Kx=3A-4N
  real(kind(0d0)),dimension(:),allocatable::ecKy !3*Ky=3F+2N
  real(kind(0d0)),dimension(:),allocatable::ecKz !3*Kz=2F+C
  real(kind(0d0)),dimension(:),allocatable::mu,ecL,ecN,rhoinv,kappainv
  complex(kind(0d0)),dimension(:),allocatable:: coef1,coef2,coef
  
  ! frequency on which we are working on
  
  real(kind(0d0)) :: omega

  ! artificial damping attenuation 
  real(kind(0d0)) :: omegai

  ! variables for the matrices

  complex(kind(0d0)),dimension(:,:),allocatable:: a0,a1,a2,a,c, ctmp
  real(kind(0d0)),dimension(:),allocatable :: t
  real(kind(0d0)),dimension(:),allocatable :: h1x,h1y,h1z,h2L,h2N,h3ax,h3ay,h3az,h4aL,h4aN,h5ax,h5ay,h5az,h6aL,h6aN,h3x,h3y,h3z,h4L,h4N,h5x,h5y,h5z,h6L,h6N,h7x,h7y,h7z,h8L,h8N 
  real(kind(0d0),dimension(:,:),allocatable :: h3mx,h3my,h3mz,h5mx,h5my,h5mz,h4m1L,h4m1N,h4m2L,h4m2N,h6m1L,h6m1N,h6m2L,h6m2N
  real(kind(0d0)),dimension(:),allocatable:: p1,p2,p3
  complex(kind(0d0)),dimension(:),allocatable:: g0
  complex(kind(0d0)),dimension(:),allocatable:: d0
  complex(kind(0d0)):: g0tmp(2),g0dertmp(2) ! forward

  ! vectors for stack points
  
  integer,dimension(:),allocatable:: isp,issp,ilsp,jssp,jsp, ksp,lsp

  ! vectors for the output stack points
  
  integer,dimension(:),allocatable:: istazone
  integer,dimension(:),allocatable:: ksta  ! output stack point for g
  integer,dimension(:),allocatable:: jsta   ! output stack point for d

  ! variables for gridding
  
  integer,dimension(:),allocatable:: jjdr,kkdr
  real(kind(0d0)),dimension(:),allocatable:: vmin,gridpar,dzpar
  
  ! variables for l cut off
 
  integer:: kc,lsuf,sufzone,ismall,llog
  real(kind(0d0)):: maxamp
  
  ! variables for numerical integration
  complex(kind(0d0)):: anum(4,4,10),bnum(4,4,10)
    

  ! other variables
  
  integer:: ll(12),lli(12),llj(12)
  real(kind(0d0)):: eps,l2,lsq
  real(kind(0d0)),dimension(:),allocatable:: work
  complex(kind(0d0)),dimension(:),allocatable:: z,w,cwork
  
end module variables
  
module receiversphericals

  implicit none
  
  complex(kind(0d0)),dimension(:,:,:),allocatable:: dvecr,dvecdtr,dvecdpr
  real(kind(0d0)),dimension(:,:,:),allocatable:: plmr
  
end module receiversphericals

module scatterersphericals
  
  implicit none
  
  complex(kind(0d0)),dimension(:,:,:),allocatable :: dvecp,dvecdtp,dvecdpp
  real(kind(0d0)),dimension(:,:,:),allocatable:: plmp


end module scatterersphericals

module outputs
  
  implicit none
  
  complex(kind(0d0)),dimension(:,:,:),allocatable::synn
  complex(kind(0e0)),dimension(:,:,:),allocatable::synnsngl

end module outputs

  
