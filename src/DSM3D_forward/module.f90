
module parameters

  implicit none
  



module parameters_for_KernelMaker

  implicit none

  !   Constants.
  !
  !
  !
  !     ddeltas, ddeltar - distance sampling intervals in TSGT and RSGT.
  !     nfilter          - number of narrow-band filters (0 for broad-band only)
  !     nktype           - number of different types of kernels, including
  !                        phase delay for alpha   (1) and beta (5)
  !                        group delay for alpha   (2) and beta (6)
  !                        amplitude for alpha     (3) and beta (7)
  !                        amplitude for Q_kappa   (4) and Q_mu (8)
  !
  !                       9-20: I will include in the near future!
  !
  !                             So nktype = 8
  !
  !                        phase delay for epsilon (9)
  !                        amplitude for epsilon   (10)
  !                        phase delay for delta   (11)
  !                        amplitude for delta     (12)
  !                        phase delay for gamma   (13)
  !                        amplitude for gamma     (14)
  !                        split. int. for beta    (15)
  !                        split. int. for epsilon (16)
  !                        split. int. for delta   (17)
  !                        split. int. for gamma   (18)
  !                        split. int. for gamma_c (19)
  !                        split. int. for gamma_s (20)
  !                        Def. for epsilon, delta and gamma see Chevrot (2006).
  !                        Both radial and azimuthal aniso. are implemented,
  !                        i.e. the axis of symmetry s=r, or s is horizontal with
  !                        azimuth 0 <= sym <= 360.
  !     nkvtype          - number of different types of video kernels (for the moment 2)
  !                        synthetic 0, alpha 1, beta 2 (for attenuation, anisotropy, they will soon implemented)
  !
  !
  !   Feb. & Mar. 2010 LZ:    complete debug of isotropic kernels (nktype:1-8)
  !   Jul. & Aug. 2010 LZ:    implementation of anisotropic kernels (nktype:9-20)
  !   Agu. 2010 LZ:           implementation of zero-phase filter
  !   Sep. 2010 LZ:           complete debug of anisotropic kernels
  !   Apr. 2014 NF:  for video kernels
  !


  integer, parameter :: nfilter=0 ! it is ready for the extention nfilter > 0
  integer, parameter :: nktype = 8 ! isotropic kernels
  integer, parameter :: nkvtype = 2 ! isotropic video kernels
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0
  character(120) :: SGTinfo,parentDir,eventName,stationName,phase,compo,paramWRT
  character(120) :: synnfile
  character(120) :: INFO_TSGT,INFO_RSGT,rsampletxt,modelcard
  character(120) :: Poutputdir,psvmodel,modelname,DSMconfFile
  integer :: mtype
  real(kind(0d0)) :: evla, evlo, evdepth,stla,stlo
  real(kind(0d0)) :: mt(1:6),twin(1:4)
  integer :: itwin(1:4)
  real(kind(0d0)) :: ipdistance,c_red_reci
  integer :: ibwfilt,npButterworth,ifastFFT,fmin,fmax ! filter, butterworth
  integer :: itranslat ! geodetic <-> geocentric
  character(120) :: freqid(0:nfilter)
  real(kind(0d0)) :: fclp(0:nfilter), fchp(0:nfilter)
  real(kind(0d0)) :: dph,ph1,dth,thw ! horizontal region of partial calculation
  real(kind(0d0)) :: start, end ! time window for SGT vectors
  integer :: iwindowStart, iwindowEnd ! integer version of start and end
  real(kind(0d0)) :: samplingHz, dtn

  ! for ignoring scheme
  real(kind(0d0)) :: calculrapide
  integer :: nntype
  integer, allocatable :: idecidetype(:)
  real(kind(0d0)), allocatable :: minici(:,:)
  real(kind(0d0)) :: MinLengthFor0
  integer :: nMinLengthFor0


  ! Anisotropy (here I use isotropic case sym=0.d0)
  real(kind(0d0)), parameter :: sym = 0.d0
  real(kind(0d0)) :: csym,ssym


  ! Kernel
  integer :: nr,ntot
  real(kind(0d0)) :: rmin,rmax,rdelta
  real(kind(0d0)), allocatable :: r(:)
  real(kind(0d0)), allocatable :: rhom(:),vpm(:),vsm(:),qmm(:),qkp(:)
  complex(kind(0d0)), allocatable :: jacobianFuji(:,:)
  !real(kind(0d0)), allocatable :: vpv(:),vph(:),vsv(:),vsh(:),eta(:) ! for the future extention
  real(kind(0d0)) :: rs,rr
  real(kind(0d0)), parameter :: rdep = 0.d0 ! surface
  real(Kind(0d0)), parameter :: rEarth = 6371.d0 ! kilometres


  ! DSM
  integer :: iPSVSH
  real(kind(0d0)) :: re, ratc, ratl, omegai
  real(kind(0d0)), allocatable :: omega(:)
  real(kind(0d0)) :: tlen
  real(kind(0d0)) :: rmin_, rmax_, rdelta_
  real(kind(0d0)) :: r0min, r0max, r0delta  !!! JUST FOR ONE DEPTH FOR THIS MOMENT !!
  real(kind(0d0)) :: thetamin, thetamax, thetadelta
  real(kind(0d0)), allocatable :: r_(:),r0D(:),thetaD(:)
  integer :: imin, imax
  integer :: rsgtswitch, tsgtswitch, synnswitch

  real(kind(0d0)),allocatable:: vrminD(:),vrmaxD(:),rrhoD(:,:),vpvD(:,:),vphD(:,:)
  real(kind(0d0)),allocatable:: vsvD(:,:),vshD(:,:),etaD(:,:),qmuD(:),qkappaD(:)


  integer :: nzone
  integer :: r_n,r0_n,ciista, ir_,ir0,imt,icomp,idepth,itheta, theta_n
  integer,allocatable :: ksegr(:),lrads(:)
  integer,parameter :: maxradsamples = 5
  ! FOR FOURIER TRANSFORM & BUTTERWORTH FILTERING


  integer :: lsmooth,np,np0,np1,np1_fk
  integer :: n1,m1




end module parameters_for_KernelMaker


module angles
  implicit none
  real(kind(0d0)), allocatable :: phi00(:),phi0(:),theta0(:)
  integer :: nphi0, ntheta,nphi
  real(kind(0d0)), allocatable :: phitheta(:,:),thetaphi(:,:)
  real(kind(0d0)), allocatable :: phi(:,:), theta(:,:)
  real(kind(0d0)), allocatable :: crq(:,:),srq(:,:),crq2(:,:),srq2(:,:)
  real(kind(0d0)), allocatable :: csq(:,:),ssq(:,:),csq2(:,:),ssq2(:,:)
  real(kind(0d0)), allocatable :: cqs(:,:),sqs(:,:),cqs2(:,:),sqs2(:,:)
  real(kind(0d0)), allocatable :: deltar(:,:),deltas(:,:)
  real(kind(0d0)) :: slat,slon,sdep,rlat,rlon

end module angles



module tmpSGTs
  implicit none
  ! in frequency domain
  complex(kind(0d0)), allocatable :: synnF(:,:),rsgtF(:,:),tsgtF(:,:)
  complex(kind(0d0)), allocatable :: h3(:,:),h4(:,:)
  complex(kind(0d0)), allocatable :: synnomega(:,:,:),rsgtomega(:,:,:),tsgtomega(:,:,:)
  complex(kind(0d0)), allocatable :: u_freq(:)
  ! in time domain
  real(kind(0d0)), allocatable :: t(:),u(:),u0(:,:),v(:),v0(:,:),hu(:),hu0(:,:)
  real(kind(0d0)), allocatable :: fwin(:,:)
  integer, allocatable :: nt1(:),nt2(:)
  real(kind(0d0)), allocatable :: denomv(:),denomu(:),coeff(:,:,:,:),coeffV(:,:)

end module tmpSGTs


module kernels
  implicit none

  ! parameters
  real(kind(0d0)) :: rx,xlat,xlon
  integer :: ip,ith,ir

  ! for the stockage
  real(kind(0e0)), allocatable :: ker(:,:,:,:)

  ! for calculation
  real(kind(0d0)), allocatable :: tmpker(:,:)

  ! tmparray for video (already in single precision)
  real(kind(0e0)), allocatable :: videoker(:,:,:,:) ! unlike ker(:....), we stock along ip but not along ith (theta)
  real(kind(0e0)), allocatable :: tmpvideoker(:,:,:)
  integer :: number_of_snapshots
  real(kind(0d0)) :: timeincrementV ! time interval for video
  integer :: jtstep_timeincrementV ! in integer

  ! tmparrays

  complex(kind(0d0)), dimension (:), allocatable :: tmph01,tmph02,tmph03,tmph04
  complex(kind(0d0)), dimension (:), allocatable :: tmph05,tmph06,tmph07,tmph08
  complex(kind(0d0)), dimension (:), allocatable :: tmph09,tmph10,tmph11,tmph12
  complex(kind(0d0)), dimension (:), allocatable :: tmph13,tmph14,tmph15,tmph16
  complex(kind(0d0)), dimension (:), allocatable :: tmph17,tmph18,tmph19,tmph20
  complex(kind(0d0)), dimension (:), allocatable :: tmph21,tmph22,tmph23,tmph24
  complex(kind(0d0)), dimension (:), allocatable :: tmph25,tmph26

  ! in fact those arrays are not necessary because it was from old version of ZC2011a
  ! they had 36 components (but 20 components are independent)
  ! however, my actual code seems to have some problems (in terms of signs)
  ! so I use this for the moment

  ! partial derivatives in time domain

  real(kind(0d0)), allocatable :: du(:),duf(:,:),duq(:),duqf(:,:)


end module kernels


module rotate
  real(kind(0d0)) :: cc(3,3),ct(3,3)
end module rotate
