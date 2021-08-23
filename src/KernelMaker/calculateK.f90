subroutine calculateRSGT
  use parameters
  use tmpSGTs
  use angles
  use kernels
  
  ! subroutine for calculate RSGT video for ip and ith

  implicit none

  ! number of components ! attention! these definitions are also necessary in the main code !!
  integer, parameter :: num_tsgtPSV = 20
  integer, parameter :: num_rsgtPSV = 10
  integer, parameter :: num_synnPSV = 10
  integer, parameter :: num_tsgtSH = 10
  integer, parameter :: num_rsgtSH = 5
  integer, parameter :: num_synnSH = 5
  integer, parameter :: num_h3 = 6
  integer, parameter :: num_h4 = 6
  integer :: i_sgt,ift,jt
  
  rsgtF=dcmplx(0.d0)
  tsgtF=dcmplx(0.d0)
  tmpker=0.d0

  !  interpolate SGTs
  call clsgt(deltar(ip,ith),num_rsgtPSV,rsgtF(1:num_rsgtPSV,fmin:fmax),rsgtomega(1:num_rsgtPSV,fmin:fmax,1:theta_n))
  !call clsgt(deltas(ip,ith),num_tsgtPSV,tsgtF(1:num_tsgtPSV,fmin:fmax),tsgtomega(1:num_tsgtPSV,fmin:fmax,1:theta_n))
  ! Calculate the RSGT and TSGT: The 4th-order SSGT is reduced to 2nd-order
  ! after double-dotted with the moment tensor
  call rsgt2h3freq
  !call tsgt2h4freq

 


  do i_sgt=1,num_h3
     du=0.d0
     call vectorFFT_double(fmin,fmax,np1,h3(i_sgt,fmin:fmax),du(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)


     if(ibwfilt) then
        do ift = 0,nfilter
           call bwfilt(du(iWindowStart:iWindowEnd),duf(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
        enddo
     endif

     
     
     do ift=0,nfilter
        
        do jt=iWindowStart,iWindowEnd,jtstep_timeincrementV
           !tmpvideoker(i_sgt,ift,1+(jt-nt1(ift))/jtstep_timeincrementV)=duf(ift,jt) 
           tmpvideoker(i_sgt,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)=duf(ift,jt) 
        enddo
        
     enddo
  enddo

  return

end subroutine calculateRSGT


subroutine calculateTSGT
  use parameters
  use tmpSGTs
  use angles
  use kernels
  
  ! subroutine for calculate RSGT video for ip and ith

  implicit none

  ! number of components ! attention! these definitions are also necessary in the main code !!
  integer, parameter :: num_tsgtPSV = 20
  integer, parameter :: num_rsgtPSV = 10
  integer, parameter :: num_synnPSV = 10
  integer, parameter :: num_tsgtSH = 10
  integer, parameter :: num_rsgtSH = 5
  integer, parameter :: num_synnSH = 5
  integer, parameter :: num_h3 = 6
  integer, parameter :: num_h4 = 6
  integer :: i_sgt,ift,jt

  rsgtF=dcmplx(0.d0)
  tsgtF=dcmplx(0.d0)
  tmpker=0.d0

  !  interpolate SGTs
  !call clsgt(deltar(ip,ith),num_rsgtPSV,rsgtF(1:num_rsgtPSV,fmin:fmax),rsgtomega(1:num_rsgtPSV,fmin:fmax,1:theta_n))
  call clsgt(deltas(ip,ith),num_tsgtPSV,tsgtF(1:num_tsgtPSV,fmin:fmax),tsgtomega(1:num_tsgtPSV,fmin:fmax,1:theta_n))
  ! Calculate the RSGT and TSGT: The 4th-order SSGT is reduced to 2nd-order
  ! after double-dotted with the moment tensor
  !call rsgt2h3freq
  call tsgt2h4freq
  

  do i_sgt=1,num_h4
     du=0.d0
     call vectorFFT_double(fmin,fmax,np1,h3(i_sgt,fmin:fmax),du(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)


     if(ibwfilt) then
        do ift = 0,nfilter
           call bwfilt(du(iWindowStart:iWindowEnd),duf(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
        enddo
     endif

     
     
     do ift=0,nfilter
        do jt=iWindowStart,iWindowEnd,jtstep_timeincrementV
        !do jt=nt1(ift),nt2(ift),jtstep_timeincrementV
           !tmpvideoker(i_sgt,ift,1+(jt-nt1(ift))/jtstep_timeincrementV)=duf(ift,jt) 
           tmpvideoker(i_sgt,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)=duf(ift,jt) 
        enddo
        
     enddo
  enddo



  return

end subroutine calculateTSGT




subroutine calculateKernel
  use parameters
  use tmpSGTs
  use angles
  use kernels
  
  ! subroutine for calculate tmpker for ip and ith

  implicit none

  ! number of components ! attention! these definitions are also necessary in the main code !!
  integer, parameter :: num_tsgtPSV = 20
  integer, parameter :: num_rsgtPSV = 10
  integer, parameter :: num_synnPSV = 10
  integer, parameter :: num_tsgtSH = 10
  integer, parameter :: num_rsgtSH = 5
  integer, parameter :: num_synnSH = 5
  integer, parameter :: num_h3 = 6
  integer, parameter :: num_h4 = 6


  rsgtF=dcmplx(0.d0)
  tsgtF=dcmplx(0.d0)
  tmpker=0.d0

  !  interpolate SGTs
  call clsgt(deltar(ip,ith),num_rsgtPSV,rsgtF(1:num_rsgtPSV,fmin:fmax),rsgtomega(1:num_rsgtPSV,fmin:fmax,1:theta_n))
  call clsgt(deltas(ip,ith),num_tsgtPSV,tsgtF(1:num_tsgtPSV,fmin:fmax),tsgtomega(1:num_tsgtPSV,fmin:fmax,1:theta_n))
  ! Calculate the RSGT and TSGT: The 4th-order SSGT is reduced to 2nd-order
  ! after double-dotted with the moment tensor
  call rsgt2h3freq
  call tsgt2h4freq
  
  !   Traveltime kernel for isotropic P-wave speed perturbation. If nfilter>=3, 
  !   also calculate the group-delay kernels fro phase-delay kernels.
 
  !print *, "database"
  !print *, rsgtomega(1,8,theta_n/2)
  !print *, tsgtomega(1,8,theta_n/2)
  !print *, "interpolated"

  !print *, rsgtF(1,5:10)
  !print *, tsgtF(1,5:10)
 
  if((trim(paramWRT).eq.'alpha').or.(trim(paramWRT).eq.'all').or.&
      (trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'allV')) then
     if((mtype.eq.11).or.(mtype.eq.21).or.(mtype.eq.10).or.(mtype.eq.20)) then
        tmpker(1:4,:)=0.
        call isovpfreq
     endif
  endif
     
  !   Traveltime kernel for isotropic S-wave speed perturbation. If nfilter>=3, 
  !   also calculate the group-delay kernels fro phase-delay kernels.
  
  if((trim(paramWRT).eq.'beta').or.(trim(paramWRT).eq.'all').or.&
     (trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
     if((mtype.eq.12).or.(mtype.eq.22).or.(mtype.eq.32).or.(mtype.eq.10).or.(mtype.eq.20).or.(mtype.eq.30)) then
        tmpker(5:8,:)=0.
        !if((sym.ge.0.d0).and.(sym.le.360.d0)) tmpker(15,:)=0.
        call isovsfreq
     endif
  endif


  if(trim(paramWRT).eq.'serious') then
     call isovsfreq
     call isovpfreq
  endif
     
  ! for anisotropy I will complete it later
  
  !   Kernels for epsilon when s=r (rad. aniso.) and 0<=sym<=360 (azim. aniso.).
  
  ! if((mtype.eq.10).or.(mtype.eq.20).or.(mtype.eq.30)) then
  !    tmpker(9:10,:)=0.
  !    if((sym.ge.0.d0).and.(sym.le.360.d0)) tmpker(16,:)=0.
  !    call anepsiln
  ! endif
	            	      
  !   Kernels for delta when s=r (rad. aniso.) and 0<=sym<=360 (azim. aniso.).
  
  ! if((mtype.eq.10).or.(mtype.eq.20).or.(mtype.eq.30)) then
  !    tmpker(11:12,:)=0.
  !    if((sym.ge.0.d0).and.(sym.le.360.d0)) tmpker(17,:)=0.
  !    call andelta
  ! endif
  
  !   Kernels for gamma when s=r (in rad. aniso.) and 0<=sym<=360 (azim. aniso.).
  
  ! if((mtype.eq.10).or.(mtype.eq.20).or.(mtype.eq.30)) then
  !    tmpker(13:14,:)=0.
  !    if((sym.ge.0.d0).and.(sym.le.360.d0)) then
  !       tmpker(18,:)=0.
  !       tmpker(19,:)=0.
  !       tmpker(20,:)=0.
  !       call angammacs
  !    endif
  !    if(sym.lt.0.d0) call angamma
  ! endif
  
  
  return
end subroutine calculateKernel

subroutine isovpfreq
  use parameters
  use tmpSGTs
  use kernels
  implicit none
  
  integer :: jt,ift,jtstep

  !   This subroutine calculates 3 types of kernels:
  !
  !      1: Phase-delay time sensitivity to P-wave speed
  !      3: Amplitude perturbation sensitivity to P-wave speed
  !      4: Amplitude perturbation sensitivity to Q_kappa
  !
  !   Note: phase-delay time is the same as traveltime perturbation obtained 
  !   by cross-correlation.
  !
  !   Type 2 is for group-delay time sensitivity to P-wave speed. It can be 
  !   calculated later from Type 1 by numerical differentiation wrt frequency.

    
  du=0.d0
  
  u_freq=dcmplx(0.d0)
  do jt = fmin,fmax
     u_freq(jt)=-(h3(1,jt)+h3(2,jt)+h3(3,jt))* &
          (h4(1,jt)+h4(2,jt)+h4(3,jt))
     u_freq(jt) = u_freq(jt)*dcmplx(1.d3)
  enddo
  call vectorFFT_double(fmin,fmax,np1,u_freq(fmin:fmax),du(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)
  if(ibwfilt) then
     do ift = 0,nfilter
        call bwfilt(du(iWindowStart:iWindowEnd),duf(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
     enddo
  endif

  if(trim(paramWRT).eq.'serious') then
     open (111,file=seriousfrechetfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*(iWindowEnd-iWindowStart+1)
     write(111) coeffV(1,ir)*duf(0,iWindowStart:iWindowEnd)*1.d3
  endif



  do ift=0,nfilter
     if((trim(paramWRT).eq.'alphaV').or.(trim(paramWRT).eq.'allV')) then
        !do jt=nt1(ift),nt2(ift),jtstep_timeincrementV
        do jt=iWindowStart,iWindowEnd,jtstep_timeincrementV
           !tmpvideoker(1,ift,1+(jt-nt1(ift))/jtstep_timeincrementV)=coeffV(1,ir)*duf(ift,jt) 
           tmpvideoker(1,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)=coeffV(1,ir)*duf(ift,jt)*1.d30 
           !print *, duf(ift,jt),coeffV(1,ir),tmpvideoker(1,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)
        enddo
     endif
  enddo    

  do ift=0,nfilter
     tmpker(1,ift)=tmpker(1,ift) &
          +coeff(ift,1,ir,nt1(ift))*duf(ift,nt1(ift))/2.d0
     tmpker(3,ift)=tmpker(3,ift) &
          +coeff(ift,3,ir,nt1(ift))*duf(ift,nt1(ift))/2.d0
     tmpker(4,ift)=tmpker(4,ift) &
          +coeff(ift,4,ir,nt1(ift))*duf(ift,nt1(ift))/2.d0
     do jt=nt1(ift)+1,nt2(ift)-1
        tmpker(1,ift)=tmpker(1,ift)+coeff(ift,1,ir,jt)*duf(ift,jt)
        tmpker(3,ift)=tmpker(3,ift)+coeff(ift,3,ir,jt)*duf(ift,jt)
        tmpker(4,ift)=tmpker(4,ift)+coeff(ift,4,ir,jt)*duf(ift,jt)
     enddo
     tmpker(1,ift)=tmpker(1,ift) &
          +coeff(ift,1,ir,nt2(ift))*duf(ift,nt2(ift))/2.d0
     tmpker(3,ift)=tmpker(3,ift) &
          +coeff(ift,3,ir,nt2(ift))*duf(ift,nt2(ift))/2.d0
     tmpker(4,ift)=tmpker(4,ift) &
          +coeff(ift,4,ir,nt2(ift))*duf(ift,nt2(ift))/2.d0
  enddo
  return
end subroutine isovpfreq



subroutine isovsfreq
  use parameters
  use tmpSGTs
  use kernels
  implicit none
  real(kind(0.e0)) :: par(0:nfilter,1:nt2(0)-nt1(0)+1)
  real(kind(0.e0)) :: parq(0:nfilter,1:nt2(0)-nt1(0)+1) 
  integer :: jt,ift,jtstep
  character(250) :: kerfile 
  !   This subroutine calculates 4 types of kernels:
  !
  !      5: Phase-delay time sensitivity to S-wave speed
  !      7: Amplitude perturbation sensitivity to S-wave speed
  !      8: Amplitude perturbation sensitivity to Q_mu
  !     15: SKS-splitting intensity sensitivity to isotropic S-wave speed
  !         in an anisotropic model with a horizontal axis of symmetry 
  !         oriented at azimuth sym.
  !
  !   Note: phase-delay time is the same as traveltime perturbation obtained 
  !   by cross-correlation.
  !
  !
  !    8 and 15 is not working for this version
  !
  !
  !   rewrite the duq with Fuji et al. 2010 !!!!
  !
  !
  
  !   For Types 5, 7 and 8.
  !   Convolve the two SGTs to obtained waveform partial derivative.
  
  du=0.d0
  duq=0.d0

  par=0.e0
  parq=0.e0
  u_freq=dcmplx(0.d0)
  do jt = fmin,fmax
     u_freq(jt) = (-2.d0*(h3(4,jt)*h4(4,jt) &
          +h3(5,jt)*h4(5,jt) &
          +h3(6,jt)*h4(6,jt))+h3(1,jt)*(h4(2,jt)+h4(3,jt)) &
          +h3(2,jt)*(h4(1,jt)+h4(3,jt)) &
          +h3(3,jt)*(h4(1,jt)+h4(2,jt)))
     u_freq(jt) = u_freq(jt)*dcmplx(1.d3)
  enddo
  call vectorFFT_double(fmin,fmax,np1,u_freq(fmin:fmax),du(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)
  
  u_freq=dcmplx(0.d0)


  ! old version 
  !do jt = fmin,fmax
  !   u_freq(jt) = (-2.d0*(h3(1,jt)*h4(1,jt) &
  !        +h3(2,jt)*h4(2,jt) &
  !        +h3(3,jt)*h4(3,jt))/3.d0-2.d0*(h3(4,jt)*h4(4,jt) &
  !        +h3(5,jt)*h4(5,jt)+h3(6,jt)*h4(6,jt)) &
  !        +h3(1,jt)*(h4(2,jt)+h4(3,jt))/3.d0 &
  !        +h3(2,jt)*(h4(1,jt)+h4(3,jt))/3.d0 &
  !        +h3(3,jt)*(h4(1,jt)+h4(2,jt))/3.d0)
  !   u_freq(jt) = u_freq(jt)*dcmplx(1.d3)
  !enddo
  

  ! calculate K^u_q 



  do jt = fmin, fmax
     ! first calculate K^u_mu with dmu0=1Pa
     u_freq(jt) = -4.d0*(h3(4,jt)*h4(4,jt)+h3(5,jt)*h4(5,jt)+h3(6,jt)*h4(6,jt)) &
          + 2.d0/3.d0*(h3(1,jt)*(h4(2,jt)+h4(3,jt)) + h3(2,jt)*(h4(1,jt)+h4(3,jt)) + h3(3,jt)*(h4(1,jt)+h4(2,jt))) &
          - 4.d0/3.d0*(h3(1,jt)*h4(1,jt)+h3(2,jt)*h4(2,jt)+h3(3,jt)*h4(3,jt))
     ! then apply the Jacobian of Fuji et al. 2010
     u_freq(jt) = u_freq(jt)*jacobianFuji(ir,jt)
     u_freq(jt) = u_freq(jt)*dcmplx(1.d3) 
  enddo

  call vectorFFT_double(fmin,fmax,np1,u_freq(fmin:fmax),duq(iWindowStart:iWindowEnd),omegai,tlen,iWindowStart,iWindowEnd)
  

  if(ibwfilt) then
     do ift = 0,nfilter
        call bwfilt(du(iWindowStart:iWindowEnd),duf(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
        call bwfilt(duq(iWindowStart:iWindowEnd),duqf(ift,iWindowStart:iWindowEnd),1.d0/samplingHz,(iWindowEnd-iWindowStart+1),0,npButterworth,fclp(ift),fchp(ift))
     enddo
  endif


!!$
!!$  jtstep=1
!!$  if(timeincrement.ne.0.d0) jtstep=int(timeincrement*samlingHz)
!!$
!!$  do ift=0,nfilter
!!$     if((trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
!!$        do jt=nt1(ift),nt2(ift),jtstep
!!$           par(ift,(jt-nt1(ift))/jtstep+1)=real(du(jt))
!!$          parq(ift,(jt-nt1(ift))/jtstep+1)=real(duq(jt)) 
!!$        enddo
!!$             write(tmpchar,'(I7)') int(rx*1.d3)
!!$     do j=1,7
!!$        if(tmpchar(j:j).eq.' ') tmpchar(j:j) = '0'
!!$     enddo
!!$     kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
!!$            //"."//ir//"."//ith//"."//ip//"."//"video"
!!$
!!$     open(1,file=kerfile,status='unknown',form='unformatted', &
!!$          access = 'direct', recl=kind(0e0)*(nfilter+1)*(int((nt2(ift)-nt1(ift))/jtstep)+1))
!!$     write(1,rec=1) par(0:nfilter,1:int((nt2(ift)-nt1(ift))/jtstep)+1)
!!$     close(1)
!!$
!!$     kerfile=trim(parentDir)//"/tmpvideo/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(tmpchar)&
!!$            //"."//ir//"."//ith//"."//ip//"."//"videoq"
!!$
!!$     open(1,file=kerfile,status='unknown',form='unformatted', &
!!$          access = 'direct', recl=kind(0e0)*(nfilter+1)*(int((nt2(ift)-nt1(ift))/jtstep)+1))
!!$     write(1,rec=1) parq(0:nfilter,1:int((nt2(ift)-nt1(ift))/jtstep)+1)
!!$     close(1)
!!$ 
!!$
!!$    endif
!!$  enddo
  
  if(trim(paramWRT).eq.'serious') then
     open (111,file=seriousfrechetfile,status='unknown',form='unformatted',access='direct',recl=kind(0e0)*(iWindowEnd-iWindowStart+1)
     write(111) coeffV(2,ir)*duf(0,iWindowStart:iWindowEnd)*1.d3
  endif

  
  do ift=0,nfilter
     if((trim(paramWRT).eq.'betaV').or.(trim(paramWRT).eq.'allV')) then
        !do jt=nt1(ift),nt2(ift),jtstep_timeincrementV
        do jt=iWindowStart,iWindowEnd,jtstep_timeincrementV
           !tmpvideoker(2,ift,1+(jt-nt1(ift))/jtstep_timeincrementV)=coeffV(2,ir)*duf(ift,jt) ! 
           tmpvideoker(2,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)=coeffV(2,ir)*duf(ift,jt)*1.d30
           tmpvideoker(3,ift,1+(jt-iWindowStart)/jtstep_timeincrementV)=coeffV(2,ir)*duq(jt)*1.d30
         
        enddo
     endif
  enddo    
  
  do ift=0,nfilter
     tmpker(5,ift)=tmpker(5,ift) &
          +0.5d0*coeff(ift,5,ir,nt1(ift))*duf(ift,nt1(ift))     
     tmpker(7,ift)=tmpker(7,ift) &
          +0.5d0*coeff(ift,7,ir,nt1(ift))*duf(ift,nt1(ift))
     tmpker(8,ift)=tmpker(8,ift) &
          +0.5d0*coeff(ift,8,ir,nt1(ift))*duqf(ift,nt1(ift))
     do jt=nt1(ift)+1,nt2(ift)-1
        tmpker(5,ift)=tmpker(5,ift)+coeff(ift,5,ir,jt)*duf(ift,jt)
        tmpker(7,ift)=tmpker(7,ift)+coeff(ift,7,ir,jt)*duf(ift,jt)
        tmpker(8,ift)=tmpker(8,ift)+coeff(ift,8,ir,jt)*duqf(ift,jt)
     enddo
     tmpker(5,ift)=tmpker(5,ift) &
          +0.5d0*coeff(ift,5,ir,nt2(ift))*duf(ift,nt2(ift))
     tmpker(7,ift)=tmpker(7,ift) &
          +0.5d0*coeff(ift,7,ir,nt2(ift))*duf(ift,nt2(ift))
     tmpker(8,ift)=tmpker(8,ift) &
          +0.5d0*coeff(ift,8,ir,nt2(ift))*duqf(ift,nt2(ift))
  enddo

  if((sym.lt.0.d0).or.(sym.gt.360.d0)) return

  ! I have to work for the type 15

  
  return
end subroutine isovsfreq
  


