subroutine rdsgtomega(rx,ry,num_sgt,num_psv,ipsvorsh)

  use parameters
  use tmpSGTs
  implicit none
  real(kind(0d0)) :: rx,ry,redtime ! in order to read the catalogue
  integer :: i,j,num_sgt,ipsvorsh,num_psv
  complex(kind(0e0)) :: sgtsngl(1:num_sgt,1:theta_n)
  complex(kind(0d0)) :: sgtdouble(1:num_sgt,1:theta_n)
  character(200) :: coutfile

  !print *, ipsvorsh
  ! ATTENTION C'EST PAS NORMAL
  ! rx = 5751.d0
  ! rx = 6351.d0


  ! about ipsvorsh
  ! PSV synn =  2; SH synn =  1
  ! PSV rsgt = 20; SH rsgt = 10
  ! PSV tsgt =200; SH tsgt =100

  if((ipsvorsh.eq.2).or.(ipsvorsh.eq.20).or.(ipsvorsh.eq.200)) then
     if(num_sgt.ne.num_sgt) then
        print *, "the number of SGTs is not propre."
        stop
     endif
  endif
  

  do i = fmin,fmax
     sgtsngl = cmplx(0.e0)
     sgtdouble = cmplx(0.d0)
     if(i.ne.0) then
        if(ipsvorsh.eq.2) then
           write(coutfile, '(I7,".",I7,".SYNN_PSV")') int(rx*1.d3),i

           ! just for depth600 because of BUGG in SGTpsv!!!!
           !write(coutfile, '(I7,".",I7,".SYNN_PSV")') 0,i
           !!!!!! ATTENTION !!!! PLEASE THIS IS NOT PROPRE!!!!!
           
        elseif(ipsvorsh.eq.1) then
           write(coutfile, '(I7,".",I7,".SYNN_SH")') int(rx*1.d3),i           
        elseif(ipsvorsh.eq.20) then
           write(coutfile, '(I7,".",I7,".RSGT_PSV")') int(rx*1.d3),i
        elseif(ipsvorsh.eq.10) then
           write(coutfile, '(I7,".",I7,".RSGT_SH")') int(rx*1.d3),i
        elseif(ipsvorsh.eq.200) then
           write(coutfile, '(I7,".",I7,".",I7,".TSGT_PSV")') int(rx*1.d3),int(ry*1.d3),i
        elseif(ipsvorsh.eq.100) then
           write(coutfile, '(I7,".",I7,".",I7,".TSGT_SH")') int(rx*1.d3),int(ry*1.d3),i
        endif
        
        if((ipsvorsh.eq.200).or.(ipsvorsh.eq.100)) then ! TSGT
           do j = 1,29
                 if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
           enddo
           coutfile = trim(modelname)//"."//coutfile
           coutfile = trim(PoutputDir)//"/TSGT/"//coutfile
        else
           do j = 1,21
              if (coutfile(j:j).eq.' ')coutfile(j:j) = '0'
           enddo
           coutfile = trim(modelname)//"."//coutfile
           coutfile = trim(PoutputDir)//"/RSGT/"//coutfile ! RSGT & SYNN
        endif
        open(1,file=coutfile,status='unknown',form='unformatted', &
             access = 'direct', recl=2*num_sgt*kind(0e0)*theta_n)       
        read(1,rec=1)sgtsngl(1:num_sgt,1:theta_n)    
        sgtdouble(1:num_sgt,1:theta_n) = sgtsngl(1:num_sgt,1:theta_n) 


        ! h9 in SH is 0 !
        if((ipsvorsh.eq.1).or.(ipsvorsh.eq.10)) then
            !sgtdouble(4,1:theta_n) = cmplx(0.d0)
           !sgtdouble(3:5,1:theta_n) = sgtdouble(3:5,1:theta_n)*5.d-1
           !sgtdouble(3:5,1:theta_n) = -sgtdouble(3:5,1:theta_n)
           
        endif
        
        if(ipsvorsh.eq.1) then
           sgtdouble(3:5,1:theta_n) = sgtdouble(3:5,1:theta_n)*5.d-1
           !sgtdouble(1:2,1:theta_n) = cmplx(0.d0)
        endif

        if(ipsvorsh.eq.2) then
           sgtdouble(2,1:theta_n) = sgtdouble(2,1:theta_n)*5.d-1
           sgtdouble(8:10,1:theta_n) = sgtdouble(8:10,1:theta_n)*5.d-1
        endif



        close(1)   
            
        do itheta = 1,theta_n
           redtime = -thetaD(itheta)*c_red_reci 
           if(ipsvorsh.ge.10) then
              sgtdouble(1:num_sgt,itheta) = sgtdouble(1:num_sgt,itheta) &
                   *exp(cmplx(0.d0,-2.d0*pi*dble(i)/tlen*redtime))*1.d-21
           else
              sgtdouble(1:num_sgt,itheta) = sgtdouble(1:num_sgt,itheta) &
                   *exp(cmplx(0.d0,-2.d0*pi*dble(i)/tlen*redtime))*1.d-21
           endif
                 
        enddo                
        !rsgtomega(1:1,i,1:theta_n) = rsgtomega(1,i,1:theta_n) + rsgtdouble(1,1:theta_n)
        if(ipsvorsh.eq.2) then
           synnomega(1:num_psv,i,1:theta_n)=synnomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        elseif(ipsvorsh.eq.20) then
           rsgtomega(1:num_psv,i,1:theta_n)=rsgtomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        elseif(ipsvorsh.eq.200) then
           tsgtomega(1:num_psv,i,1:theta_n)=tsgtomega(1:num_psv,i,1:theta_n)+sgtdouble(1:num_psv,1:theta_n)
        elseif(ipsvorsh.eq.1) then
           synnomega(6:10,i,1:theta_n)=synnomega(6:10,i,1:theta_n)+sgtdouble(1:5,1:theta_n)
        elseif(ipsvorsh.eq.10) then
           rsgtomega(6:10,i,1:theta_n)=rsgtomega(6:10,i,1:theta_n)+sgtdouble(1:5,1:theta_n)
        elseif(ipsvorsh.eq.100) then
           tsgtomega(7,i,1:theta_n) = tsgtomega(7,i,1:theta_n) + sgtdouble(1,1:theta_n)
           tsgtomega(12:20,i,1:theta_n)= tsgtomega(12:20,i,1:theta_n)+ sgtdouble(2:10,1:theta_n)
        endif
     endif
  enddo
      
  return
end subroutine rdsgtomega




subroutine clsgt(distance,num_sgt,sgtF,sgtomega)
  use parameters
  implicit none
  integer :: idelta,nhregion,idstart,idend,kk,j,num_sgt
  real(kind(0d0)) :: distance,tmpnakami,normalisefactorInt
  complex(kind(0d0)) :: sgtF(1:num_sgt,fmin:fmax),sgtomega(1:num_sgt,fmin:fmax,1:theta_n)


  idelta = 1+int((distance-thetaD(1))/thetadelta)
  nhregion = int(ipdistance/thetadelta/2.d0)
  call findIdeltaRegion(idelta,nhregion,theta_n,idstart,idend)
  
  normalisefactorInt = 0.d0
  
  do kk = idstart, idend
     tmpnakami = pi*((distance-thetaD(idstart))/thetadelta-dble(kk-idstart))
     if(tmpnakami.eq.0.d0) then
        sgtF(1:num_sgt,fmin:fmax) = sgtF(1:num_sgt,fmin:fmax) + sgtomega(1:num_sgt,fmin:fmax,kk)
        normalisefactorInt = normalisefactorInt+1.d0
     else
        sgtF(1:num_sgt,fmin:fmax) = sgtF(1:num_sgt,fmin:fmax)  &
             + dsin(tmpnakami)/tmpnakami*sgtomega(1:num_sgt,fmin:fmax,kk)
        normalisefactorInt = normalisefactorInt+dsin(tmpnakami)/tmpnakami
     endif
  enddo
  
  do j=fmin,fmax
     sgtF(1:num_sgt,j)=sgtF(1:num_sgt,j)/normalisefactorInt*exp(cmplx(0.d0,-2.d0*pi*dble(j)/tlen*distance*c_red_reci))
  enddo  
  
  return
end subroutine clsgt




subroutine findIdeltaRegion(ix,nhalf,nmax,istart,iend)
  implicit none
  integer ix,nhalf,nmax,istart,iend
  
  istart = ix - nhalf
  iend = ix + nhalf
  
  if(istart.lt.1) then 
     iend = iend + (1-istart)
     istart = 1
  endif
  if(iend.gt.nmax) then
     istart = istart - (iend-nmax)
     iend = nmax
  endif
  
  if(istart.lt.1) then
     istart = 1
  endif
  if(iend.gt.nmax) then
     iend = nmax
  endif
  
  return
end subroutine findIdeltaRegion


!

subroutine hilbert(n,dt,t,f,h)

  !   This subroutine calculates the Hilbert transform h(x) of a function 
  !   f(t) by a simple numerical integration. Since Hilbert transform is 
  !   expressed as an integral of f(t)/(t-x) over the entire t axis, for
  !   this simple program to work the function f(t) can be non-zero only 
  !   over a finite range of positive t. 
  !   Call this subroutine as the following: 
  !
  !       call hilbert(n,dt,t,f,h)
  
  implicit none
  real(kind(0d0)) :: dt,t(0:n),f(0:n),h(0:n)
  real(kind(0d0)) :: pi,const
  integer :: i,j,n


  pi=4.d0*datan(1.d0)
  
  h=0.d0
  
  do i=0,n
     do j=0,n
        const=1.d0
        if((j.eq.0).or.(j.eq.n)) const=0.5d0
        if(i.ne.j) h(i)=h(i)+const*f(j)*dt/(t(j)-t(i))
     enddo
     h(i)=h(i)/pi
  enddo
  
  return
end subroutine hilbert

!


subroutine synn2h3freq(ip,ith)
          
  !   This subroutine calculates the 18 independent elements of the 3rd-order
  !   SSGT from the 10 azimuth independent coefficients and the sines and 
  !   cosines of the azimuths.
	
  use parameters
  use angles
  use tmpSGTs
  implicit none
  integer :: ip,ith
  
  
  !   Vertical component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.

  if(trim(compo).eq."Z") then
     h3(1,fmin:fmax)=synnF(1,fmin:fmax)
     h3(2,fmin:fmax)=synnF(3,fmin:fmax)*crq2(ip,ith)-synnF(4,fmin:fmax)
     h3(3,fmin:fmax)=-synnF(3,fmin:fmax)*crq2(ip,ith)-synnF(4,fmin:fmax)
     h3(4,fmin:fmax)=-synnF(2,fmin:fmax)*crq(ip,ith)
     h3(5,fmin:fmax)=-synnF(2,fmin:fmax)*srq(ip,ith)
     h3(6,fmin:fmax)=synnF(3,fmin:fmax)*srq2(ip,ith)
  
  !   Radial component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.

  elseif(trim(compo).eq."R") then
     h3(1,fmin:fmax)=synnF(5,fmin:fmax)
     h3(2,fmin:fmax)=-(synnF(7,fmin:fmax))*crq2(ip,ith)-synnF(6,fmin:fmax)
     h3(3,fmin:fmax)=(synnF(7,fmin:fmax))*crq2(ip,ith)-synnF(6,fmin:fmax)
     h3(4,fmin:fmax)=(synnF(9,fmin:fmax))*crq(ip,ith)
     h3(5,fmin:fmax)=(synnF(9,fmin:fmax))*srq(ip,ith)
     h3(6,fmin:fmax)=-(synnF(7,fmin:fmax))*srq2(ip,ith)

  !   Transverse component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.

  elseif(trim(compo).eq."T") then
     h3(1,fmin:fmax)=0.d0
     h3(2,fmin:fmax)=(synnF(8,fmin:fmax))*srq2(ip,ith)
     h3(3,fmin:fmax)=-h3(2,fmin:fmax)
     h3(4,fmin:fmax)=-(synnF(10,fmin:fmax))*srq(ip,ith)
     h3(5,fmin:fmax)=(synnF(10,fmin:fmax))*crq(ip,ith)
     h3(6,fmin:fmax)=-(synnF(8,fmin:fmax))*crq2(ip,ith)
  endif
  return
end subroutine synn2h3freq

subroutine rsgt2h3freq
  use parameters
  use angles
  use tmpSGTs
  use kernels
  implicit none


  !   This subroutine calculates the 18 independent elements of the 3rd-order
  !   SSGT from the 10 azimuth independent coefficients and the sines and 
  !   cosines of the azimuths.


  
  if(trim(compo).eq."Z") then
     !   Vertical component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.
     h3(1,fmin:fmax)=rsgtF(1,fmin:fmax)
     h3(2,fmin:fmax)=rsgtF(3,fmin:fmax)*crq2(ip,ith)-rsgtF(4,fmin:fmax)
     h3(3,fmin:fmax)=-rsgtF(3,fmin:fmax)*crq2(ip,ith)-rsgtF(4,fmin:fmax)
     h3(4,fmin:fmax)=-rsgtF(2,fmin:fmax)*crq(ip,ith)
     h3(5,fmin:fmax)=-rsgtF(2,fmin:fmax)*srq(ip,ith)
     h3(6,fmin:fmax)=rsgtF(3,fmin:fmax)*srq2(ip,ith)
     
     !   Radial component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.
  elseif(trim(compo).eq."R") then
     h3(1,fmin:fmax)=rsgtF(5,fmin:fmax)
     h3(2,fmin:fmax)=-(rsgtF(7,fmin:fmax))*crq2(ip,ith)-rsgtF(6,fmin:fmax)
     h3(3,fmin:fmax)=(rsgtF(7,fmin:fmax))*crq2(ip,ith)-rsgtF(6,fmin:fmax)
     h3(4,fmin:fmax)=(rsgtF(9,fmin:fmax))*crq(ip,ith)
     h3(5,fmin:fmax)=(rsgtF(9,fmin:fmax))*srq(ip,ith)
     h3(6,fmin:fmax)=-(rsgtF(7,fmin:fmax))*srq2(ip,ith)

     !   Transverse component. Order: 1-rr, 2-tt, 3-pp, 4-rt, 5-rp, 6-tp.
  elseif(trim(compo).eq."T") then
     h3(1,fmin:fmax)=0.d0
     h3(2,fmin:fmax)=(rsgtF(8,fmin:fmax))*srq2(ip,ith)
     h3(3,fmin:fmax)=-h3(3,fmin:fmax)
     h3(4,fmin:fmax)=-(rsgtF(10,fmin:fmax))*srq(ip,ith)
     h3(5,fmin:fmax)=(rsgtF(10,fmin:fmax))*crq(ip,ith)
     h3(6,fmin:fmax)=-(rsgtF(8,fmin:fmax))*crq2(ip,ith)
  endif
  return
end subroutine rsgt2h3freq

subroutine tsgt2h4freq	
  use parameters
  use angles
  use tmpSGTs
  use kernels
  implicit none        

  !   This subroutine calculates the 36 independent elements of the 4th-order 
  !   TSGT from the 20 azimth-independent coefficients and the sines and 
  !   cosines of the azimuths. The 4th-order TSGT is then double-dotted with 
  !   the moment tensor to obtain a 6-element 2nd-order tensor.
	

  tmph01(fmin:fmax)=tsgtF(7,fmin:fmax)*cqs(ip,ith)
  tmph02(fmin:fmax)=tsgtF(7,fmin:fmax)*sqs(ip,ith)
  tmph03(fmin:fmax)=tsgtF(8,fmin:fmax)*csq(ip,ith)
  tmph04(fmin:fmax)=tsgtF(8,fmin:fmax)*ssq(ip,ith)
  tmph05(fmin:fmax)=tsgtF(9,fmin:fmax)*cqs2(ip,ith)
  tmph06(fmin:fmax)=tsgtF(10,fmin:fmax)*csq2(ip,ith)
  tmph07(fmin:fmax)=tsgtF(11,fmin:fmax)*csq2(ip,ith)
  tmph08(fmin:fmax)=tsgtF(11,fmin:fmax)*ssq2(ip,ith)
  tmph09(fmin:fmax)=tsgtF(12,fmin:fmax)*cqs2(ip,ith)
  tmph10(fmin:fmax)=tsgtF(12,fmin:fmax)*sqs2(ip,ith)
  tmph11(fmin:fmax)=tsgtF(13,fmin:fmax)*(csq(ip,ith)*cqs(ip,ith)-ssq(ip,ith)*sqs(ip,ith))
  tmph12(fmin:fmax)=tsgtF(13,fmin:fmax)*(ssq(ip,ith)*cqs(ip,ith)+csq(ip,ith)*sqs(ip,ith))
  tmph13(fmin:fmax)=tsgtF(14,fmin:fmax)*(csq(ip,ith)*cqs(ip,ith)+ssq(ip,ith)*sqs(ip,ith))
  tmph14(fmin:fmax)=tsgtF(14,fmin:fmax)*(ssq(ip,ith)*cqs(ip,ith)-csq(ip,ith)*sqs(ip,ith))
  tmph15(fmin:fmax)=tsgtF(15,fmin:fmax)*(csq2(ip,ith)*cqs(ip,ith)-ssq2(ip,ith)*sqs(ip,ith))
  tmph16(fmin:fmax)=tsgtF(15,fmin:fmax)*(ssq2(ip,ith)*cqs(ip,ith)+csq2(ip,ith)*sqs(ip,ith))
  tmph17(fmin:fmax)=tsgtF(17,fmin:fmax)*(csq(ip,ith)*cqs2(ip,ith)-ssq(ip,ith)*sqs2(ip,ith))
  tmph18(fmin:fmax)=tsgtF(17,fmin:fmax)*(ssq(ip,ith)*cqs2(ip,ith)+csq(ip,ith)*sqs2(ip,ith))
  tmph19(fmin:fmax)=tsgtF(19,fmin:fmax)*(csq2(ip,ith)*cqs2(ip,ith)-ssq2(ip,ith)*sqs2(ip,ith))
  tmph20(fmin:fmax)=tsgtF(19,fmin:fmax)*(ssq2(ip,ith)*cqs2(ip,ith)+csq2(ip,ith)*sqs2(ip,ith))
  tmph21(fmin:fmax)=tsgtF(20,fmin:fmax)*(csq2(ip,ith)*cqs2(ip,ith)+ssq2(ip,ith)*sqs2(ip,ith))
  tmph22(fmin:fmax)=tsgtF(20,fmin:fmax)*(ssq2(ip,ith)*cqs2(ip,ith)-csq2(ip,ith)*sqs2(ip,ith))
  tmph23(fmin:fmax)=tsgtF(16,fmin:fmax)*(csq2(ip,ith)*cqs(ip,ith)+ssq2(ip,ith)*sqs(ip,ith))
  tmph24(fmin:fmax)=tsgtF(16,fmin:fmax)*(ssq2(ip,ith)*cqs(ip,ith)-csq2(ip,ith)*sqs(ip,ith))
  tmph25(fmin:fmax)=tsgtF(18,fmin:fmax)*(csq(ip,ith)*cqs2(ip,ith)+ssq(ip,ith)*sqs2(ip,ith))
  tmph26(fmin:fmax)=tsgtF(18,fmin:fmax)*(ssq(ip,ith)*cqs2(ip,ith)-csq(ip,ith)*sqs2(ip,ith))
	
  !   rr component. 

  h4(1,fmin:fmax)=tsgtF(1,fmin:fmax)*mt(1)+(tmph05(fmin:fmax)-tsgtF(2,fmin:fmax))*mt(2)+&
       (-tmph05(fmin:fmax)-tsgtF(2,fmin:fmax))*mt(3) &
       -2.d0*tsgtF(5,fmin:fmax)*(cqs(ip,ith)*mt(4)+sqs(ip,ith)*mt(5)) &
       +2.d0*tsgtF(9,fmin:fmax)*sqs2(ip,ith)*mt(6)
	  
  !   tt component.
  
  h4(2,fmin:fmax)=(tmph06(fmin:fmax)-tsgtF(3,fmin:fmax))*mt(1) & 
       +(tmph21(fmin:fmax)+tmph19(fmin:fmax)-tmph07(fmin:fmax)-tmph09(fmin:fmax)+tsgtF(4,fmin:fmax))*mt(2) &
       +(-tmph21(fmin:fmax)-tmph19(fmin:fmax)-tmph07(fmin:fmax)+tmph09(fmin:fmax)+tsgtF(4,fmin:fmax))*mt(3) &
       +2.d0*(tmph23(fmin:fmax)-tmph15(fmin:fmax)+tmph01(fmin:fmax))*mt(4) &
       +2.d0*(-tmph24(fmin:fmax)-tmph16(fmin:fmax)+tmph02(fmin:fmax))*mt(5) &
       +2.d0*(-tmph22(fmin:fmax)+tmph20(fmin:fmax)-tmph10(fmin:fmax))*mt(6)
     
  !   pp component. 

  h4(3,fmin:fmax)=(-tmph06(fmin:fmax)-tsgtF(3,fmin:fmax))*mt(1) &
       +(-tmph21(fmin:fmax)-tmph19(fmin:fmax)+tmph07(fmin:fmax)-tmph09(fmin:fmax)+tsgtF(4,fmin:fmax))*mt(2) &
       +(tmph21(fmin:fmax)+tmph19(fmin:fmax)+tmph07(fmin:fmax)+tmph09(fmin:fmax)+tsgtF(4,fmin:fmax))*mt(3) &
       +2.d0*(-tmph23(fmin:fmax)+tmph15(fmin:fmax)+tmph01(fmin:fmax))*mt(4) &
       +2.d0*(tmph24(fmin:fmax)+tmph16(fmin:fmax)+tmph02(fmin:fmax))*mt(5) &
       +2.d0*(tmph22(fmin:fmax)-tmph20(fmin:fmax)-tmph10(fmin:fmax))*mt(6)
     
  !   rt component. 

  h4(4,fmin:fmax)=-tsgtF(6,fmin:fmax)*csq(ip,ith)*mt(1) &
       +(tmph25(fmin:fmax)-tmph17(fmin:fmax)+tmph03(fmin:fmax))*mt(2) &
       +(-tmph25(fmin:fmax)+tmph17(fmin:fmax)+tmph03(fmin:fmax))*mt(3) &
       +2.d0*(-tmph13(fmin:fmax)+tmph11(fmin:fmax))*mt(4) &
       +2.d0*(tmph14(fmin:fmax)+tmph12(fmin:fmax))*mt(5) &
       +2.d0*(-tmph26(fmin:fmax)-tmph18(fmin:fmax))*mt(6)
     
  !   rp component 

  h4(5,fmin:fmax)=(-tsgtF(6,fmin:fmax)*ssq(ip,ith))*mt(1) &
       +(tmph26(fmin:fmax)-tmph18(fmin:fmax)+tmph04(fmin:fmax))*mt(2) &
       +(-tmph26(fmin:fmax)+tmph18(fmin:fmax)+tmph04(fmin:fmax))*mt(3) &
       +2.d0*(-tmph14(fmin:fmax)+tmph12(fmin:fmax))*mt(4) &
       +2.d0*(-tmph13(fmin:fmax)-tmph11(fmin:fmax))*mt(5) &
       +2.d0*(tmph25(fmin:fmax)+tmph17(fmin:fmax))*mt(6)
     
  !   tp component. 

  h4(6,fmin:fmax)=tsgtF(10,fmin:fmax)*ssq2(ip,ith)*mt(1) &
       +(tmph22(fmin:fmax)+tmph20(fmin:fmax)-tmph08(fmin:fmax))*mt(2) &
       +(-tmph22(fmin:fmax)-tmph20(fmin:fmax)-tmph08(fmin:fmax))*mt(3) &
       +2.d0*(tmph24(fmin:fmax)-tmph16(fmin:fmax))*mt(4) &
       +2.d0*(tmph23(fmin:fmax)+tmph15(fmin:fmax))*mt(5) &
       +2.d0*(tmph21(fmin:fmax)-tmph19(fmin:fmax))*mt(6)

  return
end subroutine tsgt2h4freq


subroutine vectorFFT_double(imin,imax,np1,ccvec,rvec,omegai,tlen,iWindowStart,iWindowEnd)
  ! this subroutine particularly calculates the FFT of the given tensor and make a double tensor
  
  implicit none
  integer :: iWindowStart,iWindowEnd
  integer :: i,imin,imax,np1,n1,m1
  complex(kind(0d0)) :: ccvec(imin:imax)
  complex(kind(0d0)) :: cvec(0:2*np1-1)
  real(kind(0d0)) :: rvec(iWindowStart:iWindowEnd)
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai, tlen,samplingHz
  
  cvec = cmplx(0.d0)
  cvec(imin:imax)=ccvec(imin:imax)

  samplingHz = dble(2*np1)/tlen
  
  do i = imin, np1-1
     n1 = np1 +i
     m1 = np1 -i
     cvec(n1) = conjg(cvec(m1))
  enddo
  

  call cdft(4*np1,cos(pi/(2*np1)),sin(pi/(2*np1)), cvec(0:2*np1-1))
  do i = iWindowStart, iWindowEnd
     rvec(i) = dble(dble(cvec(i))*dble(exp(omegai*dble(i)/samplingHz))/tlen*1.d3)
  enddo

  
  
  return
end subroutine vectorFFT_double

subroutine tensorFFT_double(n,imin,imax,np1,ccvec,rvec,omegai,tlen,iWindowStart,iWindowEnd)
  ! this subroutine particularly calculates the FFT of the given tensor and make a double tensor
  
  implicit none
  integer :: iWindowStart,iWindowEnd
  integer :: i,j,n,imin,imax,np1,n1,m1
  complex(kind(0d0)) :: ccvec(1:n,imin:imax)
  complex(kind(0d0)) :: cvec(1:n,0:2*np1-1)
  real(kind(0d0)) :: rvec(1:n,iWindowStart:iWindowEnd)
  real(kind(0d0)), parameter :: pi = 3.141592653589793d0
  real(kind(0d0)) :: omegai, tlen,samplingHz
  
  cvec = cmplx(0.d0)
  cvec(1:n,imin:imax)=ccvec(1:n,imin:imax)

  samplingHz = dble(2*np1)/tlen
  do j = 1,n
     do i = imin, np1-1
        n1 = np1 +i
        m1 = np1 -i
        cvec(j,n1) = conjg(cvec(j,m1))
     enddo
  enddo
  
  
  
  
  do j = 1,n
     call cdft(4*np1,cos(pi/(2*np1)),sin(pi/(2*np1)), cvec(j,0:2*np1-1))
     do i = iWindowStart, iWindowEnd
        rvec(j,i) = dble(dble(cvec(j,i))*dble(exp(omegai*dble(i)/samplingHz))/tlen*1.d3)
        
     enddo
  enddo
  
  
  return
end subroutine tensorFFT_double


subroutine cdft(n, wr, wi, c)
      
  integer :: n, i, j, k, l, m
  real(kind(0d0)) :: wr, wi, a(0 : n - 1), wmr, wmi, wkr, wki 
  real(kind(0d0)) ::wdr, wdi, ss, xr, xi
  complex(kind(0d0)) :: c(0:n/2-1)
  
  do i = 0, n/2-1
     a(2*i) = dble(c(i))
     a(2*i+1) = imag(c(i))
  enddo
  

  wmr = wr
  wmi = wi
  m = n
  do while (m .gt. 4)
     l = m / 2
     wkr = 1
     wki = 0
     wdr = 1 - 2 * wmi * wmi
     wdi = 2 * wmi * wmr
     ss = 2 * wdi
     wmr = wdr
     wmi = wdi
     do j = 0, n - m, m
        i = j + l
        xr = a(j) - a(i)
        xi = a(j + 1) - a(i + 1)
        a(j) = a(j) + a(i)
        a(j + 1) = a(j + 1) + a(i + 1)
        a(i) = xr
        a(i + 1) = xi
        xr = a(j + 2) - a(i + 2)
        xi = a(j + 3) - a(i + 3)
        a(j + 2) = a(j + 2) + a(i + 2)
        a(j + 3) = a(j + 3) + a(i + 3)
        a(i + 2) = wdr * xr - wdi * xi
        a(i + 3) = wdr * xi + wdi * xr
     enddo
     do k = 4, l - 4, 4
        wkr = wkr - ss * wdi
        wki = wki + ss * wdr
        wdr = wdr - ss * wki
        wdi = wdi + ss * wkr
        do j = k, n - m + k, m
           i = j + l
           xr = a(j) - a(i)
           xi = a(j + 1) - a(i + 1)
           a(j) = a(j) + a(i)
           a(j + 1) = a(j + 1) + a(i + 1)
           a(i) = wkr * xr - wki * xi
           a(i + 1) = wkr * xi + wki * xr
           xr = a(j + 2) - a(i + 2)
           xi = a(j + 3) - a(i + 3)
           a(j + 2) = a(j + 2) + a(i + 2)
           a(j + 3) = a(j + 3) + a(i + 3)
           a(i + 2) = wdr * xr - wdi * xi
           a(i + 3) = wdr * xi + wdi * xr
        enddo
     enddo
     m = l
  enddo
  if (m .gt. 2) then
     do j = 0, n - 4, 4
        xr = a(j) - a(j + 2)
        xi = a(j + 1) - a(j + 3)
        a(j) = a(j) + a(j + 2)
        a(j + 1) = a(j + 1) + a(j + 3)
        a(j + 2) = xr
        a(j + 3) = xi
     enddo
  endif
  if (n .gt. 4) call bitrv2(n, a)

  
  do i = 0, n/2-1
     c(i) = dcmplx(a(2*i), a(2*i+1))
  enddo
  
  
end subroutine cdft



subroutine bitrv2(n, a)
  integer :: n, j, j1, k, k1, l, m, m2, n2
  real(kind(0d0)) :: a(0 : n - 1), xr, xi
  
  m = n / 4
  m2 = 2 * m
  n2 = n - 2
  k = 0
  do j = 0, m2 - 4, 4
     if (j .lt. k) then
        xr = a(j)
        xi = a(j + 1)
        a(j) = a(k)
        a(j + 1) = a(k + 1)
        a(k) = xr
        a(k + 1) = xi
     else if (j .gt. k) then
        j1 = n2 - j
        k1 = n2 - k
        xr = a(j1)
        xi = a(j1 + 1)
        a(j1) = a(k1)
        a(j1 + 1) = a(k1 + 1)
        a(k1) = xr
        a(k1 + 1) = xi
     endif
     k1 = m2 + k
     xr = a(j + 2)
     xi = a(j + 3)
     a(j + 2) = a(k1)
     a(j + 3) = a(k1 + 1)
     a(k1) = xr
     a(k1 + 1) = xi
     l = m
     do while (k .ge. l)
        k = k - l
        l = l / 2
     enddo
     k = k + l
  enddo
end subroutine bitrv2
