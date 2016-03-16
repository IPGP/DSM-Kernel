program ReadKernel
  
  use parameters ! use the same module as in KernelMaker

  implicit none
  

  integer :: ift,i,j 
  character(200) :: kertotalfile,kernelp,kernels,kernela
  character(120) :: tmpfile,kernelinf,dummy
  integer :: nphi,ntheta,nc,idum,nview,ir,ith,ip,k,npr
  real(kind(0e0)) :: fdum,slat,slon,sdep,sdep0,rlat,rlon,delta,criteria

  real(kind(0e0)), dimension (:), allocatable :: filt,xmaxt,rtmp,r0tmp
  real(kind(0e0)), dimension (:,:,:,:), allocatable :: ker
  real(kind(0e0)), dimension (:,:), allocatable :: u,v,hu,phi,theta,phig,thetag,u0

  ! for MPI ! but we don't use this for this version
  ! integer :: nproc, my_rank, ierr

  
  

  tmpfile='tmpworkingfile_for_ReadKernel'

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
  read(1,110) kernelinf
  read(1,*) nview
  read(1,*) criteria
  close(1)



  ! Inputting parameters
  call pinputKernel(kernelinf)
  ! DSM constants calculation
  

  

  ift = 0
  
  kertotalfile = trim(parentDir)//"/"//trim(stationName)//"."//trim(eventName)//"."//trim(phase)//"."//trim(compo)//"."//trim(freqid(ift))
  open(1,file=kertotalfile,status='old',form='unformatted',access='sequential')
  !read(1) nr,nphi,ntheta,nc,kc,mtype,nktype,nfilter,ift
  read(1) nr,nphi,ntheta,nc,idum,mtype,idum,idum,ift
  !read(1) dt,twin,t1,t2
  read(1) fdum,fdum,fdum,fdum
  !read(1) fclp(ift),fchp(ift),mb,mid,mbmax
  read(1) fdum,fdum,idum,idum,idum 
  !read(1) slat,slon,sdep,sdep0,rlat,rlon,delta 
  read(1) slat,slon,sdep,sdep0,rlat,rlon,delta
  write(*,'(/a,3(2x,f8.3))') 'Source location: ',slat,slon,sdep
  write(*,'(a,2(2x,f8.3))') 'Station location: ',rlat,rlon
  write(*,'(a,2x,f7.2)') 'Epicentral distance: ',delta



  !read(1) xm
  read(1) fdum,fdum,fdum,fdum,fdum,fdum
  !read(1) mt
  read(1) fdum,fdum,fdum,fdum,fdum,fdum
  
  
  allocate(phi(nphi,ntheta),theta(nphi,ntheta))
  allocate(phig(nphi,ntheta),thetag(nphi,ntheta))
  allocate(ker(nr,nphi,ntheta,0:nktype),xmaxt(nktype))
  allocate(u(3,0:nc),v(3,0:nc),hu(3,0:nc),filt(0:nc),u0(3,0:nc))
  r0_n =  int((r0max-r0min)/r0delta)+1
  allocate(r0D(1:r0_n))
  allocate(r(1:nr))
  allocate(r0tmp(1:r0_n))
  allocate(rtmp(1:nr))
  
  
  read(1) r0tmp
  read(1) rtmp
  r(:)=rtmp(:)
  r0D=r0tmp(:)
  


  read(1) phi
  read(1) theta
  read(1) ker
  !read(1) filt
  !read(1) u
  !read(1) u0
  !read(1) v
  !read(1) hu
  read(1) phig
  read(1) thetag
  !read(1) ncevt,ncsta
  !read(1) eventid(1:ncevt),staid(1:ncsta)
  close(1)
  
  !write(*,'(/a,3(1x,i4)/)') 'Info from kernel file (nr, nphi, ntheta):',nr,nphi,ntheta
  print *, 'Info from kernel file (nr, nphi, ntheta):',nr,nphi,ntheta

  print *, 'View (1-Mapview, 2-Source-Receiver plane, 3-Transverse plane): '
  !read(*,*) nview
  !read *, nview
  if(nview.eq.1) then
     !write(*,'(a/7(1x,i3,1x,f7.2))') 'Radii of radial samples (in km):',((ir,r(ir)),ir=1,nr)
     !write(*,'(a,i3)') 'Choose horizontal slice between 1 and: ',nr
     !read(*,*) k
     do ir = 1,nr-1
        if(r(ir).eq.criteria) then
           print *, r(ir)
           k=ir
        endif
     enddo     
  elseif(nview.eq.2) then
     npr=ntheta
     !write(*,'(a/7(1x,i3,1x,f7.2))') 'Theta of slices:',((ith,theta(nphi/2,ith)),ith=1,ntheta)
     !write(*,'(a,i3)') 'Choose slice between 1 and ',ntheta
     do ith=1,ntheta 
        if(theta(nphi/2,ith).eq.criteria) then
           k=ith
        endif
     enddo
  
  else
     !write(*,'(a/7(1x,i3,1x,f7.2))') 'Phi of slices:',((ip,phi(ip,ntheta/2)),ip=1,nphi)
     !write(*,'(a,i3)') 'Choose profile between 1 and ',nphi
     !read(*,*) k
     do ip=1,nphi
        if(phi(ip,ntheta/2).eq.criteria) then
           k=ip
        endif
     enddo
  endif

  ker(:,:,:,:) =ker(:,:,:,:)*1.e9

  do ip=1,nphi
     do ith=1,ntheta
        if(phig(ip,ith).gt.360.) phig(ip,ith)=phig(ip,ith)-360.
        thetag(ip,ith)=90.-thetag(ip,ith)
        
     enddo
  enddo

  
  
  do i=1,nktype
     xmaxt(i)=-1.d10
     do ip=1,nphi
        do ith=1,ntheta
           if(abs(ker(k,ip,ith,i)).gt.xmaxt(i)) xmaxt(i)=abs(ker(k,ip,ith,i))
        enddo
     enddo
  enddo

  if(nview.eq.1) then
     open(1,file='coordrs',status='unknown')
     write(1,*) nview
     write(1,*) nphi
     write(1,*) ntheta
     write(1,*) 6371.-sdep
     write(1,*) delta
     close(1)
     open(1,file='kernelp',status='unknown')
     open(2,file='kernels',status='unknown')
     open(3,file='kernela',status='unknown')
     do ith=ntheta,1,-1
        do ip=1,nphi
           !print *, ip,ith,k,j
           !	       if(phig(ip,ith).gt.180.) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') 
           !     1	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=1,4)
           if((phig(ip,ith).gt.180).and.(phig(ip,ith).le.270.)) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') &
                thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=1,4)
           
           if(phig(ip,ith).gt.270.) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') &
     	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=1,4)
           if(phig(ip,ith).le.180.) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))')  &
     	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=1,4)
           if(phig(ip,ith).gt.180.) write(2,'(f9.4,1x,f9.4,4(1x,e12.5))')  &
     	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=5,8)
           if(phig(ip,ith).le.180.) write(2,'(f9.4,1x,f9.4,4(1x,e12.5))') &
     	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=5,8)
           !if(phig(ip,ith).gt.180.) write(3,'(f9.4,1x,f9.4,10(1x,e12.5))') &
     	   !     thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=9,18)
           !if(phig(ip,ith).le.180.) write(3,'(f9.4,1x,f9.4,10(1x,e12.5))') &
     	   !     thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=9,18)
        enddo
     enddo
     close(1)
     close(2)
     close(3)
  elseif(nview.eq.2) then
     open(1,file='coordrs',status='unknown')
     write(1,*) nview
     write(1,*) nphi
     write(1,*) nr
     write(1,*) 6371.-sdep
     write(1,*) delta
     close(1)
     open(1,file='kernelp',status='unknown')
     open(2,file='kernels',status='unknown')
     open(3,file='kernela',status='unknown')
     open(4,file='kernelw1',status='unknown')
     open(5,file='kernelw2',status='unknown')
     open(6,file='kernelw3',status='unknown')
     do ir=nr,1,-1
        do ip=1,nphi
           if(phi(ip,k).gt.180) write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),phi(ip,k)-360., &
     	        theta(ip,k),(ker(ir,ip,k,j),j=1,1),ker(ir,ip,k,3), ker(ir,ip,k,4)
           if(phi(ip,k).le.180) write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),phi(ip,k), &
     	        theta(ip,k),(ker(ir,ip,k,j),j=1,1),ker(ir,ip,k,3), ker(ir,ip,k,4)
           !(ker(ir,ip,k,j),j=1)
           write(2,'(3(1x,f9.4),4(1x,e15.8))') r(ir),phi(ip,k), &
     	        theta(ip,k),(ker(ir,ip,k,j),j=5,8)
           !     1	        theta(ip,ith),ker(ir,ip,k,5),ker(ir,ip,k,7),
           !     1          ker(ir,ip,k,8)
           
           !write(3,'(3(1x,f9.4),10(1x,e15.8))') r(ir),phi(ip,k),&
     	   !     theta(ip,k),(ker(ir,ip,k,j),j=9,18)
        enddo
     enddo
     close(1)
     close(2)
     close(3)
     close(4)
  else
     open(1,file='coordrs',status='unknown')
     write(1,*) nview
     write(1,*) ntheta
     write(1,*) nr
     write(1,*) 6371.-sdep
     write(1,*) delta
     close(1)
     open(1,file='kernelp',status='unknown')
     open(2,file='kernels',status='unknown')
     open(3,file='kernela',status='unknown')
     do ir=nr,1,-1
        do ith=1,ntheta
           write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),theta(k,ith),&
  	        phi(k,ith),(ker(ir,k,ith,j),j=1,1),(ker(ir,k,ith,j),j=3,4)
           write(2,'(3(1x,f9.4),4(1x,e15.8))') r(ir),theta(k,ith), &
     	        phi(k,ith),ker(ir,k,ith,5),ker(ir,k,ith,7), &
                ker(ir,k,ith,8)
           !write(3,'(3(1x,f9.4),10(1x,e15.8))') r(ir),theta(k,ith), &
     	   !     phi(k,ith),(ker(ir,k,ith,j),j=9,18)
        enddo
     enddo
     close(1)
     close(2)
     close(3)
  endif
  
  write(*,'(a,14(1x,e9.2))') 'Maxima: ',xmaxt
  stop
end program ReadKernel
