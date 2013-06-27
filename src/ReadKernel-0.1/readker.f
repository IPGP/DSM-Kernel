        80*character infile,eventid,staid
	real, dimension (:), allocatable :: r0,r,filt,xmaxt
	real, dimension (:,:,:,:), allocatable :: ker
	real, dimension (:,:), allocatable :: u,v,hu,phi,theta,phig,thetag,u0
	real xm(7),mt(6),twin(4)

	real, dimension (:), allocatable :: totaltraveltime

	pi=4.*atan(1.)
	
	write(*,'(a)') 'Input kernel file name: '
	read(*,'(a)') infile


        open(1,file=infile,status='old',form='unformatted',
     1    access='sequential')
	
        read(1) nr,nphi,ntheta,nc,kc,mtype,nktype,nfilter,ift
	read(1) dt,twin,t1,t2
	read(1) fclp,fchp,mb,mid,mbmax
	read(1) slat,slon,sdep,sdep0,rlat,rlon,delta
	write(*,'(/a,3(2x,f8.3))') 'Source location: ',slat,slon,sdep
	write(*,'(a,2(2x,f8.3))') 'Station location: ',rlat,rlon
	write(*,'(a,2x,f7.2)') 'Epicentral distance: ',delta
	read(1) xm
	read(1) mt
	print*, nktype


	allocate(r0(nr),r(nr),phi(nphi,ntheta),theta(nphi,ntheta))
	allocate(phig(nphi,ntheta),thetag(nphi,ntheta))
	allocate(ker(nr,nphi,ntheta,0:nktype),xmaxt(nktype))
	allocate(u(3,0:nc),v(3,0:nc),hu(3,0:nc),filt(0:mbmax),u0(3,0:nc))
	
	allocate(totaltraveltime(0:nr))

	read(1) r0
	read(1) r
	read(1) phi
	read(1) theta
	read(1) ker
	read(1) filt
	read(1) u
	read(1) u0
	read(1) v
	read(1) hu
	read(1) phig
	read(1) thetag
	read(1) ncevt,ncsta
	read(1) eventid(1:ncevt),staid(1:ncsta)
	close(1)
	
	write(*,'(/a,3(1x,i4)/)') 'Info from kernel file (nr, nphi, ntheta):',nr,nphi,ntheta

	write(*,'(a)') 'View (1-Mapview, 2-Source-Receiver plane, 3-Transverse plane): '
	read(*,*) nview
	if(nview.eq.1) then
	  write(*,'(a/7(1x,i3,1x,f7.2))') 'Radii of radial samples (in km):',((ir,r(ir)),ir=1,nr)
	  write(*,'(a,i3)') 'Choose horizontal slice between 1 and: ',nr
	  read(*,*) k
	elseif(nview.eq.2) then
	  npr=ntheta
	  write(*,'(a/7(1x,i3,1x,f7.2))') 'Theta of slices:',((ith,theta(nphi/2,ith)),ith=1,ntheta)
	  write(*,'(a,i3)') 'Choose slice between 1 and ',ntheta
	  read(*,*) k
	else
	  write(*,'(a/7(1x,i3,1x,f7.2))') 'Phi of slices:',((ip,phi(ip,ntheta/2)),ip=1,nphi)
 	  write(*,'(a,i3)') 'Choose profile between 1 and ',nphi
	  read(*,*) k
	endif∑

	ker(:,:,:,:) =ker(:,:,:,:)*1.e9
	!ker(:,:,:,5) = ker(:,:,:,5)*1.e9!*1.e27!/gnormt!*1.e-3
	!ker(:,:,:,1) = ker(:,:,:,1)*1.e9!*1.e27!/gnormt!*1.e-3
	!ker(:,:,:,3) = ker(:,:,:,3)*1.e9!*1.e27!/gnorma

	totaltraveltime=0.d0
	do ir = 1,nr
	   do ith=1,ntheta
	      do ip=1,nphi
		 totaltraveltime(ir)=totaltraveltime(ir)+ker(ir,ip,ith,1)!*40.d0
		 totaltraveltime(0)=totaltraveltime(0)+ker(ir,ip,ith,1)*40.d0
	      enddo
	   enddo
	enddo

	open (1,file='traveltime',status='unknown')
	do ir=1,nr
	   write(1,*) r(ir),totaltraveltime(ir)
	enddo
	print*, totaltraveltime(0)
	close(1)

	
	ker(ir,1:nphi,ntheta/2+2,1:4)=ker(ir,1:nphi,ntheta/2+1,1:4)

	

	do ip=1,nphi
	  do ith=1,ntheta
	    if(phig(ip,ith).gt.360.) phig(ip,ith)=phig(ip,ith)-360.
	    thetag(ip,ith)=90.-thetag(ip,ith)
	    
	    ! it isn't just normal for -side which doesn't work at all for this moment 
	    ! (id gotta verify this but since just for the revision of FCZGK2012, with a explosive source,
	    !  we can do this)

	    if(ith.lt.ntheta/2+1) then
	       do ir=nr,1,-1
		  ker(ir,ip,ith,1:4) = ker(ir,ip,ntheta-ith+1,1:4)
	       enddo
	    endif
	    

	    !!!!!
	    
	    
	  end do
	end do



	do i=1,nktype
	  xmaxt(i)=-1.d10
	  do ip=1,nphi
	    do ith=1,ntheta
	      if(abs(ker(k,ip,ith,i)).gt.xmaxt(i)) xmaxt(i)=abs(ker(k,ip,ith,i))
	    end do
	  end do
	end do
	
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
	       if((phig(ip,ith).gt.180).and.(phig(ip,ith).le.270.)) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=1,4)

	        if(phig(ip,ith).gt.270.) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=1,4)


	       if(phig(ip,ith).le.180.) write(1,'(f9.4,1x,f9.4,4(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=1,4)
	      if(phig(ip,ith).gt.180.) write(2,'(f9.4,1x,f9.4,4(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=5,8)
	      if(phig(ip,ith).le.180.) write(2,'(f9.4,1x,f9.4,4(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=5,8)
	      if(phig(ip,ith).gt.180.) write(3,'(f9.4,1x,f9.4,10(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith)-360.,(ker(k,ip,ith,j),j=9,18)
	      if(phig(ip,ith).le.180.) write(3,'(f9.4,1x,f9.4,10(1x,e12.5))') 
     1	        thetag(ip,ith),phig(ip,ith),(ker(k,ip,ith,j),j=9,18)
	    end do
	  end do
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
	      if(phi(ip,k).gt.180) write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),phi(ip,k)-360.,
     1	        theta(ip,k),(ker(ir,ip,k,j),j=1,1),ker(ir,ip,k,3), ker(ir,ip,k,4)
	      if(phi(ip,k).le.180) write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),phi(ip,k),
     1	        theta(ip,k),(ker(ir,ip,k,j),j=1,1),ker(ir,ip,k,3), ker(ir,ip,k,4)
	      !(ker(ir,ip,k,j),j=1)
	      write(2,'(3(1x,f9.4),4(1x,e15.8))') r(ir),phi(ip,k),
     1	        theta(ip,k),(ker(ir,ip,k,j),j=5,8)
!     1	        theta(ip,ith),ker(ir,ip,k,5),ker(ir,ip,k,7),
!     1          ker(ir,ip,k,8)

	      write(3,'(3(1x,f9.4),10(1x,e15.8))') r(ir),phi(ip,k),
     1	        theta(ip,k),(ker(ir,ip,k,j),j=9,18)
	    end do
	  end do
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
	       write(1,'(3(1x,f9.4),3(1x,e15.8))') r(ir),theta(k,ith),
     1	        phi(k,ith),(ker(ir,k,ith,j),j=1,1),(ker(ir,k,ith,j),j=3,4)
	      write(2,'(3(1x,f9.4),4(1x,e15.8))') r(ir),theta(k,ith),
     1	        phi(k,ith),ker(ir,k,ith,5),ker(ir,k,ith,7),
     1          ker(ir,k,ith,8)
	      write(3,'(3(1x,f9.4),10(1x,e15.8))') r(ir),theta(k,ith),
     1	        phi(k,ith),(ker(ir,k,ith,j),j=9,18)
	    end do
	  end do
	  close(1)
	  close(2)
	  close(3)
	endif

	write(*,'(a,14(1x,e9.2))') 'Maxima: ',xmaxt

	stop
	end
