program combinekernels
    implicit none

    !==== VARIABLE DECLARATIONS ====

    !input/output file names:
    character(100),parameter:: gridfile = 'kernels/tmpvideo/griddata'
    character(100),parameter:: videoprefix = 'kernels/tmpvideo/videokernel.'
    character(100) tmpchar,infile,outfile

    !grid information parameters:
    integer nr,nphi,ntheta,nkvtype,nfilter,ntimes_seismo,ntimes_video
    real mtensor(6),tstart,tend
    real,allocatable,dimension(:)   :: radii,displacement
    real,allocatable,dimension(:,:) :: phis,thetas

    !2D kernel all timestep, all types, all filters (DSM output)
    real, allocatable, dimension(:,:,:,:) :: dsm_kernel

    !3D kernel for a single timestep:
    real, allocatable, dimension(:,:,:) :: kernel

    !local variables
    integer ir,itheta,iphi,it,itime

    !==== READ GRID DATA ====
    open(1,file=trim(gridfile),action='read',status='old',form='unformatted',access='sequential')

    read(1) nr,nphi,ntheta,nkvtype,nfilter,ntimes_seismo,ntimes_video
    read(1) tstart,tend
    read(1) mtensor

    allocate(radii(nr))
    allocate(displacement(ntimes_seismo))
    allocate(phis(nphi,ntheta))
    allocate(thetas(nphi,ntheta))
    allocate(dsm_kernel(nphi,nkvtype+1,nfilter+1,ntimes_video))
    allocate(kernel(nphi,ntheta,nr))

    read(1) radii
    read(1) phis
    read(1) thetas
    read(1) displacement

    close(1)

    print *,'found file with grid dimensions:',nr,nphi,ntheta
    print *,'number of timesteps::',ntimes_seismo,ntimes_video
    print *,'from times:',tstart,tend
    print *,'rmin,rmax:',radii(1),radii(nr)
    print *,'moment tensor:',mtensor

    !==== EXTRACT SINGLE TIMESTEP FOR EVERY RADIUS AND WRITE 3D KERNEL ====
    do itime=1,ntimes_video
      !read files for single r and theta
      do ir=1,nr
        do it =1,ntheta
          write(tmpchar,'(I0.3,".",I0.3)') ir,it
          infile = trim(videoprefix)//trim(tmpchar)
          open(1,file=trim(infile),status='old',form='unformatted',access='direct',&
                            recl=4*nphi*(1+nkvtype)*(1+nfilter)*ntimes_video)
          read(1,rec=1) dsm_kernel
          close(1)
          kernel(:,it,ir) = dsm_kernel(:,2,1,itime)
        enddo
      enddo
      !write 3D kernel for single timestep
      tmpchar(:) = ' '
      write(tmpchar,'(I0.3)') int(itime)
      print *,'writing kernel: ',trim(tmpchar)
      outfile = 'kernels/movie/kernel3d.'//trim(adjustl(tmpchar))
      open(1,file=trim(outfile),status='new',form='unformatted',access='direct',recl=4*nphi*ntheta*nr)
      write(1,rec=1) kernel
      close(1)
      print *,'kernel norm',itime,sum(kernel(:,:,:)*kernel(:,:,:))
    enddo

    print *,'ALL DONE'

end program
