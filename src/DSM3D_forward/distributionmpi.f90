
!  MPI coloring
!                                      2016.03. FUJI Nobuaki
!
!           Copyright Institut de Physique du Globe de Paris

! module paramsmpi
! subroutine mpicoloring


module paramsmpi
  
  implicit none
  
  ! for MPI

  integer :: mpii
  integer :: nproc,my_rank,ierr
  integer :: petot, my_rank, ierr
  integer :: filenum, mpios
  integer :: outputmemory   ! MB
  integer :: outputinterval
  real(kind(0d0)) :: memoryperomega ! MB


  ! NF does not fully understand KK's comments below:
  !     when the values to be output use memory over outputmemory MB,
  !     they are written in output files. The interval is outputinterval
  !     memoryperomega is the quantity of memory used for one omega step


  
  

end module pramsmpi



subroutine mpicoloring
  
  use paramsmpi
  use variables

  implicit none

  ! pinput will read number of colors (not yet modified NF!!!!)
  !  and here nb_colors*nb_keys should be nbbigproc !!!
  !
  ! definition of Ifrq (which color calculates which frequencies)
  ! Here I consider a trapezoid feature of lmax with respect to frequency (i)
  !
  nb_keys=nbbigproc/nb_colors
  if(nb_keys*nb_colors.ne.nbbigproc) then
     print *, "It is mendatory to set nb_colors to be a divisor of nbbigproc"
     stop
  endif
  nb_freq_color = (imax-imin+1)/nb_colors
  
  allocate(key(0:nbbigproc-1))
  allocate(color(0:nbbigproc-1))
  allocate(Ifrq(0:nb_colors-1,0:nb_freq_color-1))

  if(mybigrank.eq.0) then
     Ifrq(:,:)=0
     do i=0,nbbigproc-1
        color(i)=i/nb_keys
        key(i)=mod(i,nb_keys)
     enddo
     do j=0,nb_colors-1
        l=0
        do i=imin,imax
           if((i.ne.0).and.((mod(imax-i_color-i,2*nproc).eq.0).or.(mod(imax+i_color+1-i,2*nproc).eq.0))) then
              Ifrq(j,l)=i
              l=l+1
           endif
        enddo
     enddo
  endif



end subroutine mpicoloring
