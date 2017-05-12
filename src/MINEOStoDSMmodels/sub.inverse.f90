  !    Subroutine for inversing a regular matrix
  !
  !     gfortran main.f sub.inverse.f -llapack
  !
  !     inverse(n,a,nmax,ia)
  !             n : (in) dimension of matrix 'a'
  !        a(*,*) : (in) n*n matrix
  !          nmax : (in) working space
  !       ia(*,*) : (out) inv(a)
  !
  !************************* Subroutine using LAPACK ****************************
  !
  !     DGETRF(M,N,A,LDA,IPIV,INFO)
  !             M : number of lines
  !             N : number of columns
  !      A(LDA,*) : input->matrix A(M,N), output->[L]&[U]
  !                 (A=[P]*[L]*[U],[P]:permulation matrix)
  !           LDA : size of 1st dimension of A. LDA=M as usual (so don't care)
  !       IPIV(*) : array. * > min(M,N)
  !          INFO : if successful, return 0
  !
  !     DGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
  !             N : dimension of A(N,N)
  !      A(LDA,*) : input->N*N matrix after LU deconposition, output->inv(A) 
  !           LDA : N as usual
  !       IPIV(*) : * > N
  !   WORK(LWORK) : working array.
  !         LWORK : size of WORK. =N,ok
  !          INFO : if successful return 0
  !
  !******************************************************************************
  

subroutine inverseLU(nmax,a,ia)
  implicit none
  integer n,nmax,i,j
  double precision a(nmax,nmax),ia(nmax,nmax),work(nmax)
  integer ipiv(nmax),info
  n=nmax
  
  !     LU deconposition
  !     keep a copy of 'a' because 'dgetrf' overwrites 'a' 

  do i = 1,n
     do j = 1,n
        ia(i,j) = a(i,j)
     enddo
  enddo

  !     1st, LU deconposition of 'ia' using 'dgetrf'

  call dgetrf(n,n,ia,nmax,ipiv,info)

  !    if info < 0: some probrem about 'a', make sure of type of 'a'
  !                                                     or size of arrays

  if(info .lt. 0) then
     write(*,*) "LU decomposition error!"
     write(*,*) "maybe wrong matrix input (+_+;)"
     stop
  endif
  
  !     info > 0: division by 0. input matrix 'a' is not regular.
  
  if(info .gt. 0) then
     write(*,*) "LU decomposition error!"
     write(*,*) "Input matrix is not regular."
     stop
  endif

  !     if ever info = 0, 
  !    in the case of |ia(i,i)|~0, it's regarded as division by 0
  
  do i = 1,n
     if(abs(ia(i,i)) .lt. 1.D-14) then
        write(*,*) "LU decomposition error!"
        write(*,*) "Input matrix may be not regular."
        write(*,*) "it has eigen values near 0."
        stop
     endif
  enddo
  
  !     after LU decomposition,inverse 'ia' by using 'dgetri'
  
  call dgetri(n,ia,nmax,ipiv,work,n,info)
  
  !     in fact error handling should be finished at the step of LU decomposition
  
  if(info .ne. 0) then
     write(*,*) "inversion error!"
     write(*,*) "maybe serious problems happen."
     stop
  endif
  
  return
end subroutine inverseLU

