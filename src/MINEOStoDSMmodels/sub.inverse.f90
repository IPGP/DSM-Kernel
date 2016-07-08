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
  
  
subroutine svdinverse(m,n,a,ia,LWORK,INFO)
  implicit none
  
  ! m.le.n so U is small and VT is big
  
  integer m,n,k,j,i
  double precision :: A(m,n),IA(n,m),tmpA(n,n),origin(m,n)
  integer :: INFO, LDA, LDU, LDVT, LWORK
  !parameter(LWORK = 3*n+m)
  double precision :: S(m),B(m,m),VT(n,n)
  double precision :: maxS
  double precision :: work(LWORK)
  double precision, parameter :: eps = 1.d-9
  
  integer :: ii
  
  ! LWORK=5*m for this configuration

  !if(m.eq.12) print *, a
  !origin=a

  CALL  DGESVD( 'A', 'A', M, N, A, m, S, B, m, VT, n, &
       WORK, LWORK, INFO )
  
  !if(info.ne.0) return

  
  maxS=maxval(S)
  !print *, VT(1,:)
  !print *, VT(2,:)
  do i = 1,M
     ii=i
     if(S(i)<eps*maxS) exit
     S(i) = 1.d0/S(i)    
  enddo

  ii=ii-1

 
  
  !tmpA=0.d0
  !do i = 1,ii
  !   do j = 1,M
  !      tmpA(i,j) = B(j,i) / S(i)
  !   enddo
  !enddo

  ia=0.d0
  
  do i = 1,n
     do j = 1,m
        do k = 1,ii
           ia(i,j) = ia(i,j) + VT(k,i)*B(j,k)*S(k)         
        enddo
     enddo
  enddo
  !print *, "tmpA",tmpA
      
  !ia=0.d0
  !ia=matmul(transpose(VT(1:ii,1:N)),tmpA(1:ii,1:M))
  !ia=matmul(tmpA(1:N,1:ii),transpose(B(1:ii,1:M))
  !print *, ia
  

  return
end subroutine svdinverse


subroutine inversebug(m,a,n,PINV)

   implicit none

   
   
   external ZLANGE
   double precision :: ZLANGE
   
   integer :: i, ii,j, M, N,NMAX,K, L, LWORK, INFO


   !K = MIN(M,N) but here we assume K=M and L=N
   !L = MAX(M,N)

   ! LWORK = MAX(1,2*K+L) so this should be 2*M+N

   double precision, dimension(M,N) :: A, A1, A2, SIGMA
   double precision, dimension(N,M) :: PINV
   double precision, dimension(M,M) :: U
   double precision, dimension(N,N) :: VT
   double precision, dimension(N,N) :: BUFF
   double precision, dimension(2*M+N) :: WORK,WORK2
   double precision, dimension(5*M) :: RWORK
   double precision, dimension(M) :: S
   !integer, dimension(4) :: ISEED
   double precision :: eps

   double precision :: normA, normAPA, normPAP

   return

   K=M
   L=N

   

   eps = 1.d-8

   !  Fill A1 with random values and copy into A2
   !call ZLARNV( 1, ISEED, M*N, A1 )
   
   ! copy A into A1 and A2

   do i=1,M
      do j=1,N
         A1(i,j) = A(i,j)
         A2(i,j) = A(i,j)
      enddo
   enddo

   !print *, "A"



   !  Compute the SVD of A1
   call DGESVD( 'S', 'S', M, N, A1, M, S, U, M, VT, K, WORK, LWORK, &
        RWORK, INFO) 
   !print *, "S"
   !print *, S(1:M)
   

   do i = 1,M
      ii=i
      if(S(i)<eps) exit
   enddo

   A = 0.d0 

   do i = 1,ii
      do j = 1,M
         A(i,j) = U(j,i) /S(i)
      enddo
   enddo
   !print *,"1/S"


   ! A(1:ii,1:M) = matmul(A(1:ii,1:ii),transpose(U(1:M,1:ii)))
   !print *, " before"
   PINV=matmul(transpose(VT(1:ii,1:M)),A(1:ii,1:M))
   !print *, "PINV"

 end subroutine inversebug


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
     if(abs(ia(i,i)) .lt. 1.D-10) then
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

