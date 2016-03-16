



subroutine calbvecphi0( l,theta,plm,bvec,bvecdt,bvecdp)
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  integer  :: l,m,i,j
  real(kind(0d0)) :: theta,x,plm(1:3,0:3),fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  complex(kind(0d0)) :: bvecdt(1:3,-2:2),bvecdp(1:3,-2:2)
  real(kind(0d0)) :: plmdt,xl2
  real(kind(0d0)) :: rtxl2,coeff,rtxl22,sign1,sign2
  

  if((theta.eq.0.d0).or.(theta.eq.pi)) then
     sign1 = 1.d0
     sign2 = 1.d0

     if(theta.eq.pi) then
        sign1 = dble((-1)**l)
        sign2 = dble((-1)**(l+1)) 
     endif
     
     bvec = cmplx(0.d0)
     bvecdt = cmplx(0.d0)
     bvecdp = cmplx(0.d0)
     return
  endif



  x = dcos( theta )
  xl2 = dble(l) * dble(l+1)
  do m=0,min0(l,3)
     call calplm( l,m,x,plm(1:3,m))
  enddo 

  do m=0,min0(l,2)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     plmdt = dble(m) * x / sin( theta ) * plm(1,m) + plm(1,m+1)
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m) ) / dsin( theta) * coef * plm(1,m) 
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) = - coef * plmdt 
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     ! calculate derivatives
     bvecdt(1,m)  = dcmplx( 0.d0 )
     bvecdt(1,-m) = dcmplx( 0.d0 )
     bvecdt(2,m)  = dcmplx( 0.d0, dble(m) ) * ( plmdt / dsin(theta) &
          - x / ( 1 - x * x ) * plm(1,m) ) * coef 
     bvecdt(2,-m) = dcmplx( conjg( bvecdt(2,m) ))
     bvecdt(3,m) = ( x / dsin(theta) * plmdt - dble(m) * dble(m)/(1-x*x) *plm(1,m) &
          &           + xl2 * plm(1,m) ) * coef 
     bvecdt(3,-m) = dcmplx(conjg( bvecdt(3,m)) )
     bvecdp(1,m)  = dcmplx( 0.d0 )
     bvecdp(1,-m) = dcmplx( 0.d0 )
     bvecdp(2,m)  = - dble(m) * dble(m) / dsin(theta) * plm(1,m) * coef 
     bvecdp(2,-m) = dcmplx(conjg( bvecdp(2,m)) )
     bvecdp(3,m)  = - dcmplx( 0.d0, dble(m) ) * plmdt * coef 
     bvecdp(3,-m) = dcmplx(conjg( bvecdp(3,m)) )

     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
        bvecdt(2,-m) = - bvecdt(2,-m)
        bvecdt(3,-m) = - bvecdt(3,-m)
        bvecdp(2,-m) = - bvecdp(2,-m)
        bvecdp(3,-m) = - bvecdp(3,-m)
     endif
  enddo



  return
end subroutine calbvecphi0



subroutine calplm( l,m,x,plm )
  implicit none
  integer :: l,m,i
  real(kind(0d0)) :: x,plm(1:3),pmm,somx2,fact

  if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.d0)) pause 'bad arguments'
  if ( l.eq.m ) then
     pmm = 1.d0
     if ( m.gt.0 ) then
        somx2 = dsqrt( (1.d0-x)*(1.d0+x) )
        fact = 1.d0
        do i=1,m
           pmm = -pmm * fact * somx2
           fact = fact + 2.d0
        enddo
     endif
     plm(3) = 0.d0
     plm(2) = 0.d0
     plm(1) = pmm
  else
     plm(3) = plm(2)
     plm(2) = plm(1)
     if ( l.eq.m+1 ) then
        plm(1) = x * dble(2*m+1) * plm(2)
     else
        !print *, l,m,x
        plm(1) = (x*dble(2*l-1) * plm(2)-dble(l+m-1) * plm(3) )/dble(l-m)
        !print *, plm(1)
     endif
  endif


end subroutine calplm


subroutine calbveczero( l,bvec )
  
  implicit none
  real(kind(0d0)), parameter ::  pi=3.1415926535897932d0 
  
  integer  :: l,m,i
  real(kind(0d0)) :: fact,coef
  complex(kind(0d0)) :: bvec(1:3,-2:2)
  real(kind(0d0)) :: xl2


  xl2 = dble(l) * dble(l+1)

  do m=0,min0(l,1)
     fact = 1.d0
     if ( m.ne.0 ) then
        do i=l-m+1,l+m
           fact = fact * dble(i)
        enddo
     endif
     coef = dsqrt( dble(2*l+1)/(4.d0*pi) / fact )
     bvec(1,m)  = dcmplx( 0.d0 )
     bvec(1,-m) = dcmplx( 0.d0 )
     bvec(2,m)  = dcmplx( 0.d0, dble(m)) *  xl2 *coef / 2.d0
     bvec(2,-m) = dcmplx(conjg( bvec(2,m)) )
     bvec(3,m) =  -dcmplx(dble(m),0.d0) * xl2 * coef / 2.d0
     bvec(3,-m) = dcmplx(conjg( bvec(3,m) ))

     if ( mod(m,2).eq.1 ) then
        bvec(2,-m) = - bvec(2,-m)
        bvec(3,-m) = - bvec(3,-m)
     endif
  enddo



  return
end subroutine calbveczero
