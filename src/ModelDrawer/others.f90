
subroutine readpsvmodel(psvmodel,tmpfile)
  implicit none
  character(120) :: psvmodel, tmpfile, dummy
  open(unit=2, file=psvmodel, status='old',action='read',position='rewind')
  open(unit=1, file=tmpfile,status='unknown')
100 continue
  read(2,110) dummy
110 format(a120)
  if(dummy(1:1).eq.'#') goto 100
  if(dummy(1:3).eq.'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(2)
end subroutine readpsvmodel



subroutine calstg_for_card(r,nzone,vrmin,vrmax,rrho,vpv,vph,vsv,vsh,eta,qmu,qkappa,array)

  ! Computing the structure grid points.
  implicit none
  integer:: nzone
  real(kind(0d0)):: r,rrho(4,nzone),vpv(4,nzone),vph(4,nzone),vsv(4,nzone),vsh(4,nzone),eta(4,nzone)
  real(kind(0d0)):: qmu(nzone), qkappa(nzone),vrmin(nzone),vrmax(nzone)
  real(kind(0d0)), parameter:: rmax  = 6371.d0
  real(kind(0d0)):: rho,ecKx,ecKy,ecKz
  real(kind(0d0)):: ecL,ecN
  real(kind(0d0)):: ecA,ecC,ecF
  real(kind(0d0)):: trho,tvpv,tvph,tvsv,tvsh,teta,coef
  complex(kind(0d0)):: coef1,coef2
  integer:: izone,j
  real(kind(0d0)):: array(1:9)
  
  array = 0.d0
  do izone = 1, nzone
     if((r.gt.vrmin(izone)).and.(r.le.vrmax(izone))) then
        
        coef1 = cmplx(qmu(izone))
        coef2 = cmplx(qkappa(izone))
        trho = 0.d0
        tvpv = 0.d0
        tvph = 0.d0
        tvsv = 0.d0
        tvsh = 0.d0
        teta = 0.d0
        do j=1,4
           if ( j.eq.1 ) then
              coef = 1.d0
           else
              coef = coef * (r / rmax )
           endif
           trho  = trho  + rrho(j,izone)  * coef
           tvpv  = tvpv  + vpv(j,izone)   * coef
           tvph  = tvph  + vph(j,izone)   * coef
           tvsv  = tvsv  + vsv(j,izone)   * coef
           tvsh  = tvsh  + vsh(j,izone)   * coef
           teta  = teta  + eta(j,izone)   * coef
        enddo
        rho = trho
        ecL  = rho * tvsv * tvsv
        ecN  = rho * tvsh * tvsh
        ecA = trho * tvph * tvph
        ecC = trho * tvpv * tvpv
        ecF = teta * ( ecA - 2.d0 * ecL )
        !kappa(itmp) = ( 4.d0 * ecA + ecC  + 4.d0 * ecF - 4.d0 * ecN(itmp) ) / 9.d0
        ecKx = ecA - 4.d0 / 3.d0 * ecN
        ecKy = ecF + 2.d0 / 3.d0 * ecN
        ecKz = ( ecC + 2.d0 * ecF ) / 3.d0


        array(1) = 1.d3 * r
        array(2) = 1.d3 * rho
        array(3) = 1.d3 * tvpv
        array(4) = 1.d3 * tvsv
        array(5) = qkappa(izone)
        array(6) = qmu(izone)
        array(7) = 1.d3 * tvph
        array(8) = 1.d3 * tvsh
        array(9) = teta
     endif
  enddo
 

  
  return
end subroutine calstg_for_card

