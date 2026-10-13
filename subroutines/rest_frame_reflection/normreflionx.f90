!-----------------------------------------------------------------------
subroutine normreflionx(ear,ne,Gamma,Afe,logne,kTe,logxi,thetae,photar)
! !
! ! Returns reflionx model renormalised with xillverDCp
! !  
  use rtconstants, only: wp
  implicit none
  integer ne,ifl
  integer, parameter :: dim = 6, dimCp = 7
  real(wp) ear(0:ne),Gamma,Afe,logne,kTe,logxi,thetae,photar(ne)
  real(wp) kTbb, param(7), xillpar(dim), xillparDCp(dimCp), xillphotar(ne)
  real(wp) E,rintegral,xintegral,fac,lognex
  integer i,Cp,ilo,ihi
! Set integration bounds
  ilo = ceiling( log( 50.0_wp / ear(0) ) / log(ear(ne)/ear(0)) * real(ne, wp) )
  ilo = max( ilo , 1 )
  ihi = floor( log( 100.0_wp / ear(0) ) / log(ear(ne)/ear(0)) * real(ne, wp) )
  ihi = min( ihi , ne )
! Hardwired seed photon temperature
  kTbb = 0.05_wp
! Set reflionx parameters
  param(1) = 10.0_wp**logne
  param(2) = 10.0_wp**logxi
  param(3) = Afe
  param(4) = kTe
  param(5) = Gamma
  param(6) = kTbb
  param(7) = 0.0_wp !Redshift
  ifl = 0
! Call model
  call get_reflionx(ear, ne, param, ifl, photar)
! Call xillver equivalent
  xillpar(1) = Gamma
  xillpar(2) = Afe    !
  xillpar(3) = 15.0_wp !logne or Ecut or kTe
  xillpar(4) = logxi  !ionization par
  xillpar(5) = thetae !emission angle
  xillpar(6) = 0.0_wp !redshift
  xillparDCp(1) = Gamma  !photon index
  xillparDCp(2) = Afe    !Afe
  xillparDCp(3) = logxi  !ionization par
  xillparDCp(4) = kTe    !kTe
  lognex = logne
  lognex = min(logne,20.0_wp)
  lognex = max(logne,15.0_wp)
  xillparDCp(5) = lognex   !logne
  xillparDCp(6) = thetae !emission angle
  xillparDCp(7) = 0.0_wp !redshift
  Cp = 2
  call get_xillver(ear, ne, dim, dimCp, xillpar, xillparDCp, Cp, xillphotar)
! Integrate both spectra
  rintegral = 0.0_wp
  xintegral = 0.0_wp
  do i = ilo,ihi
     E = 0.5_wp * ( ear(i) + ear(i-1) )
     rintegral = rintegral + E * photar(i)
     xintegral = xintegral + E * xillphotar(i)
  end do
  fac = xintegral / rintegral 
! renormalise reflionx spectrum
  photar = photar * fac
  return
end subroutine normreflionx
!-----------------------------------------------------------------------

