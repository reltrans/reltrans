!-----------------------------------------------------------------------
subroutine normreflionx(ear,ne,Gamma,Afe,logne,kTe,logxi,thetae,photar)
!
! Returns reflionx model renormalised with xillverDCp
!  
  implicit none
  integer ne,ifl,j,jmax
  parameter (jmax=20)
  real ear(0:ne),Gamma,Afe,logne,kTe,logxi,thetae,photar(ne)
  real kTbb, param(7), xillpar(6), xillparDCp(7), xillphotar(ne)
  real E,dE,rintegral,xintegral,fac,lognej,lognex
  integer i,Cp,ilo,ihi
! Set integration bounds
  ilo = ceiling( log( 50.0 / ear(0) ) / log(ear(ne)/ear(0)) * real(ne) )
  ilo = max( ilo , 1 )
  ihi = floor( log( 100.0 / ear(0) ) / log(ear(ne)/ear(0)) * real(ne) )
  ihi = min( ihi , ne )
! Hardwired seed photon temperature
  kTbb = 0.05
! Set reflionx parameters
  param(1) = 10.0**logne
  param(2) = 10.0**logxi
  param(3) = Afe
  param(4) = kTe
  param(5) = Gamma
  param(6) = kTbb
  param(7) = 0.0    !Redshift
  ifl = 0
! Call model
  call get_reflionx(ear, ne, param, ifl, photar)
! Call xillver equivalent
  xillpar(1) = Gamma
  xillpar(2) = Afe    !
  xillpar(3) = 15.0   !logne or Ecut or kTe
  xillpar(4) = logxi  !ionization par
  xillpar(5) = thetae !emission angle
  xillpar(6) = 0.0    !redshift
  xillparDCp(1) = Gamma  !photon index
  xillparDCp(2) = Afe    !Afe
  xillparDCp(3) = logxi  !ionization par
  xillparDCp(4) = kTe    !Ecut or kTe
  lognex = logne
  lognex = min(logne,20.0)
  lognex = max(logne,15.0)
  xillparDCp(5) = lognex   !logne
  xillparDCp(6) = thetae !emission angle
  xillparDCp(7) = 0.0    !redshift
  Cp = 2
  call get_xillver(ear, ne, xillpar, xillparDCp, Cp, xillphotar)
! Integrate both spectra
  rintegral = 0.0
  xintegral = 0.0
  do i = ilo,ihi
     E = 0.5 * ( ear(i) + ear(i-1) )
     rintegral = rintegral + E * photar(i)
     xintegral = xintegral + E * xillphotar(i)
  end do
  fac = xintegral / rintegral 
! renormalise reflionx spectrum
  photar = photar * fac
  return
end subroutine normreflionx
!-----------------------------------------------------------------------

