!-----------------------------------------------------------------------
subroutine myreflionx(ear, ne, param, ifl, photar)
  implicit none
  integer, intent(in)  :: ne, ifl
  real,    intent(in)  :: ear(0:ne), param(7)
  real,    intent(out) :: photar(ne)
  real                 :: photer(ne)
  character (len=500)  :: filenm,strenv
  character (len=200)  :: envnm
  logical              :: needfile
  data needfile/.true./
  save needfile
  
! Get the reflionx table  
  if( needfile )then
     envnm  = 'REFLIONX_FILE'
     filenm = strenv(envnm)
     if( trim(filenm) .eq. 'none' )then
        write(*,*)"Enter reflionx file (with full path)"
        read(*,'(a)')filenm
     end if
     needfile = .false.
  end if
  
! Interpolate a spectrum from it
  call xsatbl(ear, ne, param, filenm, ifl, photar, photer) 

  return
end subroutine myreflionx
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine normreflionx(ear,ne,Gamma,Afe,logne,kTe,logxi,thetae,photar)
!
! Returns reflionx model renormalised with xillverDCp
!  
  implicit none
  integer ne,ifl,j,jmax
  parameter (jmax=20)
  real ear(0:ne),Gamma,Afe,logne,kTe,logxi,thetae,photar(ne)
  real kTbb,param(7),xillpar(7),xillparDCp(8),xillphotar(ne)
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
! Call model
  call myreflionx(ear, ne, param, ifl, photar)
! Call xillver equivalent
  xillpar(1) = Gamma
  xillpar(2) = Afe    !
  xillpar(3) = 15.0   !logne or Ecut or kTe
  xillpar(4) = logxi  !ionization par
  xillpar(5) = 0.0    !redshift
  xillpar(6) = thetae !emission angle
  xillpar(7) = -1.0   !refl_frac
  xillparDCp(1) = Gamma  !photon index
  xillparDCp(2) = Afe    !Afe
  xillparDCp(3) = kTe    !Ecut or kTe
  lognex = logne
  lognex = min(logne,20.0)
  lognex = max(logne,15.0)
  xillparDCp(4) = lognex   !logne
  xillparDCp(5) = logxi  !ionization par
  xillparDCp(6) = 0.0    !redshift
  xillparDCp(7) = thetae !emission angle
  xillparDCp(8) = -1.0   !refl_frac  
  Cp = 2
  call myxill(ear, ne, xillpar, xillparDCp, ifl, Cp, xillphotar)
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

