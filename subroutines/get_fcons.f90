!-----------------------------------------------------------------------
function get_fcons(h,spin,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,photarx,dlogE)
! Fx(r) = fcons * eps_bol(r)
! fcons is in units of erg cm^{-2} s^{-1}
  implicit none
  double precision             :: get_fcons
  integer, intent(in)          :: nex
  double precision, intent(in) :: h,spin,zcos,Gamma
  real, intent(in)             :: Dkpc,Mass,Anorm,earx(0:nex)
  real, intent(in)             :: photarx(nex),dlogE
  real, parameter              :: pi = acos(-1.0)
  double precision :: gso,dgsofac
  real             :: integral,Eintegrate 
  gso       = dgsofac(spin,h) / ( 1.0 + zcos )
  integral  = Anorm * Eintegrate(0.1,1e3,nex,earx,photarx,dlogE)
  get_fcons = 4.0 * pi * (Dkpc/Mass)**2 * gso**(Gamma-2.0)
  get_fcons = get_fcons * integral * 6.99367e+23
  !above constant is (1kpc/Rgsun)^2 * 1 keV in erg
  !Calculated to high precision keeping all decimal places
  return
end function get_fcons
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function get_lacc(h,spin,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,photarx,dlogE)
! Returns accretion luminosity as a fraction of Eddington
  implicit none
  double precision             :: get_lacc
  integer, intent(in)          :: nex
  double precision, intent(in) :: h,spin,zcos,Gamma
  real, intent(in)             :: Dkpc,Mass,Anorm,earx(0:nex)
  real, intent(in)             :: photarx(nex),dlogE
  real, parameter              :: pi = acos(-1.0)
  double precision :: gso,dgsofac
  real             :: integral,Eintegrate 
  gso       = dgsofac(spin,h) / ( 1.0 + zcos )
  integral  = Anorm * Eintegrate(earx(0),earx(nex),nex,earx,photarx,dlogE)
  get_lacc  = 4.0 * pi * Dkpc**2 * gso**(Gamma-2.0)
  get_lacc  = get_lacc * integral * 1.21097e-4 / Mass
  get_lacc  = get_lacc * 2.0
  return
end function get_lacc
!-----------------------------------------------------------------------
