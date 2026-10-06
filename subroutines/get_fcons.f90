!-----------------------------------------------------------------------
function get_fcons(h,spin,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,photarx,dlogE)
! Fx(r) = fcons * eps_bol(r)
! fcons is in units of erg cm^{-2} s^{-1}
  use rtconstants, only: wp
  implicit none
  real(wp)             :: get_fcons
  integer, intent(in)          :: nex
  real(wp), intent(in) :: h,spin,zcos,Gamma
  real(wp), intent(in)             :: Dkpc,Mass,Anorm,earx(0:nex)
  real(wp), intent(in)             :: photarx(nex),dlogE
  real(wp), parameter              :: pi = acos(-1.0_wp)
  real(wp) :: gso,dgsofac
  real(wp)         :: integral,Eintegrate 
  gso       = dgsofac(spin,h) / ( 1.0_wp + zcos )
  integral  = Anorm * Eintegrate(0.1_wp,1e3_wp,nex,earx,photarx,dlogE)
  get_fcons = 4.0_wp * pi * (Dkpc/Mass)**2 * gso**(Gamma-2.0_wp)
  get_fcons = get_fcons * integral * 6.99367e23_wp
  !above constant is (1kpc/Rgsun)^2 * 1 keV in erg
  !Calculated to high precision keeping all decimal places
  return
end function get_fcons
!-----------------------------------------------------------------------
