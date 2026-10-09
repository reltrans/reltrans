!-----------------------------------------------------------------------
function get_lacc(h,spin,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,photarx,dlogE)
! Returns accretion luminosity as a fraction of Eddington
  use rtconstants, only: wp
  implicit none
  real(wp)             :: get_lacc
  integer, intent(in)          :: nex
  real(wp), intent(in) :: h,spin,zcos,Gamma
  real(wp), intent(in)             :: Dkpc,Mass,Anorm,earx(0:nex)
  real(wp), intent(in)             :: photarx(nex),dlogE
  real(wp), parameter              :: pi = acos(-1.0_wp)
  real(wp) :: gso,dgsofac
  real(wp)         :: integral,Eintegrate 
  gso       = dgsofac(spin,h) / ( 1.0_wp + zcos )
  integral  = Anorm * Eintegrate(earx(0),earx(nex),nex,earx,photarx,dlogE)
  get_lacc  = 4.0_wp * pi * Dkpc**2 * gso**(Gamma-2.0_wp)
  get_lacc  = get_lacc * integral * 1.21097e-4_wp / Mass
  get_lacc  = get_lacc * 2.0_wp
  ! write(*,*) 'in get_lacc', get_lacc, integral, Mass, Anorm, gso
  return
end function get_lacc
!-----------------------------------------------------------------------
