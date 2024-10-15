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
  ! write(*,*) 'in get_lacc', get_lacc, integral, Mass, Anorm, gso
  return
end function get_lacc
!-----------------------------------------------------------------------
