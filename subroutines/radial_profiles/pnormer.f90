!-----------------------------------------------------------------------
function pnormer(b1,b2,boost)
  use rtconstants, only: wp
  implicit none
  real(wp) pnormer,b1,b2,boost,pi
  integer i,imax
  parameter (imax=1000)
  real(wp) integral,mu,pfunc_raw
  pi = acos(-1.0_wp)
  integral = 0.0_wp
  do i = 1,imax
     mu       = -1.0_wp + 2.0_wp*real(i-1, wp)/real(imax, wp)
     integral = integral + pfunc_raw(mu,b1,b2,boost)
  end do
  integral = integral * 2.0_wp / real(imax, wp)
  pnormer  = 1.0_wp / ( 2.0_wp*pi*integral )
  return
end function pnormer
!-----------------------------------------------------------------------
