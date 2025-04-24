!-----------------------------------------------------------------------
function pnormer(b1,b2,boost)
  implicit none
  double precision pnormer,b1,b2,boost,pi
  integer i,imax
  parameter (imax=1000)
  double precision integral,mu,pfunc_raw
  pi = acos(-1.d0)
  integral = 0.0
  do i = 1,imax
     mu       = -1.0 + 2.0*dble(i-1)/dble(imax)
     integral = integral + pfunc_raw(mu,b1,b2,boost)
  end do
  integral = integral * 2.0 / dble(imax)
  pnormer  = 1.0 / ( 2.0*pi*integral )
  return
end function pnormer
!-----------------------------------------------------------------------
