!-----------------------------------------------------------------------
function mylogne(r,rin)
! Calculates log10(ne). Don't let r = rin
  use rtconstants, only: wp
  implicit none
  real(wp) mylogne,r,rin
  mylogne = 1.5_wp * log10(r) - 2.0_wp * log10( 1.0_wp - sqrt(rin/r) )
  return
end function mylogne 
!-----------------------------------------------------------------------
