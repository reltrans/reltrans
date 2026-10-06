!-----------------------------------------------------------------------
function interper(rlp,ylp,ndelta,re,kk)
! Interpolates the array ylp between ylp(kk-1) and ylp(kk)
! ylp is a function of rlp and rlp(kk-1) .le. re .le. rlp(kk)
  use rtconstants, only: wp
  implicit none
  integer ndelta,kk
  real(wp) interper,rlp(ndelta),ylp(ndelta),re
  interper = (ylp(kk)-ylp(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
  interper = interper + ylp(kk)
  return
end function interper
!-----------------------------------------------------------------------
