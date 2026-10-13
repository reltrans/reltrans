!-----------------------------------------------------------------------
function get_index(rlp,ndelta,re,rmin,npts)
  use rtconstants, only: wp
  implicit none
  integer get_index,ndelta,npts,kk
  real(wp) rlp(ndelta),re,rmin
  kk = 2
  do while( ( rlp(kk) .le. re .or. rlp(kk-1) .lt. rmin ) .and. kk .lt. npts )
     kk = kk + 1
  end do
  get_index = kk
  return
end function get_index
!-----------------------------------------------------------------------
