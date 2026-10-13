!-----------------------------------------------------------------------
subroutine str2real(str, real_num, stat)
  use rtconstants, only: wp
  implicit none
  character (len=*) str
  integer stat
  real(wp) real_num
  read(str,*,iostat=stat)  real_num
  return
end subroutine str2real
!-----------------------------------------------------------------------
