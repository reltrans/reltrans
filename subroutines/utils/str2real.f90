!-----------------------------------------------------------------------
subroutine str2real(str, real_num, stat)
  implicit none
  character (len=*) str
  integer stat
  real real_num
  read(str,*,iostat=stat)  real_num
  return
end subroutine str2real
!-----------------------------------------------------------------------
