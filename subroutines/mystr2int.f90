!-----------------------------------------------------------------------
      subroutine mystr2int(str,int,stat)
      implicit none
      character (len=*) str
      integer int,stat
      read(str,*,iostat=stat)  int
      return
      end subroutine
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      subroutine mystr2real(str, real_num, stat)
      implicit none
      character (len=*) str
      integer stat
      real real_num
      read(str,*,iostat=stat)  real_num
      return
      end subroutine

!-----------------------------------------------------------------------
