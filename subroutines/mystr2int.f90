!-----------------------------------------------------------------------
      subroutine mystr2int(str,int,stat)
      implicit none
      character (len=*) str
      integer int,stat
      read(str,*,iostat=stat)  int
      return
      end subroutine
!-----------------------------------------------------------------------
