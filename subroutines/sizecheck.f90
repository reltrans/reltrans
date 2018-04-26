!-----------------------------------------------------------------------
      subroutine sizecheck(me,mex)
      implicit none
      integer me,mex
      if( me .gt. mex )then
        me = mex
        write(*,*)"Warning! Too many zones, set to maximum allowed value"
     end if
     return
     end subroutine sizecheck
!-----------------------------------------------------------------------
