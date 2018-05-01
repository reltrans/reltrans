!-----------------------------------------------------------------------
      subroutine chclose(a,amin,amax,adim,ilo,ihi)
!     purpouse: given the lowest and the highest values of a lin spaced array,
!      it returns the closest indices of the array 
!     note: if you want to use an array like v(0:dim) you must pass to the function
!      dim+1 as a adim. then subtract -1 to the result indices
      implicit none
      double precision :: a,amin,amax,el
      integer adim,ihi,ilo
      if (a.lt.amin.or.a.gt.amax) then
         write(*,*) "value", a
         write(*,*) "min", amin
         write(*,*) "max", amax
         
         write(*,*) 'error: in chclose the value is out of the array'
         goto 666
      endif

      if (adim.eq.1) then
         write(*,*)'warning: in chclose the array dimension is 1'
         ihi = 1
         ilo = 1
      endif

      if (adim.lt.1) then
       write(*,*) 'warning: in chclose the array dimension is less 1'
       goto 666
      endif
      el = (a-amin)*float(adim-1)/(amax-amin) + 1.0
      ihi = min( ceiling( el ) , adim )
      ilo = max( floor( el )   , 1    )      
      return
 666  stop
      end
!-----------------------------------------------------------------------
