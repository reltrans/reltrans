!-----------------------------------------------------------------------
      subroutine ch_ind_val(min,max,dim,ind,l)
! given the index (ind) of a linear array in which min, max values and the
! dimension are known, the subroutine gives back the value of the array at
! that particular index. (this is useful only if you don't have the array)
! note: if the array starts from 0 this function isn't valid       
      implicit none
      double precision :: min,max,l
      integer :: dim,ind

      if(dim.eq.1)then
         l = min
      else if (dim.lt.1)then
       write(*,*) "error: in rch_ind_val the array dim is less than one" 
      else if (dim.gt.1) then
         l = float(ind-1)*(max-min)/float(dim-1) + min
      endif
      
      end
!-----------------------------------------------------------------------
