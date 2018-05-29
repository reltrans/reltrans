!-----------------------------------------------------------------------
      subroutine Rmyinterp2_re_im(l,llo,lhi,line_dim,line1_lo_re,line1_lo_im,line1_hi_re &
           ,line1_hi_im,line2_lo_re,line2_lo_im,line2_hi_re,line2_hi_im,line1_re &
           ,line1_im,line2_re,line2_im)
! this subroutive calculates two interpolations between - line1_lo and       
! line1_hi - and  - line2_lo and line2_hi -. The outputs are line1 and line2.
! This happens for real and imaginary part. l is the varible of the interpolation.
! llo and lhi are respectively the closest low and hight values to l. 
! note: if you used chclose to extract the index of the closest values you       
!      also need ch_ind_val to find the corresponding array values.
      implicit none
      integer line_dim,i
      double precision :: l,llo,lhi,da
      real :: line1_lo_re(line_dim),line1_lo_im(line_dim),line1_hi_re(line_dim),line1_hi_im(line_dim)
      real :: grad1_re,grad1_im,cons1_re,cons1_im,grad2_re,grad2_im,cons2_re,cons2_im 
      real :: line2_lo_re(line_dim),line2_lo_im(line_dim),line2_hi_re(line_dim),line2_hi_im(line_dim) 
      real :: line1_re(line_dim),line1_im(line_dim),line2_re(line_dim),line2_im(line_dim)
      if ( lhi .lt. llo ) then
         write(*,*)
         write(*,'(a)', advance="no") 'error in myinterp: low index is' 
         write(*,*) ' bigger than high index'
         goto 666
      else if( lhi .eq. llo )then
         do i=1,line_dim
            line1_re(i) = line1_lo_re(i)
            line1_im(i) = line1_lo_im(i)
            line2_re(i) = line2_lo_re(i)
            line2_im(i) = line2_lo_im(i)
         enddo
        return
      end if
      da = lhi - llo
      do i = 1,line_dim
         grad1_re = ( line1_hi_re(i) - line1_lo_re(i) ) / real(da)
         cons1_re = line1_lo_re(i) - grad1_re*real(llo)
         line1_re(i) = grad1_re*real(l) + cons1_re

         grad1_im = ( line1_hi_im(i) - line1_lo_im(i) ) / real(da)
         cons1_im = line1_lo_im(i) - grad1_im*real(llo)
         line1_im(i) = grad1_im*real(l) + cons1_im


         grad2_re = ( line2_hi_re(i) - line2_lo_re(i) ) / real(da)
         cons2_re = line2_lo_re(i) - grad2_re*real(llo)
         line2_re(i) = grad2_re*real(l) + cons2_re

         grad2_im = ( line2_hi_im(i) - line2_lo_im(i) ) / real(da)
         cons2_im = line2_lo_im(i) - grad2_im*real(llo)
         line2_im(i) = grad2_im*real(l) + cons2_im

      end do
      return
 666  stop
      end
!-----------------------------------------------------------------------
