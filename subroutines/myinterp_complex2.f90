!-----------------------------------------------------------------------
      subroutine myinterp_complex2(l,llo,lhi,line_dim,line1_lo,line1_hi &
           ,line2_lo,line2_hi,line1,line2)
! this subroutive calculates the interpolation between line_lo(llo) and       
! line_hi(lhi). the output is line(l). l is the varible of the interpolation,
! llo and lhi are respectively the closest low and hight values to l. 
! note: if you used chclose to extract the index of the closest values you       
!      also need ch_ind_val to find the corresponding array values.
      implicit none
      integer line_dim,i
      double precision ::l,llo,lhi
      complex :: line1_lo(line_dim),line1_hi(line_dim),line1(line_dim),grad1,cons1 
      complex :: line2_lo(line_dim),line2_hi(line_dim),line2(line_dim),grad2,cons2 
      double precision :: da
      if ( lhi .lt. llo ) then
         write(*,*)
         write(*,'(a)', advance="no") 'error in myinterp: low index is' 
         write(*,*) ' bigger than high index'
         goto 666
      else if( lhi .eq. llo )then
         do i=1,line_dim
            line1(i) = line1_lo(i)
            line2(i) = line2_lo(i)
         enddo
        return
      end if
      da = lhi - llo
      do i = 1,line_dim
        grad1 = ( line1_hi(i) - line1_lo(i) ) / real(da)
        cons1 = line1_lo(i) - grad1*real(llo)
        line1(i) = grad1*real(l) + cons1
        grad2 = ( line2_hi(i) - line2_lo(i) ) / real(da)
        cons2 = line2_lo(i) - grad2*real(llo)
        line2(i) = grad2*real(l) + cons2
      end do
      return
 666  stop
      end
!-----------------------------------------------------------------------
