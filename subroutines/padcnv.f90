!-----------------------------------------------------------------------
      subroutine padcnv(dyn,ne,diskline,photard,conv)
! Convolves the array diskline(1:ne) with the array photard(1:ne) using
! FFTs, the result is recorded in conv(1:ne).
! This code uses extensive zero padding, and also applies a gate to
! get rid of high energy noise.
! Parameter dyn sets the dynamic range allowed in the output array.
! Anything smaller than dyn * ( the maximum value of conv ) will
! be set to zero. The value dyn = 1e-7 works very well.
      implicit none
      integer ne,i
      real diskline(ne),photard(ne),conv(ne),padconv(4*ne)
      real padline(4*ne),padphot(4*ne),photmax,dyn

! Fill padded arrays
      padline = 0.0
      padphot = 0.0
      do i = 1,ne
        padline(i+2*ne) = diskline(i)
        padphot(i+2*ne) = photard(i)
      end do

! Call the convolution code
      call FTcnv(4*ne,padline,padphot,padconv)

! Populate output array
      photmax = 0.0
      do i = 1,ne
        conv(i) = padconv(i+5*ne/2)
        photmax = max( photmax , conv(i) )
      end do

! Clean any residual edge effects
      do i = 1,ne
        if( abs(conv(i)) .lt. abs(dyn*photmax) ) conv(i) = 0.0
      end do
      
      return
      end subroutine padcnv
!-----------------------------------------------------------------------
