!-----------------------------------------------------------------------
      subroutine mygauss(ear, ne, photar)
      use rtconstants, only: wp
      implicit none
      integer ne,i
      real(wp) ear(0:ne),photar(ne),E
      do i = 1,ne
        E = 0.5_wp * ( ear(i) + ear(i-1) )
        photar(i) = exp( -(E-6.4_wp)**2.0_wp/(2.0_wp*(0.02_wp)**2.0_wp) )
        photar(i) = photar(i) * ( ear(i) - ear(i-1) )
      end do
      return
      end
!-----------------------------------------------------------------------
