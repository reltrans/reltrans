!-----------------------------------------------------------------------
      subroutine mygauss(ear, ne, param, ifl, photar)
      implicit none
      integer ne,ifl,i
      real ear(0:ne),param(*),photar(ne),E
      do i = 1,ne
        E = 0.5 * ( ear(i) + ear(i-1) )
        photar(i) = exp( -(E-6.4)**2.0/(2.*(0.02)**2.) )
        photar(i) = photar(i) * ( ear(i) - ear(i-1) )
      end do
      return
      end
!-----------------------------------------------------------------------
