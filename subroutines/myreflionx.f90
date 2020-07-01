!-----------------------------------------------------------------------
      subroutine myreflionx(ear, ne, param, ifl, photar)
      implicit none
      integer, intent(in)  :: ne, ifl
      real,    intent(in)  :: ear(0:ne), param(7)
      real,    intent(out) :: photar(ne)
      real                 :: photer(ne)
      character (len=200)  :: filenm


      filenm = '/Users/gullo/Software/reflionx/reflionx_HD_nthcomp_v2.fits'
      call xsatbl(ear, ne, param, filenm, ifl, photar, photer) 

      return
    end subroutine myreflionx
!-----------------------------------------------------------------------
