!-----------------------------------------------------------------------
      subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
      use rtconstants, only: wp
      implicit none

      integer         , intent(in)    :: nro, nphi
      real(wp), intent(in)    :: rnmin, rnmax, mueff
      real(wp), intent(out)   :: domega(nro)
      real(wp), intent(inout) :: rn(nro)
      real(wp), parameter     :: pi = acos(-1.0_wp)
      real(wp) rar(0:nro), dlogr, logr
      integer i
      rar(0) = rnmin
      dlogr  = log10( rnmax/rnmin ) / real(nro, wp)
      do i = 1,NRO
        logr = log10(rnmin) + real(i, wp) * dlogr
        rar(i)    = 10.0_wp**logr
        rn(i)     = 0.5_wp * ( rar(i) + rar(i-1) )
        domega(i) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.0_wp * pi / real(nphi, wp)
      end do
      return
      end subroutine getrgrid
!-----------------------------------------------------------------------
