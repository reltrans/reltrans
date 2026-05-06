!-----------------------------------------------------------------------
      subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
      implicit none

      integer         , intent(in)    :: nro, nphi
      double precision, intent(in)    :: rnmin, rnmax, mueff
      double precision, intent(out)   :: domega(nro)
      double precision, intent(inout) :: rn(nro)
      double precision, parameter     :: pi = acos(-1.0)
      double precision rar(0:nro), dlogr, logr
      integer i
      rar(0) = rnmin
      dlogr  = log10( rnmax/rnmin ) / dble(nro)
      do i = 1,NRO
        logr = log10(rnmin) + dble(i) * dlogr
        rar(i)    = 10.d0**logr
        rn(i)     = 0.5 * ( rar(i) + rar(i-1) )
        domega(i) = rn(i) * ( rar(i) - rar(i-1) ) * mueff * 2.d0 * pi / dble(nphi)
      end do
      return
      end subroutine getrgrid
!-----------------------------------------------------------------------
