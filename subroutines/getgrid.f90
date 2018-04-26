!-----------------------------------------------------------------------
      subroutine getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
! Calculates an r-grid that will be used to define impact parameters
      implicit none
      integer nro,nphi,i
      double precision rnmin,rnmax,mueff,rn(nro),domega(nro)
      double precision rar(0:nro),dlogr,logr,pi
      pi     = acos(-1.d0)
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
