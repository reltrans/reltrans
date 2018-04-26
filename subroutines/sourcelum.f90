!-----------------------------------------------------------------------
      subroutine sourcelum(nex,earx,contx,mass)
! Calculates implies source luminosity in units of the Eddington limit       
      integer nex,i
      real earx(0:nex),contx(nex),integral,E,dE,lambda,mass
      integral = 0.0
      do i = 1,nex
        E  = 0.5 * ( earx(i) + earx(i-1) )
        dE =         earx(i) - earx(i-1)
        integral = integral + E * contx(i) * dE
      end do
      lambda = 9.123e-3 * (6.0/Mass) * integral
      write(*,*)"Ls/Ledd=norm*(D/6kpc)*",lambda
      return
      end subroutine sourcelum
!-----------------------------------------------------------------------
