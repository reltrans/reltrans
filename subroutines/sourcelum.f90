!-----------------------------------------------------------------------
      subroutine sourcelum(nex,earx,contx,mass,gso,gamma)
! Calculates implies source luminosity in units of the Eddington limit       
      integer nex,i
      real earx(0:nex),contx(nex),integral,E,lambda,mass
      real gso,Gamma,F
! contx(i) is photar(i), so no need to multiply by dE
      integral = 0.0
      F        = 0.0
      do i = 1,nex
        E  = 0.5 * ( earx(i) + earx(i-1) )
        integral = integral + E * contx(i)
        if( E .gt. 13.6e-3 .and. E .gt. 13.6 ) F = F + E * contx(i)
      end do
      lambda = 1.436e-3 * gso**(gamma-2.0) * integral / mass
      write(*,*)"\mathcal{F} = norm * ",F*1.6e-9,"erg/cm^2/s"
      write(*,*)"Ls/Ledd = norm * (D/kpc)**2 *",lambda
! Ls = A * 4*pi*D**2 * gso**(Gamma-2) * I; units erg/s
! A        = reltrans norm; units = cm^{-2}
! D        = Distance; units = cm
! I        = the above integral; units = erg/s
! D        = ( D / kpc ) * 3e21 cm
! I        = ( integral / keV/s ) 1.6e-9 erg/s
! Ledd     = 1.26e38 (M/Msun) erg/s
! Ls/Ledd  = A * 4*pi * (D/kpc)**2 * 9e42 * gso**(Gamma-2) * integral * 1.6e-9 / [ 1.26e38 (M/Msun) ]
!          = 1.436e-3 * A * (D/kpc)**2 * gso**(Gamma-2) * integral / mass
      return
      end subroutine sourcelum
!-----------------------------------------------------------------------
