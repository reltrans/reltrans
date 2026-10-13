!-----------------------------------------------------------------------
subroutine sourcelum(nex,earx,contx,mass,gso,gamma)
!> Calculates implied source luminosity in units of the Eddington limit       
    use rtconstants, only: wp
    integer nex,i
    real(wp) earx(0:nex),contx(nex),integral,E,lambda,mass
    real(wp) gso,Gamma,F
! contx(i) is photar(i), so no need to multiply by dE
    integral = 0.0_wp
    F        = 0.0_wp
    do i = 1,nex
        E  = 0.5_wp * ( earx(i) + earx(i-1) )
        integral = integral + E * contx(i)
        if( E .gt. 13.6e-3_wp .and. E .gt. 13.6_wp ) F = F + E * contx(i)
    end do
    ! write(*,*) "for lambda calculation", gso, gamma, integral, mass
    lambda = 1.5217e-3_wp * gso**(gamma-2.0_wp) * integral / mass
    !write(*,*)"\mathcal{F} = norm * ",integral*1.6e-9,"erg/cm^2/s"
    write(*,*)"Ls/Ledd = norm * (D/kpc)**2 *",lambda
    write(*,*)"Lacc/Ledd = norm * (D/kpc)**2 *",2.0_wp*lambda
! Ls = A * 4*pi*D**2 * gso**(Gamma-2) * I; units erg/s
! A        = reltrans norm; units = cm^{-2}
! D        = Distance; units = cm
! I        = the above integral; units = erg/s
! D        = ( D / kpc ) * 3.086e21 cm
! I        = ( integral / keV/s ) 1.6e-9 erg/s
! Ledd     = 1.26e38 (M/Msun) erg/s
! Ls/Ledd  = A * 4*pi * (D/kpc)**2 * 9e42 * gso**(Gamma-2) * integral * 1.6e-9 / [ 1.26e38 (M/Msun) ]
!          = 1.5217e-3 * A * (D/kpc)**2 * gso**(Gamma-2) * integral / mass
    return
end subroutine sourcelum
!-----------------------------------------------------------------------
