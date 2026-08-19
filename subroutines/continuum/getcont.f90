!-----------------------------------------------------------------------
    subroutine getcont(Cp, earx, nex, Gamma, Cutoff_s, Cutoff_obs,             &
                       logxi, logne, zcos, contx)
!> Calculates continuum spectrum calling nthComp with the correct normalisation
!> based on the xillver spectrum
!>   Arg:
!> earx: energy grid
!> nex: number of grid points
!> Gamma: continuum spectrum inclination
!> Cutoff_s: high energy cut-off or electron temperature (source frame)
!> Cutoff_obs: high energy cut-off or electron temperature (observer frame)
!> logxi: ionisation parameter
!> logne: density
!> zcos: host galaxy redshift
!> (output) contx: continuum spectrum
!> Derivation of renormalisation constant:
!> First convert to incident flux in units of [keV/cm^2/s]
!> contx = contx * 10**(logne + logxi) / (4.0 * pi) / ergsev
!> Then divide out ad hoc factor of 1e20 included in the xillver tables to ease comparison with reflionx
!> contx = contx / 1e20
!> Then divide by the integral in the rest frame from 0.1-1e3 keV
!> contx = contx / Icomp
!> The divide by xi and by ne/1e15
!> contx = contx / (10**(logxi + logne - 15))
      use gr_continuum, only: gso
      implicit none
      integer, intent(in)           :: nex, Cp
      real   , intent(in)           :: earx(0:nex), Cutoff_s, Cutoff_obs, logxi, logne
      real   , intent(out)          :: contx(nex)
      double precision , intent(in) :: Gamma, zcos

      real   , parameter  :: pi = acos(-1.0),ergsev  = 1.602197e-9 ! Convert keV to ergs
      integer :: i, ifl, j
      real    :: nth_par(5), photer(nex), E, Icomp
      real    :: contx_base(nex), renorm, dlogE

      Icomp = 0.0

      if (Cp .eq. 2) then
!So far this works only with kTe, so only with nthComp continuum model

         !First call with no redshift
         nth_par(1) = real(Gamma)
         nth_par(2) = Cutoff_s
         nth_par(3) = 0.05
         nth_par(4) = 1.0
         nth_par(5) = 0.0
         Ifl=1
         call donthcomp(earx, nex, nth_par, ifl, contx_base, photer)

         !Now call with correct redshift
         nth_par(1) = real(Gamma)
         nth_par(2) = Cutoff_s
         nth_par(3) = 0.05
         nth_par(4) = 1.0
         nth_par(5) = ( 1.0 + zcos ) / real( gso(1) ) - 1.0
         Ifl=1
         call donthcomp(earx, nex, nth_par, ifl, contx, photer)
         
         !Calculate the ratio between the two at 1 keV
         dlogE  = log10(earx(nex)/earx(0)) / real(nex)
         i      = ceiling( -log10(earx(0)) / dlogE )
         j      = i + nint( log10((1.+zcos)/gso(1)) / dlogE )
         renorm = contx(i) / contx_base(j)

         !Take out the ridiculous nthcomp renormalisation
         contx = contx / renorm
         !Note: this is needed because nthcomp renormalises a redshifted spectrum to still
         !have a flux of 1 (or norm within xspec) at 1keV. We need to undo that
         !catastrophic monstrosity.
         
!The continuum needs to be renormalised according to the illuminating flux that was considered in xillver
!Plus we divide by a factor that depends on ionisation and density to agree with the first versions of reltrans
         do i = 1, nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            if (E .ge. 0.1 .and. E .le. 1e3) then
               Icomp = Icomp + ((earx(i) + earx(i-1)) * 0.5 * contx_base(i))
            endif
         enddo
         ! Renormalise (see header for derivation)
         contx = contx / ( 4.0*pi * ergsev * 1e5 * Icomp )
         
      else
         !First calculate normalisation constant in source frame
         do i = 1, nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            contx(i) = E**(-1.0*real(Gamma)+1) * exp(-E/(Cutoff_s))
            if (E .ge. 0.1 .and. E .le. 1e3) then
               Icomp = Icomp + ((earx(i) + earx(i-1)) * 0.5 * contx(i))
            end if
         end do
         !Now calculate spectrum in observer frame
         do i = 1, nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            contx(i) = E**(-1.0*real(Gamma)+1) * exp(-E/(Cutoff_obs))
         end do
         ! Renormalise (see header for derivation)
         contx = contx / ( 4.0*pi * ergsev * 1e5 * Icomp )
         
      end if
         
      return
    end subroutine getcont
