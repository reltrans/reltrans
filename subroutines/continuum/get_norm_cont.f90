!-----------------------------------------------------------------------
    function get_norm_cont(Gamma, Ecut_obs, logxi, logne)
      implicit none

      real   , intent(IN) :: logxi, logne, Gamma, Ecut_obs
      integer, parameter  :: nex = 50
      real   , parameter  :: pi = acos(-1.0), ergsev  = 1.602197e-12 ! Convert eV to ergs
      integer             :: i, ifl
      real                :: Icomp, inc_flux, Emin, Emax, dloge 
      real                :: nth_par(5), energy(0:nex), spec(nex), photer(nex)
      real                :: get_norm_cont
            
!continuum paramters 
      nth_par(1) = Gamma
      nth_par(2) = Ecut_obs
      nth_par(3) = 0.05d0
      nth_par(4) = 1.d0
      nth_par(5) = 0.d0

!define the energy grid
      Emin = 1e-1
      Emax = 1000    
      dloge = log10( Emax / Emin ) / float(nex)
      do i = 0, nex
         energy(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
      end do

      Ifl=1      
      call donthcomp(energy, nex, nth_par, ifl, spec, photer)

      Icomp = 0.0
      do i = 1, nex
         Icomp = Icomp + ((energy(i) + energy(i-1)) * 0.5 * spec(i)) 
      enddo
      inc_flux = 10**(logne + logxi) / (4.0 * pi) / (ergsev * 1000.0) !calculate incident flux in units  [keV/cm^2/s]
      get_norm_cont = inc_flux/ Icomp / 1e20

      
    end function get_norm_cont
