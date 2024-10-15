!-----------------------------------------------------------------------
      subroutine getcont(earx, nex, Gamma, Ecut_obs, logxi, logne, contx)
!!! Calculates continuum spectrum calling nthComp with the correct normalisation
!!!based on the xillver spectrum 
!!!  Arg:
        !  earx: energy grid
        !  nex: number of grid points
        !  Gamma: continuum spectrum inclination
        !  Ecut_obs: high energy cut-off or electron temperature
        !  logxi: ionisation parameter
        !  logne: density 
        !  (output) contx: continuum spectrum
!!! Last change: Gullo 2023 Nov; adapted to match the xillver tables call 
        !
      implicit none
      integer, intent(in)           :: nex
      real   , intent(in)           :: earx(0:nex), Ecut_obs, logxi, logne
      real   , intent(out)          :: contx(nex)
      double precision , intent(in) :: Gamma

      real   , parameter  :: ergsev  = 1.602197d-12 ! Convert eV to ergs
      integer :: i, ifl
      real    :: nth_par(5), photer(nex), E
      real    :: get_norm_cont
      
!So far this works only with kTe, so only with nthComp continuum model
      nth_par(1) = real(Gamma)
      nth_par(2) = Ecut_obs
      nth_par(3) = 0.05d0
      nth_par(4) = 1.d0
      nth_par(5) = 0.d0
      Ifl=1

      write(*,*) 'continuum parameters', nth_par, logxi, logne
      call donthcomp(earx, nex, nth_par, ifl, contx, photer)
!the continuum needs to be renormalised according to the illuminating flux that was considered in xillver get_norm_cont does the job
! Plus we divide by a factor that depends on ionisation and density to agree with the first versions of reltrans     
      contx = contx * get_norm_cont(real(Gamma), Ecut_obs, logxi, logne)
      contx = contx / 10**(logxi + logne - 15)
      return
    end subroutine getcont
!-----------------------------------------------------------------------


! !-----------------------------------------------------------------------
! subroutine getcont(nex,earx,Gamma,Ecut_obs,Cp,contx)
!     !!! Calculates continuum spectrum and sets xillver parameters !!!
!     !!!  Arg:
!     !  earx: energy grid
!     !  nex: number of grid points
!     !  Gamma: continuum spectrum inclination
!     !  Afe: iron abundance (not important for the continuum)
!     !  Ecut_obs: high energy cut-off or electron temperature
!     !  dens: log(density) of the disc
!     !  logxi: ionisation parameter (not important for the continuum)
!     !  Cp: sets which xillver spectrum
!     ! (output) contx: continuum spectrum
!     !!! Last change: Adam 2021 March; based on Gullo's getcont code from 2020 Jul
!     !
!     implicit none
!     integer, intent(in)  :: nex, Cp
!     real   , intent(in)  :: earx(0:nex), Ecut_obs
!     real   , intent(out) :: contx(nex)
!     real                 :: xillpar(7), xillparDCp(8)
!     double precision , intent(in) :: Gamma
!     integer :: ifl   
!     !Set the parameters for the continuum and not only 

!     !First determine which xillver and so which array of parameters 7 or 8
!     ! Remember that only xillverDCp has 8 par all the others have 7 but with different interpretation of parameter 3 (xillver(3)) 
!     !to clean up the code and not pass a million arguments we just take default values for ne/csi/theta/Afe, it does not impact
!     !the continuum anyway 
!     if (Cp .eq. 2) then
!         xillparDCp(1) = real( Gamma )
!         xillparDCp(2) = 1.
!         xillparDCp(3) = Ecut_obs
!         xillparDCp(4) = 15.
!         xillparDCp(5) = 0.0
!         xillparDCp(6) = 0.0   
!         xillparDCp(7) = 30.0 
!         xillparDCp(8) = 0.0       !reflection fraction of 0
!     else if (Cp .lt. 0 ) then
!         xillpar(3) = Ecut_obs
!     else if (Cp .eq. 1) then
!         xillpar(3) = 15.
!     endif

!     xillpar(1) = real( Gamma )
!     xillpar(2) = 1.
!     xillpar(4) = 0.
!     xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
!     xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
!     xillpar(7) = 0.0       !reflection fraction of 0
!     ifl        = 1

!     call get_xillver(earx,nex,xillpar,xillparDCp,ifl,Cp,contx)
!     return
! end subroutine getcont
! !-----------------------------------------------------------------------
