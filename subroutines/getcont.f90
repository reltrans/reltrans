!-----------------------------------------------------------------------
      subroutine getcontDCp(nex, earx, Gamma, Afe, kTe, dens, logxi, contx, xillparDCp)
! Calculates continuum spectrum
      implicit none
      integer,          intent(in)  :: nex
      integer                       :: ifl 
      real,             intent(in)  :: earx(0:nex), Afe, logxi, dens, kTe
      real,             intent(out) :: xillparDCp(8),  contx(nex)
      double precision, intent(in)  :: Gamma
        
      !First continuum
      xillparDCp(1) = real( Gamma )
      xillparDCp(2) = Afe
      xillparDCp(3) = kTe
      xillparDCp(4) = dens
      xillparDCp(5) = logxi
      xillparDCp(6) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillparDCp(7) = 30.0       !inclination angle (doesn't matter for continuum)
      xillparDCp(8) = 0.0       !reflection fraction of 0

      ifl        = 1

      call myxillDCp(earx, nex, xillparDCp, ifl, contx)
      return
    end subroutine getcontDCp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine getcont(nex, earx, Gamma, Afe, kTe, Ecut_obs, dens, logxi, Cp, contx, xillpar, xillparDCp)
!!! Calculates continuum spectrum and sets xillver parameters !!!
!!!  Arg:
        !  earx: energy grid
        !  nex: number of grid points
        !  Gamma: continuum spectrum inclination
        !  Afe: iron abundance (not important for the continuum)
        !  kTe: electron temperature (illuminating spectrum: nthcomp)
        !  Ecut_obs: high energy cut-off (illuminating spectrum: powerlaw)
        !  dens: density of the disc
        !  logxi: ionisation parameter (not important for the continuum)
        !  Cp: sets which xillver spectrum
        ! (output) contx: continuum spectrum
        ! (output) xillpar: xillver parameters
        ! (output) xillparDCp: xillverDCp parameters
!!! Last change: Gullo - 2020 Jul
        
      implicit none
      integer, intent(in)  :: nex, Cp
      real   , intent(in)  :: earx(0:nex), Afe, kTe, Ecut_obs, logxi, dens
      real   , intent(out) :: contx(nex), xillpar(7), xillparDCp(8)
      double precision , intent(in) :: Gamma
      integer :: ifl   
!Set the parameters for the continuum and not only 

!First determine which xillver and so which array of parameters 7 or 8
! Remember that only xillverDCp has 8 par all the others have 7 but with different interpretation of parameter 3 (xillver(3))  
      if (Cp .eq. 2) then
         xillparDCp(1) = real( Gamma )
         xillparDCp(2) = Afe
         xillparDCp(3) = kTe
         xillparDCp(4) = dens
         xillparDCp(5) = logxi
         xillparDCp(6) = 0.0   
         xillparDCp(7) = 30.0 
         xillparDCp(8) = 0.0       !reflection fraction of 0
      else if (Cp .eq. -1) then
         xillpar(3) = Ecut_obs
      else if (Cp .eq. -2) then
         xillpar(3) = kTe
      else if (Cp .eq. 1) then
         xillpar(3) = dens
      endif

      xillpar(1) = real( Gamma )
      xillpar(2) = Afe
      xillpar(4) = logxi
      xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
      xillpar(7) = 0.0       !reflection fraction of 0
      ifl        = 1

      call myxill(earx, nex, xillpar, xillparDCp, ifl, Cp, contx)
      return
      end subroutine getcont
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine ad_getcont(nex, earx, Gamma, Afe, Ecut_obs, dens, logxi, Cp, contx)
!!! Calculates continuum spectrum and sets xillver parameters !!!
!!!  Arg:
        !  earx: energy grid
        !  nex: number of grid points
        !  Gamma: continuum spectrum inclination
        !  Afe: iron abundance (not important for the continuum)
        !  Ecut_obs: high energy cut-off or electron temperature
        !  dens: log(density) of the disc
        !  logxi: ionisation parameter (not important for the continuum)
        !  Cp: sets which xillver spectrum
        ! (output) contx: continuum spectrum
!!! Last change: Adam 2021 March; based on Gullo's getcont code from 2020 Jul
        !
      implicit none
      integer, intent(in)  :: nex, Cp
      real   , intent(in)  :: earx(0:nex), Afe, Ecut_obs, logxi, dens
      real   , intent(out) :: contx(nex)
      real                 :: xillpar(7), xillparDCp(8)
      double precision , intent(in) :: Gamma
      integer :: ifl   
!Set the parameters for the continuum and not only 

!First determine which xillver and so which array of parameters 7 or 8
! Remember that only xillverDCp has 8 par all the others have 7 but with different interpretation of parameter 3 (xillver(3))  
      if (Cp .eq. 2) then
         xillparDCp(1) = real( Gamma )
         xillparDCp(2) = Afe
         xillparDCp(3) = Ecut_obs
         xillparDCp(4) = dens
         xillparDCp(5) = logxi
         xillparDCp(6) = 0.0   
         xillparDCp(7) = 30.0 
         xillparDCp(8) = 0.0       !reflection fraction of 0
      else if (Cp .lt. 0 ) then
         xillpar(3) = Ecut_obs
      else if (Cp .eq. 1) then
         xillpar(3) = dens
      endif

      xillpar(1) = real( Gamma )
      xillpar(2) = Afe
      xillpar(4) = logxi
      xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
      xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
      xillpar(7) = 0.0       !reflection fraction of 0
      ifl        = 1

      call myxill(earx, nex, xillpar, xillparDCp, ifl, Cp, contx)
      return
      end subroutine ad_getcont
!-----------------------------------------------------------------------
