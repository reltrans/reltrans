!-----------------------------------------------------------------------
    subroutine get_xillver(ear, ne, dim, dimCp, param_xillPL, param_xillCp, Cp, photar)
!!! This routine calls the correct version of xillver based on the Cp !!!
!!!   Arg:
      !  ear: energy grid
      !  ne: number of grid points
      !  param6: array of xillver parameters for xillver, xillverD
      !  param7: array of xillver parameters for xillverDCp
      !  ifl : parameter to call xspec model
      !  Cp : chooses xillver model
      !      -1 xillver      1e15 density and powerlaw illumination  
      !       1 xillverD     high density and powerlaw illumination
      !       2 xillverDCp   high density and nthcomp  illumination
      ! photar: (output) xillver energy spectrum

      use rtconstants, only: wp
      use xillver_tables
      use xspec_interface, only: table_model
      implicit none
      integer, intent(in)  :: ne, Cp, dim, dimCp
      real(wp) , intent(in)  :: ear(0:ne), param_xillPL(dim), param_xillCp(dimCp)
      real(wp) , intent(out) :: photar(ne)

      integer             :: ifl, i
      real(wp)            :: m

      ifl = 0
      if( Cp .eq. -1 )then         !xillver
         ! Pass the null terminator for C compatability, as xsatbl is an
         ! external C function.
         call table_model(ear, ne, param_xillPL,                               &
              [trim(pathname_xillver), char(0)], ifl, photar)
      else if ( Cp .eq. 2 )then    !xillverDCp
         ! Pass the null terminator for C compatability, as xsatbl is an
         ! external C function.
         call table_model(ear, ne, param_xillCp,                               &
            [trim(pathname_xillverDCp), char(0)], ifl, photar)
      else
         write(*,*) 'No xillver model available for this configuration'
         stop 
      end if


! The internal energy grid runs from 0.01 keV while the xillver tables
! are computes from 0.071770 keV (lower edge of the first bin) to 1000.1 keV
! When we ask xsatbl to return the spectrum in an energy range where
! it is not computed it returns zero. This causes a sudden drop in the
! spectral flux which can cause some problem to the FFT and the emulator
! Here we avoid it with a patch:
! 0.071770 keV falls inside internal bin 641 (0.07175–0.07197 keV), thus
! bin 642 is the first bin with non-zero flux
      m = 0.1_wp
      do i = 1, 641
         photar(i) = m * real(i, wp) + photar(642) - (m * 642.0_wp)
         if (photar(i) .lt. 0.0_wp) then
            photar(i) = 0.0_wp
         endif
      end do

      return
    end subroutine get_xillver
!-----------------------------------------------------------------------
