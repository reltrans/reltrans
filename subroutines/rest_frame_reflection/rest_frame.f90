!-----------------------------------------------------------------------
subroutine rest_frame(model_args, arrays, Gamma, logne, Cutoff, logxi,         &
    thetae, photar)
!
!  model_args%Cp : chooses reflection model
!      -1 xillver      1e15 density and powerlaw illumination
!       1 xillverD     high density and powerlaw illumination
!       2 xillverDCp   high density and nthcomp  illumination
!       0 reflionxDCp  reflionx high density and nthcomp  illumination
!
!       The first 2 have the same number of parameters xillpar(6),
!       in the first one Cutoff is a parameter (either energy or temperature)
!       and the density is fixed to 10^15
!       in the second one the density is a parameter and Cutoff is fixed to
!       300 keV
!       The Cp = 2 has one more parameter both Cutoff and density are
!       parameters

!       Last change: Gullo - 2022 Oct

   use common_types
   use conv_mod, only: nex
   implicit none
   type(t_model_arguments), intent(in) :: model_args
   type(t_arrays),          intent(in) :: arrays
   real, intent(in)    :: Gamma, logne, Cutoff, logxi, thetae
   real, intent(out)   :: photar(nex)
   integer, parameter  :: dim = 6, dimCp = 8
   real                :: xillpar(dim), xillparDCp(dimCp)

   if( model_args%Cp .ne. 0 )then
      write(*,*) "DEBUG: model_args%Cp", model_args%Cp
      !The model is a xillver model
      !Fill parameter arrays
      xillpar(1) = Gamma                !Power law index
      xillpar(2) = model_args%Afe       !Iron abundance
      xillpar(3) = logxi                !ionization par
      xillpar(4) = Cutoff               !Ecut or kTe
      if( model_args%Cp .eq. 1 )then
         xillpar(4) = logne             !logne
      end if
      xillpar(5) = thetae               !emission angle
      xillpar(6) = 0.0                  !redshift
      xillparDCp(1) = Gamma             !photon index
      xillparDCp(2) = model_args%Afe    !Afe
      xillparDCp(3) = logxi             !ionization par
      xillparDCp(4) = Cutoff            !kTe
      xillparDCp(5) = logne             !logne
      xillparDCp(6) = thetae            !emission angle
      xillparDCp(7) = 0.0               !redshift

      call get_xillver(arrays%earx, nex, dim, dimCp, xillpar, xillparDCp,      &
          model_args%Cp, photar)

   else
      !The model is reflionx
      call normreflionx(arrays%earx, nex, Gamma, model_args%Afe, logne,        &
          Cutoff, logxi, thetae, photar)
   end if

   return
 end subroutine rest_frame
!-----------------------------------------------------------------------


