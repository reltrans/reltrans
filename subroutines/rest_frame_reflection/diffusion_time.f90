module m_diffusion_time
!> Diffusion (reprocessing) time of the reflected photons.
!>
!> Reflected photons random-walk inside the disc atmosphere before escaping,
!> so the reflection responds to the illuminating continuum with an extra
!> delay tau(E) on top of the light-travel time. In the Fourier domain this
!> multiplies the rest-frame reflection spectrum by exp(i 2 pi f tau(E)),
!> which makes it complex; the convolution with the relativistic kernel is
!> then done with complex-to-complex FFTs (conv_mod: spectrum_FT_c2c and
!> conv_FT_c2c).
!>
!> Currently the delay is energy independent and set by the diff_time
!> parameter (seconds, rest frame); the (1+z) factor converts it to the
!> observer frame, consistent with the light-travel delays in strans.
!>
!> The delay is returned on the internal (rest-frame) energy grid, so that an
!> energy-dependent diffusion time (e.g. from absorption opacities) can be
!> added here later without touching the convolution.
    use rtconstants, only: wp
    use common_types, only: t_model_arguments
    implicit none
    private
    public :: is_diffusion_time_active, get_diffusion_delay

    ! G M_sun / c^3 in seconds (same constant used for the frequency grid in config_frequency).
    real(wp), parameter :: GM_over_C3 = 4.92695275718945e-06_wp

contains

    logical function is_diffusion_time_active(model_args, dc) result(active)
!> Whether the diffusion time has to be included in the convolution. 
!> Inputs:
!>     model_args : t_model_arguments, model parameters
!>     dc         : integer, 1 for the time-averaged spectrum, 0 otherwise
!> Outputs:
!>     active     : logical
        type(t_model_arguments), intent(in) :: model_args
        integer                , intent(in) :: dc
        active = (dc .eq. 0) .and. (model_args%diff_time .gt. 0.0_wp)
    end function is_diffusion_time_active

    subroutine get_diffusion_delay(model_args, tau_diff)
!> Diffusion time of the reflected photons on the internal (rest-frame)
!> energy grid, in units of Rg/c in the observer frame (the same units as
!> the light-travel delays and frequencies used in the kernel calculation).
!> Inputs:
!>     model_args : t_model_arguments, model parameters (diff_time, Mass,
!>                  zcos)
!> Outputs:
!>     tau_diff   : real(:), diffusion delay in Rg/c for each energy bin
        type(t_model_arguments), intent(in)  :: model_args
        real(wp)               , intent(out) :: tau_diff(:)

        tau_diff = (1.0_wp + model_args%zcos) * model_args%diff_time /         &
            (GM_over_C3 * model_args%Mass)
    end subroutine get_diffusion_delay

end module m_diffusion_time
