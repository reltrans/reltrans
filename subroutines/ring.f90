! The radial profiles are currently all averaged from the illumination profiles.
! This averaging is performed by the `ring_average` function below, and
! populates the `dyn_gr` module arrays. This is to minimise the number of
! changes that need to be make to reltrans at this time.
!
! The following functions are therefore not fully self-consistent in the time
! resolved spectra:
!
!  - radial_profiles/logxiraw.f90
!  - radial_profiles/xiraw.f90
!  - radial_profiles/radfunctions_dens.f90
!  - radial_profiles/radfunctions_dist.f90
module ring_corona
implicit none

    ! The number of bins in coronal ring azimuth to sum over.
    integer, parameter :: ring_nphi = 50
    ! The weighting factor needed when summing over the different azimuthal
    ! bins.
    double precision, parameter :: ring_weight = 1.0d0 / float(ring_nphi)

    ! This normalisation is determined so that the emissivity profiles that
    ! kerrz calculates match what reltrans calculates using the internal method.
    double precision, parameter :: ring_em_normalisation = 1.0d0 / (4d0 * acos(-1d0))

contains

    double precision function mean_r_annulus(rmax, rmin, n, i) result(mean_r)
        double precision, intent(in) :: rmax, rmin
        integer, intent(in) :: n, i
        mean_r = (rmax/rmin)**(real(i-1) / real(n))
        mean_r = mean_r + (rmax/rmin)**(real(i) / real(n))
        mean_r = mean_r * rmin * 0.5
    end function mean_r_annulus

    subroutine ring_radial_gradients(config, model_args, arrays)
        !> The ring-corona specific version of `radfuncs_dens`. The main
        !> difference here is that the parameterisation cannot be in terms of
        !> `cosd`, as that is not uniquely defined for each radius on the disc
        !> for a ring like corona. Instead, as kerrz has calculated the
        !> emissivity profile out to a large radial distance, those values are
        !> used, and the Newtonian extrapolation is performed over disc radius.
        !>
        !> Side effects are all arrays in `radial_grids`:
        !>
        !>      logxir, gsdr, logner, dfer_arr
        use common_types, only: t_config, t_model_arguments, t_arrays
        use radial_grids, only: logxir, gsdr, logner, dfer_arr
        use env_variables, only : adensity
        use kerrz, only: krz_EmissivityTrace, emissivity_values_at,            &
            krz_ContinuumRingPoint, ring_continuum_at
        use rtconstants, only: pi
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(in) :: model_args
        type(t_arrays), intent(in) :: arrays
        double precision :: radius, total_ionisation, g_so, phi, ionisation
        double precision :: g_sd_weighted, logxinorm, lognenorm, cos_incident_angle
        type(krz_EmissivityTrace) :: source_to_disc
        type(krz_ContinuumRingPoint) :: source_to_obs
        integer :: i, phi_i
        ! TODO: this is a function and should be in a module:
        double precision :: mylogne, dgsofac, xiraw

        do i = 1, config%xe
            radius = mean_r_annulus(config%rnmax, model_args%rin, config%xe, i)
            total_ionisation = 0.0
            g_sd_weighted = 0.0

            ! Calculate the raw-density (only matters for high density
            ! reltransD)
            logner(i) = adensity * mylogne(radius, model_args%rin)

            ! Loop over every azimuthal bin of the ring-corona:
            do phi_i = 1, ring_nphi
                phi = 2 * pi * phi_i * ring_weight
                source_to_disc = emissivity_values_at(radius, phi)
                source_to_obs = ring_continuum_at(phi)
                g_so = source_to_obs%energyshift

                cos_incident_angle = cos(source_to_disc%local_theta)

                ! TODO: double check this later, but I'm pretty sure the `xiraw`
                ! function is just the emissivity profile.
                ionisation = source_to_disc%em * ring_weight                   &
                    * ring_em_normalisation

                ! Correction to account for the radial dependence of incident
                ! angle, and for the g factors
                ! TODO: check where contx_int gets calculated.
                ionisation = ionisation * arrays%contx_int(1)                  &
                    * (g_so)**(model_args%Gamma - 2)      &
                    / (sqrt(2.0) * cos_incident_angle)
                total_ionisation = total_ionisation + ionisation

                ! Accumulate the energyshift factor weighted by the ionisation:
                g_sd_weighted = g_sd_weighted + ionisation * source_to_disc%g
            end do

            ! Now compute the azimuthally averaged values:
            gsdr(i) = g_sd_weighted / total_ionisation

            logxir(i) = log10(total_ionisation) - logner(i)
        end do

        ! Then we can determine the max and min ionisation, which allows it to
        ! be renormalised:
        logxinorm = maxval(logxir)
        lognenorm = maxval(logner)
        logxir = logxir - (logxinorm - dble(model_args%logxi))
        logner = logner - (lognenorm - dble(model_args%lognep))

    end subroutine ring_radial_gradients

    subroutine ring_cross_spectrum_with_continuum(config, model_args, arrays,  &
        nex)
        ! The Fourier transform of the disc and continuum looks as follows:
        !
        !     s(nu, E) = A(nu) [ C(nu, E) + R(nu, E) ]
        !
        use common_types, only: t_config, t_model_arguments, t_arrays
        use kerrz, only: krz_ContinuumRingPoint, ring_continuum_at
        use rtconstants, only: pi
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(in) :: model_args
        type(t_arrays), intent(inout) :: arrays
        integer, intent(in) :: nex
        type(krz_ContinuumRingPoint) :: source_to_obs
        integer :: j, i, phi_i
        double precision :: f, E, fac, phase, phi, phi_
        complex :: W0, W1, W2, W3, Stemp, c_exp

        do j=1,config%nf
            f = config%flo *                                                   &
                (config%fhi/config%flo)**((real(j)-0.5d0) / real(config%nf))
            do i=1,nex
                E = 0.5d0 * ( arrays%earx(i) + arrays%earx(i-1) )
                ! Set up transfer functions
                W0 = model_args%boost *                                        &
                    cmplx(arrays%ReW0(1,i,j),arrays%ImW0(1,i,j))
                W1 = (1-config%DC) * model_args%boost *                        &
                    cmplx(arrays%ReW1(1,i,j),arrays%ImW1(1,i,j))
                W2 = (1-config%DC) * model_args%boost *                        &
                    cmplx(arrays%ReW2(1,i,j),arrays%ImW2(1,i,j))
                W3 = config%ionvar * (1-config%DC) * model_args%boost *        &
                    cmplx(arrays%ReW3(1,i,j),arrays%ImW3(1,i,j))

                ! This represents the disc response:
                Stemp = model_args%g(1) * (W1 + W2)
                Stemp = Stemp + W0 + W3

                arrays%ReSraw(i,j) = arrays%ReSraw(i,j) + real(Stemp)
                arrays%ImSraw(i,j) = arrays%ImSraw(i,j) + aimag(Stemp)

                ! Then sum the coronal contribution:
                do phi_i = 1, ring_nphi
                    phi = 2 * pi * phi_i * ring_weight
                    source_to_obs = ring_continuum_at(phi)

                    phase = 2 * pi * source_to_obs%delta_t * f
                    c_exp = cmplx(cos(phase), sin(phase))
                    c_exp = cmplx(1.0, 0.0)

                    fac = log(source_to_obs%energyshift/((1d0+model_args%zcos)*E))

                    Stemp = (1 + fac * model_args%g(1)) * c_exp *              &
                        arrays%ring_continuums(i,phi_i) * ring_weight

                    ! Separate into real/imaginary parts for compatibility with the
                    ! rest of the code
                    arrays%ReSraw(i,j) = arrays%ReSraw(i,j) + real(Stemp)
                    arrays%ImSraw(i,j) = arrays%ImSraw(i,j) + aimag(Stemp)
                 end do
             end do
        end do
    end subroutine ring_cross_spectrum_with_continuum

    subroutine ring_calculate_continuum(config, model_args, arrays)
        use common_types, only: t_config, t_model_arguments, t_arrays
        use kerrz, only: krz_ContinuumRingPoint, ring_continuum_at
        use rtconstants, only: pi
        use gr_continuum, only: gso, lens
        use conv_mod, only: nex
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(inout) :: model_args
        type(t_arrays), intent(inout) :: arrays
        integer :: phi_i
        double precision :: old_gso, aggregate_cutoff, phi
        type(krz_ContinuumRingPoint) :: source_to_obs
        ! Store the old `gso` value so it can be restored afterwards. This is
        ! the 'average' source-to-disc energyshift factor.
        old_gso = gso(1)

        do phi_i = 1, ring_nphi
            phi = 2 * pi * phi_i * ring_weight
            source_to_obs = ring_continuum_at(phi)

            ! TODO: make gso a parameter of getcont
            ! gso(1) = source_to_obs%energyshift

            model_args%Cutoff_obs = model_args%Cutoff_s * gso(1)               &
                / real(1.d0 + model_args%zcos)

            ! Sum the aggregate cutoff, so there is some effective 'seen'
            ! cutoff that can be assigned back later.
            aggregate_cutoff = aggregate_cutoff                                &
                + model_args%Cutoff_obs * ring_weight

            call getcont(2, arrays%earx, nex, model_args%Gamma,                &
                model_args%Cutoff_s, model_args%Cutoff_obs, model_args%logxi,  &
                model_args%lognep, model_args%zcos,                            &
                arrays%ring_continuums(:,phi_i))

            arrays%ring_continuums(:,phi_i) = lens(1)                          &
                * (gso(1) / (real(1.d0 + model_args%zcos)))                    &
                * arrays%ring_continuums(:,phi_i)

            ! This aggregate is stil needed for the ionisation `xilimits`.
            arrays%contx(:,1) = arrays%contx(:,1)                              &
                + arrays%ring_continuums(:,phi_i) * ring_weight
        end do

        gso(1) = old_gso
        ! model_args%Cutoff_obs = aggregate_cutoff
    end subroutine ring_calculate_continuum

end module ring_corona

