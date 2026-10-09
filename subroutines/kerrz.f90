module kerrz
    use kerrz_interface
    use rtconstants, only: pi, wp
    use iso_c_binding, only: c_double

    implicit none

    ! This effective infinity r coordinate is for compatability with what YNOGK
    ! used to do.
    real(wp), parameter :: R_AT_INFINITY = 1.0e5_wp

    type :: LamppostContinuum
        ! This is |∂cosδ / ∂cosθ|, the lensing factor.
        real(wp) :: lensing_factor
        ! The cosine of the angle on the local sky, with δ=0 pointing towards
        ! the black hole.
        real(wp) :: cos_delta
        ! The corona-to-observer time.
        real(wp) :: time
        ! The impact parameters for this geodesic.
        real(wp) :: alpha, beta
    end type LamppostContinuum

contains

    type(krz_TraceResult) function trace_impact_parameters(metric, mu_obs,     &
        alpha, beta, distance) result(res)
        !> Trace a photon on an image plane parameterised by some impact
        !> parameters `alpha` and `beta`. The image plane is assumed to be at
        !> infinity, though the exact distance can be set with the optional
        !> `distance` parameter.
        type(krz_KerrMetric), intent(in) :: metric
        real(wp), intent(in) :: mu_obs, alpha, beta
        real(wp), optional, intent(in) :: distance
        real(wp) :: dist
        type(krz_InitialConditions) :: ic
        type(krz_FourVector) :: x_obs

        ! The default distance is `R_AT_INFINITY` unless provided by the user
        if (present(distance)) dist = distance
        if (.not. present(distance)) dist = R_AT_INFINITY

        x_obs = krz_FourVector(t = 0.0_wp, r = dist, th = acos(mu_obs),        &
            ph = 0.0_wp)

        ! TODO: remove this once kerrz has fully face-on implemented
        if (abs(x_obs%th) < 1e-3_wp) then
            x_obs%th = 1e-3_wp
        end if

        ic = krz_fromImpactParameters(metric, x_obs, alpha, beta)
        res = krz_traceToAngle(metric, ic, pi / 2.0_wp)
    end function trace_impact_parameters

    type(krz_TraceResult) function trace_lamppost(metric, h, delta_s)          &
        result(res)
        !> Trace a photon from a lamppost at some height `h` to the disc with an
        !> initial inclination angle `delta_s` in the local sky of the lamppost.
        !>
        !> As in Ingram et al. 2019, the convention is that `delta_s = 0` traces
        !> directly downwards towards the black hole, and `delta_s = pi` directly
        !> upwards towards infinity.
        type(krz_KerrMetric), intent(in) :: metric
        real(wp), intent(in) :: h, delta_s
        type(krz_InitialConditions) :: ic
        type(krz_FourVector) :: x
        type(krz_OrthonormalFrame) :: frame

        ! Theta in the four-position cannot be directly on the spin axis
        ! currently in kerrz so it is marginally offset here to avoid errors.
        ! TODO: cache the frame and pass it in as an argument to avoid
        ! recomputation.
        x = krz_FourVector(t=0.0_wp, r=h, th=1.0e-5_wp, ph=0.0_wp)
        frame = krz_stationaryFrame(metric, x)
        ic = krz_fromSkyAngles(metric, frame, delta_s - pi, 0.0_wp)
        res = krz_traceToAngle(metric, ic, pi / 2.0_wp)
    end function trace_lamppost

    type(LamppostContinuum) function trace_lensing(metric, h, r_obs, mu_obs)   &
        result(cont)
        type(krz_KerrMetric), intent(in) :: metric
        real(wp), intent(in) :: h, r_obs, mu_obs
        type(krz_ContinuumLamppost) :: continuum
        type(krz_FourVector) :: x

        x = krz_FourVector(t=0.0_wp, r=r_obs, th=acos(mu_obs), ph=0.0_wp)

        ! TODO: remove this once kerrz has fully face-on implemented
        if (abs(x%th) < 1e-3_wp) then
            x%th = 1e-3_wp
        end if

        continuum = krz_traceContinuumLamppost(metric, x, h, 0.0_wp)

        ! Note the angle mapping to be consistent with the Reltrans convention.
        ! Also the sign change on beta.
        cont = LamppostContinuum(lensing_factor=1.0_wp/continuum%dcosd_dcosth, &
            cos_delta = cos(pi - continuum%angle_delta),                       &
            time = continuum%res%x_final%t, alpha = continuum%alpha,           &
            beta = -continuum%beta)
    end function trace_lensing

    ! These subroutines are defined for the test suite:
    subroutine test_kerrz_trace(spin, mu_obs, alpha, beta, t, r, theta, phi)   &
        bind(C, name="test_kerrz_trace")
        real(c_double), intent(in) :: spin, mu_obs, alpha, beta
        real(c_double), intent(out) :: t, r, theta, phi
        type(krz_KerrMetric) :: metric
        type(krz_TraceResult) :: res
        metric = krz_KerrMetric_init(1.0_wp, spin)
        res = trace_impact_parameters(metric, mu_obs, alpha, beta)
        t = res%x_final%t
        r = res%x_final%r
        theta = res%x_final%th
        phi = res%x_final%ph
    end subroutine test_kerrz_trace

    subroutine test_kerrz_trace_lamppost(spin, h, delta_s, t, r, theta, phi)   &
        bind(C, name="test_kerrz_trace_lamppost")
        real(c_double), intent(in) :: spin, h, delta_s
        real(c_double), intent(out) :: t, r, theta, phi
        type(krz_KerrMetric) :: metric
        type(krz_TraceResult) :: res
        metric = krz_KerrMetric_init(1.0_wp, spin)
        res = trace_lamppost(metric, h, delta_s)
        t = res%x_final%t
        r = res%x_final%r
        theta = res%x_final%th
        phi = res%x_final%ph
    end subroutine test_kerrz_trace_lamppost

    subroutine test_kerrz_lensing(spin, h, r_obs, mu_obs, lensing_factor,      &
        cos_delta, time) bind(C, name="test_kerrz_lensing")
        real(c_double), intent(in) :: spin, h, r_obs, mu_obs
        real(c_double), intent(out) :: lensing_factor, cos_delta, time
        type(krz_KerrMetric) :: metric
        type(LamppostContinuum) :: cont
        metric = krz_KerrMetric_init(1.0_wp, spin)
        cont = trace_lensing(metric, h, r_obs, mu_obs)
        lensing_factor = cont%lensing_factor
        cos_delta = cont%cos_delta
        time = cont%time
    end subroutine test_kerrz_lensing

end module kerrz
