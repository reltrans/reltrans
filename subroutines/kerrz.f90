module kerrz
    use kerrz_interface
    use rtconstants, only: pi

    implicit none

    ! This effective infinity r coordinate is for compatability with what YNOGK
    ! used to do.
    double precision, parameter :: R_AT_INFINITY = 18000000.0d0

    ! This is a singleton metric so that it does not need to be recalculated often.
    ! TODO: move this to `t_config` or related.
    type(krz_KerrMetric) :: kerr_metric

contains

    type(krz_TraceResult) function trace_impact_parameters(mu_obs, alpha, beta,&
        distance) result(res)
        !> Trace a photon on an image plane parameterised by some impact
        !> parameters `alpha` and `beta`. The image plane is assumed to be at
        !> infinity, though the exact distance can be set with the optional
        !> `distance` parameter.
        double precision, intent(in) :: mu_obs, alpha, beta
        double precision, optional, intent(in) :: distance
        double precision :: dist
        type(krz_InitialConditions) :: ic
        type(krz_FourVector) :: x_obs

        ! The default distance is `R_AT_INFINITY` unless provided by the user
        if (present(distance)) dist = distance
        if (.not. present(distance)) dist = R_AT_INFINITY

        x_obs = krz_FourVector(t = 0.0d0, r = dist, th = acos(mu_obs),         &
            ph = 0.0d0)

        ! TODO: remove this once kerrz has fully face-on implemented
        if (abs(x_obs%th) < 1d-4) then
            x_obs%th = 1d-4
        end if

        ic = krz_fromImpactParameters(kerr_metric, x_obs, alpha, beta)
        res = krz_traceToAngle(kerr_metric, ic, pi / 2.0)
    end function trace_impact_parameters

    type(krz_TraceResult) function trace_lamppost(h, delta_s)       &
        result(res)
        !> Trace a photon from a lamppost at some height `h` to the disc with an
        !> initial inclination angle `delta_s` in the local sky of the lamppost.
        !>
        !> As in Ingram et al. 2019, the convention is that `delta_s = 0` traces
        !> directly downwards towards the black hole, and `delta_s = pi` directly
        !> upwards towards infinity.
        double precision, intent(in) :: h, delta_s
        type(krz_InitialConditions) :: ic
        type(krz_FourVector) :: x
        type(krz_OrthonormalFrame) :: frame

        ! Theta in the four-position cannot be directly on the spin axis
        ! currently in kerrz so it is marginally offset here to avoid errors.
        ! TODO: cache the frame and pass it in as an argument to avoid
        ! recomputation.
        x = krz_FourVector(t=0.0d0, r=h, th=1.0d-5, ph=0.0d0)
        frame = krz_stationaryFrame(kerr_metric, x)
        ic = krz_fromSkyAngles(kerr_metric, frame, delta_s - pi, 0.0d0)
        res = krz_traceToAngle(kerr_metric, ic, pi / 2.0)
    end function trace_lamppost

    ! These subroutines are defined for the test suite:
    subroutine test_kerrz_trace(spin, mu_obs, alpha, beta, t, r, theta, phi)   &
        bind(C, name="test_kerrz_trace")
        double precision, intent(in) :: spin, mu_obs, alpha, beta
        double precision, intent(out) :: t, r, theta, phi
        type(krz_TraceResult) :: res
        kerr_metric = krz_KerrMetric_init(1.0d0, spin)
        res = trace_impact_parameters(mu_obs, alpha, beta)
        t = res%x_final%t
        r = res%x_final%r
        theta = res%x_final%th
        phi = res%x_final%ph
    end subroutine test_kerrz_trace

    subroutine test_kerrz_trace_lamppost(spin, h, delta_s, t, r, theta, phi)   &
        bind(C, name="test_kerrz_trace_lamppost")
        double precision, intent(in) :: spin, h, delta_s
        double precision, intent(out) :: t, r, theta, phi
        type(krz_TraceResult) :: res
        kerr_metric = krz_KerrMetric_init(1.0d0, spin)
        res = trace_lamppost(h, delta_s)
        t = res%x_final%t
        r = res%x_final%r
        theta = res%x_final%th
        phi = res%x_final%ph
    end subroutine test_kerrz_trace_lamppost

end module kerrz
