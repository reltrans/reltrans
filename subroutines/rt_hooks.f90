! C-callable hooks for the B1 investigation. Not for upstream.

! Trace a batch of image-plane rays. status = 0 if the ray reaches the disc plane.
subroutine rt_trace_batch(n, spin, mu0, alpha, beta, r, taudo, status) bind(C, name="rt_trace_batch")
    use iso_c_binding
    use kerrz, only: kerr_metric, krz_KerrMetric_init, krz_TraceResult,        &
        trace_impact_parameters, KRZ_STATUS_NONE
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: spin, mu0, alpha(n), beta(n)
    real(c_double), intent(out) :: r(n), taudo(n)
    integer(c_int), intent(out) :: status(n)
    double precision, parameter :: d = 1.0d5
    type(krz_TraceResult) :: res
    integer :: i
    kerr_metric = krz_KerrMetric_init(1.0d0, spin)
    do i = 1, n
        res = trace_impact_parameters(mu0, alpha(i), beta(i), d)
        if (res%status == KRZ_STATUS_NONE) then
            status(i) = 0
            r(i) = res%x_final%r
            taudo(i) = res%x_final%t - d
        else
            status(i) = 1
            r(i) = -1.d0
            taudo(i) = 0.d0
        end if
    end do
end subroutine rt_trace_batch

! Per-sample quantities for a single lamppost, reproducing sum_multiple_lampposts.
! Requires a preceding model call (state saved by rt_store_and_override).
! Outputs: g, dfe = emissivity * (g/(1+z))^(2+Gamma) per unit solid angle,
! tau, lgsd = log(gsd), mue.
subroutine rt_sample_quantities(n, nonrel, alpha, beta, re, taudo, g, dfe, tau, lgsd, mue) &
        bind(C, name="rt_sample_quantities")
    use iso_c_binding
    use rt_state
    use dyn_gr
    use radial_grids, only: pnorm
    use gr_continuum, only: tauso
    use emissivities
    use kerrz, only: kerr_metric, krz_KerrMetric_init
    implicit none
    integer(c_int), intent(in) :: n, nonrel
    real(c_double), intent(in) :: alpha(n), beta(n), re(n), taudo(n)
    real(c_double), intent(out) :: g(n), dfe(n), tau(n), lgsd(n), mue(n)
    double precision :: newtex, dglpfacthick, demang, interper, pfunc_raw, dlgfacthick
    integer :: get_index
    integer :: i, kk
    double precision :: tausd, cosfac, mus, ptf, em, h, zc
    real :: gsd
    kerr_metric = krz_KerrMetric_init(1.0d0, sv_model%a)
    h = sv_model%h(1)
    zc = sv_model%zcos
    do i = 1, n
        g(i) = dlgfacthick(sv_model%a, sv_model%muobs, alpha(i), re(i), sv_mudisk)
        kk = get_index(rlp(:, 1), ndelta, re(i), sv_risco, npts(1))
        if (nonrel .ne. 0) then
            ! as in sum_multiple_lampposts: sin0, sindisk are uninitialised
            ! (zero) there, so the azimuthal term vanishes
            tau(i) = sqrt(re(i)**2 + (h - sv_model%honr * re(i))**2)          &
                - re(i) * (sv_model%muobs * sv_mudisk) + h * sv_model%muobs
            tau(i) = (1.d0 + zc) * tau(i)
        else
            tausd = interper(rlp(:, 1), tlp(:, 1), ndelta, re(i), kk)
            tau(i) = (1.d0 + zc) * (tausd + taudo(i) - tauso(1))
        end if
        cosfac = interper(rlp(:, 1), dcosdr(:, 1), ndelta, re(i), kk)
        mus = interper(rlp(:, 1), cosd(:, 1), ndelta, re(i), kk)
        if (kk .eq. npts(1)) then
            cosfac = newtex(rlp(:, 1), dcosdr(:, 1), ndelta, re(i), h, sv_model%honr, kk)
            mus = newtex(rlp(:, 1), cosd(:, 1), ndelta, re(i), h, sv_model%honr, kk)
        end if
        ptf = pnorm * pfunc_raw(-mus, sv_model%b1, sv_model%b2, sv_model%qboost)
        gsd = real(dglpfacthick(re(i), sv_model%a, h, sv_mudisk))
        em = determine_emissivity(re(i), sv_model%a, sv_model%Gamma, cosfac, ptf, gsd)
        dfe(i) = em * (g(i) / (1.d0 + zc))**(2. + sv_model%Gamma)
        lgsd(i) = dble(real(log(gsd)))
        mue(i) = demang(sv_model%a, sv_model%muobs, re(i), alpha(i), beta(i))
    end do
end subroutine rt_sample_quantities

! Flat-space (straight line) mapping used for the outer disc.
subroutine rt_flat_batch(n, alpha, beta, re) bind(C, name="rt_flat_batch")
    use iso_c_binding
    use rt_state
    implicit none
    integer(c_int), intent(in) :: n
    real(c_double), intent(in) :: alpha(n), beta(n)
    real(c_double), intent(out) :: re(n)
    double precision :: phie
    integer :: i
    do i = 1, n
        call drandphithick(alpha(i), beta(i), sv_model%muobs, sv_mudisk, re(i), phie)
    end do
end subroutine rt_flat_batch

subroutine rt_state_info(ne, nf, xe, a, muobs, rin, rout, h, mudisk, risco, mueff, zcos, gam) &
        bind(C, name="rt_state_info")
    use iso_c_binding
    use rt_state
    implicit none
    integer(c_int), intent(out) :: ne, nf, xe
    real(c_double), intent(out) :: a, muobs, rin, rout, h, mudisk, risco, mueff, zcos, gam
    ne = sv_ne; nf = sv_nf; xe = sv_xe
    a = sv_model%a; muobs = sv_model%muobs; rin = sv_model%rin; rout = sv_model%rout
    h = sv_model%h(1); mudisk = sv_mudisk; risco = sv_risco; mueff = sv_mueff
    zcos = sv_model%zcos; gam = sv_model%Gamma
end subroutine rt_state_info

subroutine rt_get_freqs(nf, fi) bind(C, name="rt_get_freqs")
    use iso_c_binding
    use rt_state
    implicit none
    integer(c_int), intent(in) :: nf
    real(c_double), intent(out) :: fi(nf)
    fi = sv_fi(1:nf)
end subroutine rt_get_freqs

! Copy out the last kernels computed by rtrans (before any override).
subroutine rt_get_kernel(ne, nf, xe, w0, w1, w2, w3) bind(C, name="rt_get_kernel")
    use iso_c_binding
    use rt_state
    implicit none
    integer(c_int), intent(in) :: ne, nf, xe
    complex(c_float_complex), intent(out) :: w0(ne, nf, xe), w1(ne, nf, xe), w2(ne, nf, xe), w3(ne, nf, xe)
    w0 = sv_w0; w1 = sv_w1; w2 = sv_w2; w3 = sv_w3
end subroutine rt_get_kernel

subroutine rt_set_kernel(ne, nf, xe, w0, w1) bind(C, name="rt_set_kernel")
    use iso_c_binding
    use rt_state
    implicit none
    integer(c_int), intent(in) :: ne, nf, xe
    complex(c_float_complex), intent(in) :: w0(ne, nf, xe), w1(ne, nf, xe)
    if (allocated(ov_w0)) deallocate(ov_w0, ov_w1)
    allocate(ov_w0(ne, nf, xe), ov_w1(ne, nf, xe))
    ov_w0 = w0; ov_w1 = w1
    ov_active = .true.
end subroutine rt_set_kernel

subroutine rt_clear_kernel() bind(C, name="rt_clear_kernel")
    use rt_state
    implicit none
    ov_active = .false.
end subroutine rt_clear_kernel
