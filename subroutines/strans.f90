! vim: cc=80 wrap tw=80
module m_rtrans
    use rtconstants, only: wp
    use common_types
    implicit none

    type :: t_rtrans_args
        type(t_config), pointer :: conf => null()
        type(t_model_arguments), pointer :: model => null()
        type(t_arrays), pointer :: arrays => null()

        ! Observer reflection fractions
        real(wp), pointer :: frobs(:) => null()

        real(wp), pointer :: dFe(:) => null()
        real(wp), pointer :: fi(:) => null()

        ! Number of energy bins. This is passed to `strans` from the original
        ! callsite (e.g. XSPEC or similar)
        integer :: ne

        ! Computed values: the ISCO, the cosine angle on the disc, the radial
        ! bin width, the effective cosine angle
        real(wp) :: r_isco, mudisk, dlogr, mueff

        ! These are set as part of initialising the impulse response matrix
        real(wp) :: dlogt = 0.0_wp, dg = 0.0_wp
    end type t_rtrans_args

contains

    ! Bind the arguments into `t_rtrans_args`, and compute some of the other
    ! values. This initialises and sets up all fields of the structure.
    subroutine bind_arguments(args, config, model_args, arrays, frobs, dFe,    &
         fi, ne)
        type(t_rtrans_args), intent(inout) :: args
        type(t_config), target, intent(in) :: config
        type(t_model_arguments), target, intent(in) :: model_args
        type(t_arrays), target, intent(in) :: arrays
        real(wp), target, intent(in) :: frobs(model_args%nlp)
        real(wp), target, intent(in) :: dFe(model_args%nlp),                   &
             fi(config%nf)
         integer, intent(in) :: ne
        ! functions
        real(wp) :: disco
        args%arrays => arrays
        args%model => model_args
        args%conf => config
        args%frobs => frobs
        args%dFe => dFe
        args%fi => fi
        args%ne = ne

        ! computed values
        args%r_isco = disco(args%model%a)
        args%mudisk = args%model%honr / sqrt(args%model%honr**2 + 1.0_wp)
        args%dlogr = log10(args%conf%rnmax / args%model%rin) /                 &
            real(args%conf%xe - 1, wp)
        args%mueff = max(args%model%muobs, 0.3_wp)
    end subroutine bind_arguments

    ! Zeros all the output arrays
    subroutine outputs_zero_arrays(args)
        type(t_rtrans_args), intent(inout) :: args
        args%arrays%ker_W0 = 0.0_wp
        args%arrays%ker_W1 = 0.0_wp
        args%arrays%ker_W2 = 0.0_wp
        args%arrays%ker_W3 = 0.0_wp
        args%frobs = 0.0_wp
        args%dFe = 0.0_wp
    end subroutine outputs_zero_arrays

end module m_rtrans

!-----------------------------------------------------------------------
subroutine rtrans(config, model_args, arrays, dset, d, ne, frobs, frrel)
    ! Code to calculate the transfer function for an accretion disk.
    ! This code first does full GR ray tracing for a camera with impact parameters < bmax
    ! It then also does straight line ray tracing for impact parameters >bmax
    ! It adds both up to produce a transfer function for a disk extending from rin to rout
    !
    ! This routine populates the transfer functions in `arrays`
    !
    ! Non-standard arguments:
    !
    ! d: Distance of the source
    ! dset: dset=1 means calculate ionization from distance, dset=0 means ignore
    ! distance
    ! ne: Number of energy bins
    ! frobs: Observer's reflection fraction
    ! frrel: Reflection fraction defined by relxilllp.
    use dyn_gr, only: rlp, cosd, dcosdr, ndelta, npts, rlp, tlp, ndelta
    use radial_grids, only: dfer_arr, pnorm
    use gr_continuum, only: gso, lens, tauso, gso, cosdelta_obs
    use saved_variables
    use impulseresponse, only: response_zero_edges, response_allocate
    use common_types, only: t_config, t_model_arguments, t_arrays
    use m_rtrans
    use raytracing, only: trace_disk_observer, getdcos, getlens
    use rtconstants, only: pi, wp
    implicit none

    type(t_config), intent(inout) :: config
    type(t_model_arguments), intent(in) :: model_args
    real(wp), intent(in) :: d
    integer, intent(in) :: ne, dset

    type(t_arrays), intent(inout) :: arrays
    real(wp), intent(inout) :: frobs(model_args%nlp),                          &
         frrel(model_args%nlp)

    type(t_rtrans_args) :: args

    integer m
    real(wp) cosdout(model_args%nlp)
    real(wp) domega(config%nro)

    real(wp) sin0
    integer fbin
    real(wp) rfunc, scal, velocity(3), sysfref
    real(wp) rnmin, rn(config%nro)
    real(wp) fi(config%nf), dgsofac, sindisk
    real(wp) :: dFe(model_args%nlp)
    real(wp) pnormer, pfunc_raw, ang_fac
    real(wp) rnn(config%nro), domegan(config%nro)
    logical dotrace

    ! Setup the output arrays
    call bind_arguments(args, config, model_args, arrays, frobs, dFe, fi, ne)
    ! Zero the outputs
    dfer_arr = 0.0_wp
    call outputs_zero_arrays(args)

    ! Settings/initialization
    scal = 1.0_wp
    velocity = 0.0_wp
    sindisk = sqrt(1.0_wp - args%mudisk**2)

    if (args%conf%calculate_impulse_response) then
        ! This also sets up the time and energy axes
        ! and assigns the dlogt and dg global variables
        call response_allocate(args%ne, args%conf%nt, args%dlogt, args%dg)
    end if

    !set up saving the impulse response function if user desires
    !note: the ideal parameters to plot the transfer function are nro~=7000,nphi~=7000,nt~=2e9,nex~=2e10

    !get the GR ray-tracing CONTINUUM parameters which are stored in the module gr_continuum
    if (args%model%nlp .eq. 1) then
       gso(1) = dgsofac(args%model%a, args%model%h(1))
       call getlens(args%model%metric, args%model%h(1), args%model%muobs,      &
            lens(1), tauso(1), cosdelta_obs(1))
       if (tauso(1) .ne. tauso(1)) stop "tauso is NaN"
    else
       !here the observed cutoffs are set from the temperature in the source frame
       do m = 1, args%model%nlp
          gso(m) = dgsofac(args%model%a, args%model%h(m))
          call getlens(args%model%metric, args%model%h(m), args%model%muobs,   &
               lens(m), tauso(m), cosdelta_obs(m))
          if (tauso(m) .ne. tauso(m)) stop "tauso is NaN"
       enddo
    endif

    ! Set up observer's camera ( alpha = rn sin(phin), beta = mueff rn cos(phin) )
    ! to do full GR ray tracing with
    rnmin = rfunc(args%model%a, args%model%muobs)
    !Grid to do in full GR
    call getrgrid(rnmin, args%conf%rnmax, args%mueff, args%conf%nro,           &
         args%conf%nphi, rn, domega)
    !Grid for Newtonian approximation
    call getrgrid(args%conf%rnmax, args%model%rout, args%mueff,                &
         args%conf%nron, args%conf%nphin, rnn, domegan)

    ! Trace rays in full GR for the small camera (ie with relativistic effects) from the osberver to the disk,
    ! which is why it doesnt depend on h
        dotrace = .false.
        if (abs(spinsav - args%model%a) .gt. tiny(args%model%a)) then
            dotrace = .true.
        endif
        if (abs(musav - args%model%muobs) .gt. tiny(args%model%muobs)) then
            dotrace = .true.
        endif
        if (abs(routsav - args%model%rout) .gt. tiny(args%model%rout)) then
            dotrace = .true.
        endif
        if (abs(mudsav-args%mudisk) .gt. tiny(args%mudisk)) dotrace = .true.
        if (dotrace) then
            call trace_disk_observer(args%model%metric, args%conf%nro,         &
                args%conf%nphi, rn, args%mueff, args%model%muobs, args%r_isco, &
                args%model%rout, d)
            spinsav = args%model%a
            musav = args%model%muobs
            routsav = args%model%rout
            mudsav = args%mudisk
       end if

    ! Set frequency array
    do fbin = 1, args%conf%nf
        args%fi(fbin) = args%conf%flo * (args%conf%fhi /                       &
            args%conf%flo)**((real(fbin, wp) - 0.5_wp) / real(args%conf%nf, wp))
    end do
    if (args%conf%fhi .lt. tiny(args%conf%fhi)) args%fi(1) = 0.0_wp

    !initialize radius grid, angles, and transfer functions
    sin0 = sqrt(1.0_wp-args%model%muobs**2)

    ! Calculate dcos/dr and time lags vs r for the lamppost model
    call getdcos(args%model%metric, args%model%h, ndelta, args%model%nlp,      &
         args%model%rout, npts, rlp, dcosdr, tlp, cosd, cosdout)

    ! set continuum normalisations depending on model flavour
    if (dset .eq. 0)then
        pnorm = 1.0_wp / (4.0_wp * pi)
    else
        pnorm = pnormer(args%model%b1, args%model%b2, args%model%qboost)
    end if

    ! the only arguments that change here are .false., nro, nphi, rn, domega
    ! the first call is for the relativistic version
    call sum_impulse_components(.false., args%conf%nro, args%conf%nphi, rn,    &
         domega, args)
    ! then for the non-relativistic flat-space version
    call sum_impulse_components(.true., args%conf%nron, args%conf%nphin,       &
         rnn, domegan, args)

    do m = 1, args%model%nlp
        ! Calculate 4pi p(theta0,phi0) = ang_fac
        ang_fac = 4.0_wp * pi * pnorm * pfunc_raw(-cosdelta_obs(m),            &
            args%model%b1, args%model%b2, args%model%qboost)
        ! Adjust the lensing factor (easiest way to keep track)
        lens(m) = lens(m) * ang_fac
        ! Calculate the relxill reflection fraction for one columncosdout
        frrel(m) = sysfref(args%model%rin, rlp(:, m), cosd(:, m), ndelta,      &
            cosdout(m))
        !Finish calculation of observer's reflection fraction
        args%frobs(m) = args%frobs(m) / dgsofac(args%model%a,                  &
            args%model%h(m)) / lens(m)
    end do

    if (args%conf%calculate_impulse_response) then
        ! Deal with edge effects
        call response_zero_edges()
    end if

    return
end subroutine rtrans
!-----------------------------------------------------------------------

! Clamps an integer value between some high and low value.
integer function clamp_i(v, low, high) result(o)
    integer, intent(in) :: v, low, high
    o = max(low, v)
    o = min(o, high)
end function clamp_i

! This is an attempt to cleanup the strans function to put common code into a
! common subroutine so that there are fewer edits needed to add args behaviours
subroutine sum_impulse_components(non_relativistic, r_length, phi_length,      &
     r_grid, domega, args)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use rtconstants, only: pi, wp
    use emissivities
    use m_rtrans
    implicit none
    logical, intent(in) :: non_relativistic

    integer, intent(in) :: r_length, phi_length
    ! lamppost heights
    real(wp), intent(in) :: r_grid(r_length)
    real(wp), intent(in) :: domega(r_length)

    ! functions
    real(wp) :: dlgfacthick
    integer :: clamp_i

    type(t_rtrans_args), intent(inout) :: args

    real(wp) :: re, alpha, beta, phie, phin
    ! photon time from/to
    ! tauso is in `gr_continuum`
    real(wp) :: taudo, g
    integer :: ri, i, j, gbin, rbin
    ! Set to true once the disc has been seen by the ray-tracing techniques.
    ! This is to avoid a bug where a disc with the outer radius truncated below
    ! the maximum of r_grid would cause the ray-tracing to terminate
    ! prematurely.
    logical :: disc_seen
    ! Used in the loop to check if all over the points along a particular r_grid
    ! radius on the image plane hit the accretion disc.
    logical :: at_least_one_hit

    ! loop over all photon directions (l), disk radii (i), disk azimuth (j), and
    ! calculate the contribution to the
    ! transfer function/convolution kernel in energy (gbin), frequency (fbin),
    ! emission angle (mubin), disk radial
    ! bin (rbin) from the m-th/nl-th lamp post

    ! TODO: for ring-like corona, pre-load the correct time-dependent emissivity
    ! profile here, before the loop over observer coordinates

    disc_seen = .false.
    do ri = 1, r_length
        ! Loop in reverse order over the radial grid, starting at the
        ! outermost radius.
        i = r_length - (ri - 1)
        at_least_one_hit = .false.
        do j = 1, phi_length
            phin = (j-0.5_wp) * 2.0_wp * pi / real(phi_length, wp)
            alpha = r_grid(i) * sin(phin)
            beta = -r_grid(i) * cos(phin) * args%mueff

            ! If the ray hits the disk, calculate flux and time lag
            if (non_relativistic) then
                call drandphithick(alpha, beta, args%model%muobs,              &
                     args%mudisk, re, phie)
            else
                if (pem1(j, i) .le. 0.0_wp) then
                    ! Did not intersect with the accretion disc.
                    cycle
                endif
                re = re1(j, i)
                taudo = taudo1(j, i)
            endif

            if (re .lt. args%model%rin .or. re .gt. args%model%rout) then
                ! Not in the disc domain, skip this photon and move on to the
                ! next j.
                cycle
            endif

            ! Mark this photon as hit within disc boundaries.
            at_least_one_hit = .true.

            ! disc to observer energy shift
            g = dlgfacthick(args%model%a, args%model%muobs, alpha, re,         &
                args%mudisk)

            ! Work out energy bin
            gbin = clamp_i(ceiling(log10(g / (1.0_wp + args%model%zcos)) /     &
                args%conf%dloge) + args%ne / 2, 1, args%ne)

            ! Work out radial bin
            rbin = clamp_i(ceiling(log10(re / args%model%rin) /                &
                args%dlogr), 1, args%conf%xe)

            if (args%conf%ring_like) then
                call sum_ringlike_corona(i, non_relativistic, r_length,        &
                     re, alpha, beta, taudo, g, domega, gbin, rbin, args)
            else
                call sum_multiple_lampposts(i, non_relativistic, r_length,     &
                     re, alpha, beta, taudo, g, domega, gbin, rbin, args)
            endif
        end do

        if (disc_seen) then
            if (.not. at_least_one_hit) then
                ! Break out of the loop early to avoid unnecessary ray-tracing.
                exit
            end if
        else
            if (at_least_one_hit) then
                ! The disc has now been seen.
                disc_seen = .true.
            end if
        end if
    end do
end subroutine sum_impulse_components

subroutine sum_ringlike_corona(i, non_relativistic, r_length,                  &
     re, alpha, beta, taudo, g, domega, gbin, rbin, args)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use rtconstants, only: pi, wp
    use emissivities
    use impulseresponse, only: time_axis, response
    use m_rtrans
    implicit none

    logical, intent(in) :: non_relativistic
    integer, intent(in) :: i, r_length, gbin, rbin
    real(wp), intent(in) :: domega(r_length)
    real(wp), intent(in) :: alpha, beta, taudo, g, re

    type(t_rtrans_args), intent(inout) :: args

    ! functions
    real(wp) :: dareafac, demang, dglpfacthick
    integer :: clamp_i

    real(wp) :: sin0, mue, cosfac, sindisk, phie
    integer :: gbin_resp, tbin, mubin, fbin
    real(wp) :: gsd, normfac
    complex(wp) :: cexp

    real(wp) :: tausd, tau, emissivity

    ! this is a fixed number for now, representing the number of bins in azimuth
    integer, parameter :: r_nphi = 50
    real(wp), parameter :: dphi = 2 * pi / real(r_nphi, wp)
    ! index counting which phi bin we are currently considering
    integer :: phi_i
    real(wp) :: phi

    if (args%model%nlp .ne. 1) then
        print *, "panic: expected only one corona for ring-like corona"
        error stop 1
    end if

    ! Add to reflection fraction
    args%frobs(1) = args%frobs(1) + 2.0_wp * g**3 * gsd * cosfac /             &
        dareafac(re, args%model%a) * domega(i)

    ! Calculate flux from pixel
    gsd = dglpfacthick(re, args%model%a, args%model%h(1), args%mudisk)

    normfac = (g/(1.0_wp+args%model%zcos))**(2.0_wp+args%model%Gamma)*domega(i)

    ! the observed energy bin for the response matrix
    gbin_resp = clamp_i(ceiling(g/args%dg), 1, args%ne)

    ! calculate emission angle and work out which mue bin to add to
    mue = demang(args%model%a, args%model%muobs, re, alpha, beta)
    mubin = ceiling(mue * real(args%conf%me, wp))

    ! loop over all azmithal bins
    do phi_i = 1, r_nphi
        phi = phi_i * dphi

        ! the source to disc time of the current azimuthal bin
        call get_emissivity_time(re, phi, emissivity, tausd)
        ! normalise
        emissivity = emissivity / real(r_nphi, wp)

        if (non_relativistic) then
            ! TODO: for the non-relativistic case, can likely also consider the
            ! corona to be a lamppost, since the spread of time values will be
            ! small, likely of order the size of the ring
            ! - this should be checked, else a full non-relatvistic version that
            !   loops over each azimuth used
            tau = sqrt(re**2 + (args%model%h(1) - args%model%honr *            &
                re)**2) - re * (sin0 * sindisk * cos(phie) + args%model%muobs  &
                * args%mudisk) + args%model%h(1) * args%model%muobs
            tau = (1.0_wp+args%model%zcos)*tau
        else
            tau = (1.0_wp + args%model%zcos) * (tausd + taudo - tauso(1))
        endif

        args%dFe(1) = (emissivity * (g / (1.0_wp + args%model%zcos))**         &
            (2.0_wp + args%model%Gamma) * domega(i))

        dfer_arr(rbin) = dfer_arr(rbin) + args%dFe(1)

        ! Add to the transfer function integral
        do fbin = 1, args%conf%nf
            cexp = cmplx(cos(2.0_wp * pi * tau * args%fi(fbin)),               &
                sin(2.0_wp * pi * tau * args%fi(fbin)), kind=wp)

            args%arrays%ker_W0(1, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W0(1, gbin, fbin, mubin, rbin) +               &
                args%dFe(1) * cexp

            ! TODO: the below are all particular to the lamppost corona, and do
            ! not apply to the ring-like corona currently

            args%arrays%ker_W1(1, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W1(1, gbin, fbin, mubin, rbin) +               &
                log(gsd) * args%dFe(1) * cexp

            args%arrays%ker_W2(1, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W2(1, gbin, fbin, mubin, rbin) + emissivity *  &
                normfac * cexp

            args%arrays%ker_W3(1, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W3(1, gbin, fbin, mubin, rbin) + emissivity *  &
                normfac * cexp
        end do

        if (args%conf%calculate_impulse_response) then
            ! find the appropriate energy and time bins
            tbin = clamp_i(ceiling(log10(tau / time_axis(0)) /                 &
                args%dlogt), 1, args%conf%nt)
            ! kernel of the impulse response function
            response(gbin_resp, tbin) = response(gbin_resp, tbin) + args%dFe(1)
        end if
    end do
end subroutine sum_ringlike_corona

subroutine sum_multiple_lampposts(i, non_relativistic, r_length,               &
     re, alpha, beta, taudo, g, domega, gbin, rbin, args)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use rtconstants, only: pi, wp
    use emissivities
    use impulseresponse, only: time_axis, response
    use m_rtrans
    implicit none

    logical, intent(in) :: non_relativistic
    integer, intent(in) :: i, r_length, gbin, rbin
    real(wp), intent(in) :: domega(r_length)
    real(wp), intent(in) :: alpha, beta, taudo, g, re

    type(t_rtrans_args), intent(inout) :: args

    integer :: tbin, fbin

    ! functions
    real(wp) :: newtex, dglpfacthick, demang, interper, dareafac,              &
         pfunc_raw
    integer :: get_index, clamp_i

    real(wp) :: sin0, sindisk
    real(wp) :: phie
    real(wp) :: cosfac, mus, ptf
    real(wp) :: kfac, normfac, emisfac
    real(wp) :: thetafac(args%model%nlp), gsd(args%model%nlp)
    ! photon time from/to
    ! tauso is in `gr_continuum`
    real(wp) :: tausd, mue
    real(wp) :: tau(args%model%nlp), emissivity(args%model%nlp)
    integer :: m, mubin, kk, gbin_resp
    complex(wp) :: cexp

    do m = 1, args%model%nlp

        kk = get_index(rlp(:, m), ndelta, re, args%r_isco, npts(m))

        ! Time lag between direct and reflected photons
        if (non_relativistic) then
            tau(m) = sqrt(re**2 + (args%model%h(m) - args%model%honr *         &
                re)**2) - re * (sin0 * sindisk * cos(phie) + args%model%muobs  &
                * args%mudisk) + args%model%h(1) * args%model%muobs
            tau(m) = (1.0_wp+args%model%zcos)*tau(m)
        else
            ! Interpolate (or extrapolate) the time function
            tausd = interper(rlp(:, m), tlp(:, m), ndelta, re, kk)
            tau(m) = (1.0_wp+args%model%zcos)*(tausd+taudo-tauso(1))
        endif

        ! Interpolate |dcos\delta/dr| function
        cosfac = interper(rlp(:, m), dcosdr(:, m), ndelta, re, kk)
        mus = interper(rlp(:, m), cosd(:, m), ndelta, re, kk)

        if (kk .eq. npts(m)) then
            cosfac = newtex(rlp(:, m), dcosdr(:, m), ndelta, re,               &
                args%model%h(m), args%model%honr, kk)
            mus = newtex(rlp(:, m), cosd(:, m), ndelta, re,                    &
                args%model%h(m), args%model%honr, kk)
        end if

        ! Calculate angular emissivity
        ptf = pnorm * pfunc_raw(-mus, args%model%b1, args%model%b2,            &
            args%model%qboost)

        ! Calculate flux from pixel
        gsd(m) = dglpfacthick(re, args%model%a, args%model%h(m),               &
            args%mudisk)

        ! TODO: write into emissivity(:, m) where the emissivity
        ! now holds the times as well from the time-dependent emissivity
        ! functions
        emissivity(m) = determine_emissivity(re, args%model%a,                 &
            args%model%Gamma, cosfac, ptf, gsd(m))

        ! calculate extra factors that go into the transfer functions
        ! for double lps
        if (args%model%nlp .gt. 1) then
            thetafac(m) = emissivity(m) * gso(m)**(args%model%Gamma - 2.0_wp)  &
                * gsd(m)**(2.0_wp - args%model%Gamma)
        else
            ! single lamp post case, double check this later
            thetafac(m) = 1.0_wp
        endif
        ! Add to reflection fraction
        args%frobs(m) = args%frobs(m) + 2.0_wp * g**3 * gsd(m) * cosfac /      &
            dareafac(re, args%model%a) * domega(i)
    enddo

    ! This loop combines information from the two lampposts, and so
    ! cannot be merged with the above, as it needs to known information
    ! simultaneously

    ! TODO: extend to add a loop over the emissivity time, except that
    ! it is not compatible with the existing use, as the dimensions of
    ! the emissivity would be different
    do m = 1, args%model%nlp
        args%dFe(m) = emissivity(m) * (g / (1.0_wp + args%model%zcos))**       &
            (2.0_wp + args%model%Gamma) * domega(i)

        ! Add to the radial dependence of the transfer function TBD MAKE
        ! SURE THIS IS RIGHT

        dfer_arr(rbin) = dfer_arr(rbin) + args%dFe(m)
        ! Calculate emission angle and work out which mue bin to add to
        mue = demang(args%model%a, args%model%muobs, re, alpha, beta)
        mubin = ceiling(mue * real(args%conf%me, wp))
        !calculate the extra factors for w2/3
        if (args%model%nlp .gt. 1) then
            emisfac = (emissivity(1) + args%model%eta_0 * emissivity(2)) /&
                (1.0_wp + args%model%eta_0)

            kfac = (emissivity(1) + args%model%eta_0 * emissivity(2)) /        &
                (thetafac(1) + args%model%eta_0 * thetafac(2))
        else
            emisfac = emissivity(1)
            kfac = emissivity(1)
            ! single lamp post case, double check this later
        endif

        normfac = (g / (1.0_wp + args%model%zcos))**(2.0_wp +                  &
            args%model%Gamma) * domega(i)
        ! Add to the transfer function integral
        do fbin = 1, args%conf%nf
            cexp = cmplx(cos(2.0_wp * pi * tau(m) * args%fi(fbin)),            &
                sin(2.0_wp * pi * tau(m) * args%fi(fbin)), kind=wp)

            args%arrays%ker_W0(m, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W0(m, gbin, fbin, mubin, rbin) +               &
                args%dFe(m) * cexp

            args%arrays%ker_W1(m, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W1(m, gbin, fbin, mubin, rbin) +               &
                log(gsd(m)) * args%dFe(m) * cexp

            ! tbd redo these transfer functions
            args%arrays%ker_W2(m, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W2(m, gbin, fbin, mubin, rbin) + emisfac *     &
                normfac * cexp

            args%arrays%ker_W3(m, gbin, fbin, mubin, rbin) =                   &
                args%arrays%ker_W3(m, gbin, fbin, mubin, rbin) + kfac *        &
                thetafac(m) * normfac * cexp
        end do
    end do

    if (args%conf%calculate_impulse_response) then
        do m = 1, args%model%nlp
            ! find the appropriate energy and time bins
            gbin_resp = clamp_i(ceiling(g/args%dg), 1, args%ne)
            tbin = clamp_i(ceiling(log10(tau(m) / time_axis(0)) /              &
                args%dlogt), 1, args%conf%nt)
            ! kernel of the impulse response function
            response(gbin_resp, tbin) = response(gbin_resp, tbin) + args%dFe(m)
        end do
    end if
end subroutine sum_multiple_lampposts

!-----------------------------------------------------------------------
function newtex(rlp, dcosdr, ndelta, re, h, honr, kk)
! Extrapolates using Newtonian value
  use rtconstants, only: wp
  implicit none
  integer ndelta, kk
  real(wp) newtex, rlp(ndelta), dcosdr(ndelta), re, h, honr
  newtex = dcosdr(kk) * ((h-honr*rlp(kk))**2 + rlp(kk)**2)**1.5_wp / rlp(kk)
  newtex = newtex * re / ((h-honr*re)**2 + re**2)**1.5_wp
  return
end function newtex
!-----------------------------------------------------------------------
