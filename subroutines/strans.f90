! vim: cc=80 wrap tw=80
include 'subroutines/emissivity.f90'
include 'subroutines/impulseresponse.f90'

module m_rtrans
    implicit none

    ! The arguments passed to various parts of the transfer function / impulse
    ! response summation
    type :: t_impulse_args
        ! nf: number of frequency bins
        ! ne: number of energy bins
        ! me: number of mue bins
        ! xe: number of logr bins
        ! nlp: number of lampposts
        integer :: nf, ne, me, xe, nlp

        ! dloge: energy delta energy
        real :: dloge

        ! use ring-like corona emissivity
        logical :: ring_like = .false.

        ! Gamma: photon index
        ! mu0: cosine of the observer's inclination angle
        double precision :: Gamma, mu0

        ! rin / rout: inner and outer disc radius
        ! rmin: minimum possible disc radius (ISCO)
        double precision :: rin, rout, rmin

        ! zcos: cosmological redshift
        ! honr: scale height
        double precision :: mudisk, zcos, honr

        ! b1: linear coefficient of angular emissivity function
        ! b2: quadratic coefficient of angular emissivity function
        ! qboost: asymmetry parameter of angular emissivity function
        double precision :: b1, b2, qboost

        ! rnmax: max radius for which GR ray-tracing is used
        ! spin: spin of the black hole
        double precision :: rnmax, spin, dlogr, eta_0

        double precision :: mueff

        ! these are set as part of initialising the impulse response matrix
        double precision :: dlogt, dg

        ! h(nlp): lamppost heights
        double precision, pointer :: h(:) => null()
        double precision, pointer :: fi(:) => null()

        ! number of time bins
        integer :: nt = 2**9
    end type t_impulse_args

    type :: t_outputs
        ! w0(nlp,ne,nf,me,xe): linear transfer function
        ! w1, w2: power law index pivoting
        ! w3: ionisation pivoting
        ! all of these share the same dimensions
        complex, pointer :: w0(:, :, :, :, :) => null(), w1(:, :, :, :,        &
             :) => null(), w2(:, :, :, :, :) => null(), w3(:, :, :, :,         &
             :) => null()
        ! frobs(nlp): observer reflection fractions
        double precision, pointer :: frobs(:) => null()
        double precision, pointer :: dFe(:) => null()
    end type t_outputs
contains

    ! Bind views to various input arrays to the `t_impulse_args` derived type
    subroutine impulse_bind_views(self, h, fi)
        type(t_impulse_args), intent(inout) :: self
        double precision, target, intent(in) :: h(self%nlp)
        double precision, target, intent(in) :: fi(self%nf)
        self%h => h
        self%fi => fi
    end subroutine impulse_bind_views

    ! Bind views to various output arrays to the `t_outputs` derived type
    subroutine outputs_bind_views(self, w0, w1, w2, w3, frobs, dFe)
        type(t_outputs), intent(inout) :: self
        complex, target, intent(in) :: w0(:, :, :, :, :), w1(:, :, :, :, :),   &
             w2(:, :, :, :, :), w3(:, :, :, :, :)
        double precision, target, intent(in) :: frobs(:), dFe(:)
        self%w0 => w0
        self%w1 => w1
        self%w2 => w2
        self%w3 => w3
        self%frobs => frobs
        self%dFe => dFe
    end subroutine outputs_bind_views

    ! Zeros all the output arrays
    subroutine outputs_zero_arrays(self)
        type(t_outputs), intent(inout) :: self
        self%w0 = 0.
        self%w1 = 0.
        self%w2 = 0.
        self%w3 = 0.
        self%frobs = 0.0
        self%dFe = 0.0
    end subroutine outputs_zero_arrays

end module m_rtrans

!-----------------------------------------------------------------------
subroutine rtrans(verbose, dset, nlp, spin, h, mu0, Gamma, rin, rout, honr,    &
     d, rnmax, zcos, b1, b2, qboost, eta_0, fcons, nro, nphi, ne, dloge, nf,   &
     fhi, flo, me, xe, ker_W0, ker_W1, ker_W2, ker_W3, frobs, frrel)
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
    use blcoordinate, only: pi
    use radial_grids, only: dfer_arr, pnorm
    use gr_continuum, only: gso, lens, tauso, gso, cosdelta_obs
    use saved_variables
    use impulseresponse
    use m_rtrans
    implicit none
    integer nro, nphi, ne, nf, me, xe, dset, nlp
    double precision spin, h(2), mu0, Gamma, rin, rout, zcos, fhi, flo, honr
    double precision b1, b2, qboost
    double precision fcons, cosdout(nlp)
    real dloge

    integer m
    double precision domega(nro), d, dFe(nlp)
    !double precision rlp_column(ndelta),dcosdr_column(ndelta),tlp_column(ndelta),cosd_column(ndelta)
    double precision cos0, sin0
    integer fbin
    double precision disco, rfunc, scal, velocity(3), sysfref
    double precision rnmax, rnmin, rn(nro), mueff
    double precision fi(nf), dgsofac, sindisk, frobs(nlp), frrel(nlp)
    double precision pnormer, pfunc_raw, ang_fac
    integer nron, nphin, verbose
    double precision rnn(nro), domegan(nro)
    double precision eta_0
    logical dotrace

    !new stuff - move back above once it's implemented properly
    complex :: ker_W0(nlp, ne, nf, me, xe), ker_W1(nlp, ne, nf, me, xe),       &
         ker_W2(nlp, ne, nf, me, xe), ker_W3(nlp, ne, nf, me, xe)

    ! functions
    integer :: get_env_int

    type(t_impulse_args) :: args
    type(t_outputs) :: ret

    ! Setup the common arguments structure
    args = t_impulse_args(nf = nf, ne = ne, me = me, xe = xe, nlp = nlp,       &
        Gamma = Gamma, rmin = disco(spin), mu0 = mu0, rin = rin,               &
        rout = rout, mudisk = honr / sqrt(honr**2 + 1.d0), zcos = zcos,        &
        honr = honr, qboost = qboost, b1 = b1, b2 = b2, rnmax = rnmax,         &
        spin = spin, eta_0 = eta_0, mueff = mueff, dloge = dloge,              &
        ! TODO: these are uninitialized
        ring_like = (1 .eq. get_env_int("REV_RING", 0)), dlogr = 0.0,          &
            dlogt = 0.0, dg = 0.0)
    call impulse_bind_views(args, h, fi)

    ! Setup the output arrays
    ret = t_outputs()
    call outputs_bind_views(ret, ker_W0, ker_W1, ker_W2, ker_W3, frobs, dFe)

    ! Settings/initialization
    nron = 100
    nphin = 100
    scal = 1.d0
    velocity = 0.d0
    sindisk = sqrt(1.d0 - args%mudisk**2)

    ! Zero the outputs
    dfer_arr = 0.
    call outputs_zero_arrays(ret)

    ! This also sets up the time and energy axes
    ! and assigns the dlogt and dg global variables
    call response_allocate(args%ne, args%nt, args%dlogt, args%dg)

    !set up saving the impulse response function if user desieres
    !note: the ideal parameters to plot the transfer function are nro~=7000,nphi~=7000,nt~=2e9,nex~=2e10

    !get the GR ray-tracing CONTINUUM parameters which are stored in the module gr_continuum
    if (args%nlp .eq. 1) then
       gso(1) = real(dgsofac(args%spin, args%h(1)))
       call getlens(args%spin, args%h(1), args%mu0, lens(1), tauso(1),         &
            cosdelta_obs(1))
       if (tauso(1) .ne. tauso(1)) stop "tauso is NaN"
    else
       !here the observed cutoffs are set from the temperature in the source frame
       do m = 1, args%nlp
          gso(m) = real(dgsofac(args%spin, args%h(m)))
          call getlens(args%spin, args%h(m), args%mu0, lens(m), tauso(m),      &
               cosdelta_obs(m))
          if (tauso(m) .ne. tauso(m)) stop "tauso is NaN"
       enddo
    endif

    ! Set up observer's camera ( alpha = rn sin(phin), beta = mueff rn cos(phin) )
    ! to do full GR ray tracing with
    args%mueff = max(args%mu0, 0.3d0)
    rnmin = rfunc(args%spin, args%mu0)
    !Grid to do in full GR
    call getrgrid(rnmin, args%rnmax, args%mueff, nro, nphi, rn, domega)
    !Grid for Newtonian approximation
    call getrgrid(args%rnmax, args%rout, args%mueff, nron, nphin, rnn, domegan)

    ! Trace rays in full GR for the small camera (ie with relativistic effects) from the osberver to the disk,
    !which is why it doesnt depend on h
        dotrace = .false.
       if( abs(spinsav-args%spin)  .gt. tiny(args%spin)   ) dotrace = .true.
       if( abs(musav-args%mu0)     .gt. tiny(args%mu0)    ) dotrace = .true.
       if( abs(routsav-args%rout)  .gt. tiny(args%rout)   ) dotrace = .true.
       if( abs(mudsav-args%mudisk) .gt. tiny(args%mudisk) ) dotrace = .true.
       if (dotrace)then
            call GRtrace(nro, nphi, rn, args%mueff, args%mu0, args%spin,       &
                 args%rmin, args%rout, args%mudisk, d)
            spinsav = args%spin
            musav = args%mu0
            routsav = args%rout
            mudsav = args%mudisk
       end if

    ! Set frequency array
    do fbin = 1, nf
        args%fi(fbin) = flo * (fhi/flo)**((float(fbin)-0.5d0)/dble(nf))
    end do
    if (fhi .lt. tiny(fhi)) args%fi(1) = 0.0d0

    !initialize radius grid, angles, and transfer functions
    args%dlogr = log10(args%rnmax/args%rin) / real(args%xe-1)
    cos0 = args%mu0
    sin0 = sqrt(1.0-cos0**2)

    ! Calculate dcos/dr and time lags vs r for the lamppost model
    call getdcos(args%spin, args%h, args%mudisk, ndelta, args%nlp,             &
         args%rout, npts, rlp, dcosdr, tlp, cosd, cosdout)

    ! set continuum normalisations depending on model flavour
    if (dset .eq. 0)then
        pnorm = 1.d0 / (4.d0 * pi)
    else
        pnorm = pnormer(args%b1, args%b2, args%qboost)
    end if

    ! the only arguments that change here are .false., nro, nphi, rn, domega
    ! the first call is for the relativistic version
    call sum_impulse_components(args, .false., nro, nphi, rn, domega, ret)
    ! then for the non-relativistic flat-space version
    call sum_impulse_components(args, .true., nron, nphin, rnn, domegan, ret)

    do m = 1, args%nlp
        ! Calculate 4pi p(theta0,phi0) = ang_fac
        ang_fac = 4.d0 * pi * pnorm* pfunc_raw(-cosdelta_obs(m), args%b1,      &
            args%b2, args%qboost)
        ! Adjust the lensing factor (easiest way to keep track)
        lens(m) = lens(m) * ang_fac
        ! Calculate the relxill reflection fraction for one columncosdout
        frrel(m) = sysfref(args%rin, rlp(:, m), cosd(:, m), ndelta, cosdout(m))
        !Finish calculation of observer's reflection fraction
        ret%frobs(m) = ret%frobs(m) / dgsofac(args%spin, args%h(m)) / lens(m)
    end do

    ! Deal with edge effects
    call response_zero_edges()

    return
end subroutine rtrans
!-----------------------------------------------------------------------

integer function clamp_i(v, low, high) result(o)
    integer, intent(in) :: v, low, high
    o = max(low, v)
    o = min(o, high)
end function clamp_i

! This is an attempt to cleanup the strans function to put common code into a
! common subroutine so that there are fewer edits needed to add new behaviours
!
! TODO: this should really be using a derived type to pass arguments around, but
! that's more refactoring for some(time|one) else.
subroutine sum_impulse_components(args, non_relativistic, r_length,            &
     phi_length, r_grid, domega, ret)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use constants
    use emissivities
    use impulseresponse
    use m_rtrans
    implicit none
    logical, intent(in) :: non_relativistic
    type(t_impulse_args), intent(in) :: args

    integer, intent(in) :: r_length, phi_length
    ! lamppost heights
    double precision, intent(in) :: r_grid(r_length)
    double precision, intent(in) :: domega(r_length)

    ! functions
    double precision :: dlgfacthick
    integer :: clamp_i

    type(t_outputs), intent(inout) :: ret

    double precision :: re, alpha, beta, phie, phin
    ! photon time from/to
    ! tauso is in `gr_continuum`
    double precision :: taudo, g
    integer :: i, j, gbin, rbin

    ! loop over all photon directions (l), disk radii (i), disk azimuth (j), and
    ! calculate the contribution to the
    ! transfer function/convolution kernel in energy (gbin), frequency (fbin),
    ! emission angle (mubin), disk radial
    ! bin (rbin) from the m-th/nl-th lamp post

    ! TODO: for ring-like corona, pre-load the correct time-dependent emissivity
    ! profile here, before the loop over observer coordinates

    do i = 1, r_length
        do j = 1, phi_length
            phin = (j-0.5) * 2.d0 * pi / dble(phi_length)
            alpha = r_grid(i) * sin(phin)
            beta = -r_grid(i) * cos(phin) * args%mueff

            ! If the ray hits the disk, calculate flux and time lag
            if (non_relativistic) then
                call drandphithick(alpha, beta, args%mu0, args%mudisk, re, phie)
            else
                if (pem1(j, i) .le. 0.0d0) then
                    cycle
                endif
                re = re1(j, i)
                taudo = taudo1(j, i)
            endif

            if (re .lt. args%rin .and. re .gt. args%rout) then
                ! not in the disc domain, skip this photon
                cycle
            endif

            ! disc to observer energy shift
            g = dlgfacthick(args%spin, args%mu0, alpha, re, args%mudisk)

            ! Work out energy bin
            gbin = clamp_i(ceiling(log10(g                                     &
                /(1.d0+args%zcos)) / args%dloge) + args%ne/2,                  &
                1, args%ne)

            ! Work out radial bin
            rbin = clamp_i(ceiling(log10(re/args%rin) / args%dlogr), 1, args%xe)

            if (args%ring_like) then
                call sum_ringlike_corona(args, i, non_relativistic,            &
                     r_length, phi_length, re, alpha, beta, taudo, g,          &
                     r_grid, domega, gbin, rbin, ret)
            else
                call sum_multiple_lampposts(args, i, non_relativistic,         &
                     r_length, phi_length, re, alpha, beta, taudo, g,          &
                     r_grid, domega, gbin, rbin, ret)
            endif
        end do
    end do
end subroutine sum_impulse_components

subroutine sum_ringlike_corona(args, i, non_relativistic, r_length,            &
     phi_length, re, alpha, beta, taudo, g, r_grid, domega, gbin, rbin,        &
     ret)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use constants
    use emissivities
    use impulseresponse
    use m_rtrans
    implicit none

    logical, intent(in) :: non_relativistic
    type(t_impulse_args), intent(in) :: args
    integer, intent(in) :: i, r_length, phi_length, gbin, rbin
    double precision, intent(in) :: r_grid(r_length)
    double precision, intent(in) :: domega(r_length)
    double precision, intent(in) :: alpha, beta, taudo, g, re

    type(t_outputs), intent(inout) :: ret

    ! functions
    double precision :: dareafac, demang, dglpfacthick
    integer :: clamp_i

    double precision :: sin0, mue, cosfac, sindisk, phie
    integer :: gbin_resp, tbin, mubin, fbin
    real :: gsd, normfac
    complex :: cexp

    double precision :: tausd, tau, emissivity

    ! this is a fixed number for now, representing the number of bins in azimuth
    integer, parameter :: r_nphi = 50
    double precision, parameter :: dphi = 2 * pi / float(r_nphi)
    ! index counting which phi bin we are currently considering
    integer :: phi_i
    double precision :: phi

    if (args%nlp .ne. 1) then
        print *, "panic: expected only one corona for ring-like corona"
        error stop 1
    end if

    ! Add to reflection fraction
    ret%frobs(1) = ret%frobs(1)+ 2.0*g**3*gsd*cosfac/dareafac(re,              &
        args%spin)*domega(i)

    ! Calculate flux from pixel
    gsd = dglpfacthick(re, args%spin, args%h(1), args%mudisk)

    normfac = real((g/(1.d0+args%zcos))**(2.+args%Gamma)*domega(i))

    ! the observed energy bin for the response matrix
    gbin_resp = clamp_i(ceiling(g/args%dg), 1, args%ne)

    ! calculate emission angle and work out which mue bin to add to
    mue = demang(args%spin, args%mu0, re, alpha, beta)
    mubin = ceiling(mue * dble(args%me))

    ! loop over all azmithal bins
    do phi_i = 1, r_nphi
        phi = phi_i * dphi

        ! the source to disc time of the current azimuthal bin
        call get_emissivity_time(re, phi, emissivity, tausd)
        ! normalise
        emissivity = emissivity / float(r_nphi)

        if (non_relativistic) then
            ! TODO: for the non-relativistic case, can likely also consider the
            ! corona to be a lamppost, since the spread of time values will be
            ! small, likely of order the size of the ring
            ! - this should be checked, else a full non-relatvistic version that
            !   loops over each azimuth used
            tau = sqrt(re                                                      &
                **2                                                            &
                +(args%h(1)                                                    &
                -args%honr                                                     &
                *re)**2)                                                       &
                - re                                                           &
                *(sin0*sindisk*cos(phie)                                       &
                +args%mu0*args%mudisk)+ args%h(1)*args%mu0
            tau = (1.d0+args%zcos)*tau
        else
            tau = (1.d0 + args%zcos) * (tausd + taudo - tauso(1))
        endif

        ret%dFe(1) = (emissivity                                               &
            * (g/(1.d0+args%zcos))**(2.+args%Gamma) * domega(i))

        dfer_arr(rbin) = dfer_arr(rbin) + ret%dFe(1)

        ! Add to the transfer function integral
        do fbin = 1, args%nf
            cexp = cmplx(cos(real(2.d0*pi*tau*args%fi(fbin))),                 &
                sin(real(2.d0*pi*tau*args%fi(fbin))))

            ret%w0(1, gbin, fbin, mubin, rbin) = ret%w0(1, gbin, fbin,         &
                mubin, rbin)+ real(ret%dFe(1))*cexp

            ! TODO: the below are all particular to the lamppost corona, and do
            ! not apply to the ring-like corona currently

            ret%w1(1, gbin, fbin, mubin, rbin) = ret%w1(1, gbin, fbin,         &
                mubin, rbin)+ real(log(gsd))*real(ret%dFe(1))*cexp

            ret%w2(1, gbin, fbin, mubin, rbin) = ret%w2(1, gbin, fbin,         &
                mubin, rbin)+ emissivity*normfac*cexp

            ret%w3(1, gbin, fbin, mubin, rbin) = ret%w3(1, gbin, fbin,         &
                mubin, rbin)+ emissivity*normfac*cexp
        end do

        ! find the appropriate energy and time bins
        tbin = clamp_i(ceiling(log10(tau / time_axis(0)) / args%dlogt), 1,     &
            args%nt)

        ! kernel of the impulse response function
        response(gbin_resp, tbin) = response(gbin_resp, tbin) + ret%dFe(1)
    end do
end subroutine sum_ringlike_corona

subroutine sum_multiple_lampposts(args, i, non_relativistic, r_length,         &
     phi_length, re, alpha, beta, taudo, g, r_grid, domega, gbin, rbin,        &
     ret)
    use dyn_gr
    use radial_grids
    use gr_continuum
    use constants
    use emissivities
    use impulseresponse
    use m_rtrans
    implicit none

    logical, intent(in) :: non_relativistic
    type(t_impulse_args), intent(in) :: args
    integer, intent(in) :: i, r_length, phi_length, gbin, rbin
    double precision, intent(in) :: r_grid(r_length)
    double precision, intent(in) :: domega(r_length)
    double precision, intent(in) :: alpha, beta, taudo, g, re

    type(t_outputs), intent(inout) :: ret

    ! time grid bits, should be passed in
    integer :: tbin, fbin

    ! functions
    double precision :: newtex, dglpfacthick, demang, interper, dareafac,      &
         pfunc_raw
    integer :: get_index, clamp_i

    double precision :: sin0, sindisk
    double precision :: phie
    double precision :: cosfac, mus, ptf
    real :: kfac, normfac, emisfac
    real :: thetafac(args%nlp), gsd(args%nlp)
    ! photon time from/to
    ! tauso is in `gr_continuum`
    double precision :: tausd, mue
    double precision :: tau(args%nlp), emissivity(args%nlp)
    integer :: m, mubin, kk, gbin_resp
    complex :: cexp

    do m = 1, args%nlp

        kk = get_index(rlp(:, m), ndelta, re, args%rmin, npts(m))

        ! Time lag between direct and reflected photons
        if (non_relativistic) then
            tau(m) = sqrt(re                                                   &
                **2                                                            &
                +(args%h(m)                                                    &
                -args%honr                                                     &
                *re)**2)                                                       &
                - re                                                           &
                *(sin0*sindisk*cos(phie)                                       &
                +args%mu0*args%mudisk)+ args%h(1)*args%mu0
            tau(m) = (1.d0+args%zcos)*tau(m)
        else
            ! Interpolate (or extrapolate) the time function
            tausd = interper(rlp(:, m), tlp(:, m), ndelta, re, kk)
            tau(m) = (1.d0+args%zcos)*(tausd+taudo-tauso(1))
        endif

        ! Interpolate |dcos\delta/dr| function
        cosfac = interper(rlp(:, m), dcosdr(:, m), ndelta, re, kk)
        mus = interper(rlp(:, m), cosd(:, m), ndelta, re, kk)

        if (kk .eq. npts(m)) then
            cosfac = newtex(rlp(:, m), dcosdr(:, m), ndelta, re, args%h(m),    &
                args%honr, kk)
            mus = newtex(rlp(:, m), cosd(:, m), ndelta, re, args%h(m),         &
                args%honr, kk)
        end if

        ! Calculate angular emissivity
        ptf = pnorm * pfunc_raw(-mus, args%b1, args%b2, args%qboost)

        ! Calculate flux from pixel
        gsd(m) = dglpfacthick(re, args%spin, args%h(m), args%mudisk)

        ! TODO: write into emissivity(:, m) where the emissivity
        ! now holds the times as well from the time-dependent emissivity
        ! functions
        emissivity(m) = determine_emissivity(re, args%spin, args%gamma,        &
            cosfac, ptf, gsd(m))

        ! calculate extra factors that go into the transfer functions
        ! for double lps
        if (args%nlp .gt. 1) then
            thetafac(m) = emissivity(m)                                        &
                *gso(m)**(args%Gamma-2.)* gsd(m)**(2.-args%Gamma)
        else
            ! single lamp post case, double check this later
            thetafac(m) = 1.
        endif
        ! Add to reflection fraction
        ret%frobs(m) = ret%frobs(m)+ 2.0*g**3*gsd(m)*cosfac/dareafac(re,       &
            args%spin)*domega(i)
    enddo

    ! This loop combines information from the two lampposts, and so
    ! cannot be merged with the above, as it needs to known information
    ! simultaneously

    ! TODO: extend to add a loop over the emissivity time, except that
    ! it is not compatible with the existing use, as the dimensions of
    ! the emissivity would be different
    do m = 1, args%nlp
        ret%dFe(m) = emissivity(m)                                             &
            * (g/(1.d0+args%zcos))**(2.+args%Gamma)* domega(i)

        ! Add to the radial dependence of the transfer function TBD MAKE
        ! SURE THIS IS RIGHT

        dfer_arr(rbin) = dfer_arr(rbin) + ret%dFe(m)
        ! Calculate emission angle and work out which mue bin to add to
        mue = demang(args%spin, args%mu0, re, alpha, beta)
        mubin = ceiling(mue * dble(args%me))
        !calculate the extra factors for w2/3
        if (args%nlp .gt. 1) then
            emisfac = (emissivity(1)+args%eta_0*emissivity(2))/ (1.+args%eta_0)

            kfac = (emissivity(1)                                              &
                +args%eta_0                                                    &
                *emissivity(2))/ (thetafac(1)+args%eta_0*thetafac(2))
        else
            emisfac = emissivity(1)
            kfac = emissivity(1)
            ! single lamp post case, double check this later
        endif

        normfac = real((g/(1.d0+args%zcos))**(2.+args%Gamma)*domega(i))
        ! Add to the transfer function integral
        do fbin = 1, args%nf
            cexp = cmplx(cos(real(2.d0*pi*tau(m)*args%fi(fbin))),              &
                sin(real(2.d0*pi*tau(m)*args%fi(fbin))))

            ret%w0(m, gbin, fbin, mubin, rbin) = ret%w0(m, gbin, fbin,         &
                mubin, rbin)+ real(ret%dFe(m))*cexp

            ret%w1(m, gbin, fbin, mubin, rbin) = ret%w1(m, gbin, fbin,         &
                mubin, rbin)+ real(log(gsd(m)))*real(ret%dFe(m))*cexp

            ! tbd redo these transfer functions
            ret%w2(m, gbin, fbin, mubin, rbin) = ret%w2(m, gbin, fbin,         &
                mubin, rbin)+ emisfac*normfac*cexp

            ret%w3(m, gbin, fbin, mubin, rbin) = ret%w3(m, gbin, fbin,         &
                mubin, rbin)+ kfac*thetafac(m)*normfac*cexp
        end do
    end do

    do m = 1, args%nlp
        ! find the appropriate energy and time bins
        gbin_resp = clamp_i(ceiling(g/args%dg), 1, args%ne)
        tbin = clamp_i(ceiling(log10(tau(m) / time_axis(0)) / args%dlogt),     &
            1, args%nt)

        ! kernel of the impulse response function
        response(gbin_resp, tbin) = response(gbin_resp, tbin) + ret%dFe(m)
    end do
end subroutine sum_multiple_lampposts

!-----------------------------------------------------------------------
function newtex(rlp, dcosdr, ndelta, re, h, honr, kk)
! Extrapolates using Newtonian value
  implicit none
  integer ndelta, kk
  double precision newtex, rlp(ndelta), dcosdr(ndelta), re, h, honr
  newtex = dcosdr(kk) * ((h-honr*rlp(kk))**2 + rlp(kk)**2)**1.5 / rlp(kk)
  newtex = newtex * re / ((h-honr*re)**2 + re**2)**1.5
  return
end function newtex
!-----------------------------------------------------------------------
