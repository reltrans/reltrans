module m_genreltrans
    use common_types
    implicit none

    interface
        ! The interface around the XSPEC C_tbabs function.
        ! It will call the C_tbabs symbol in the libXSFunctions shared library.
        subroutine c_tbabs(earx, nex, params, Ifl, absorbx, photerx, str)      &
            bind(C, name = "C_tbabs")
            use iso_c_binding, only: c_double, c_int, c_char
            real(c_double), intent(in) :: earx(nex+1)
            real(c_double), intent(in) :: params(1)
            real(c_double), intent(inout) :: absorbx(nex), photerx(nex)
            integer(c_int), value, intent(in) :: nex, Ifl
            character(kind = c_char), intent(in) :: str(*)
        end subroutine c_tbabs
    end interface

contains

    ! Call the tbabs function from the XSPEC model library.
    !
    ! This wrapper temporarily performs a runtime cast on all of the arrays to
    ! double precision, as the C_tbabs function that we will eventually call
    ! expects double precision, whilst much of reltrans still uses `real`, which
    ! maps to single precision.
    !
    ! Once the precision has been modified, this function can be simplified.
    subroutine tbabs(earx, nex, nh, Ifl, absorbx, photerx)
        real, intent(in) :: earx(0:nex), nh
        real, intent(inout) :: absorbx(nex), photerx(nex)
        integer, intent(in) :: nex, Ifl

        double precision :: d_earx(0:nex), d_absorbx(nex), d_photerx(nex),     &
            d_params(1)
        integer i

        do i = 0, nex
            d_earx(i) = earx(i)
        end do

        d_params(1) = nh
        call c_tbabs(d_earx, nex, d_params, Ifl, d_absorbx, d_photerx, "")

        do i = 1, nex
            absorbx(i) = d_absorbx(i)
            photerx(i) = d_photerx(i)
        end do
    end subroutine tbabs

    ! Allocate the global arrays that reltrans needs
    subroutine setup_global_arrays(config, nlp)
        use dyn_gr, only: ndelta, cosd, dcosdr, rlp, tlp, npts, re1, pem1, taudo1
        use gr_continuum, only: tauso, gso, cosdelta_obs
        use radial_grids, only: logxir, gsdr, dfer_arr, logner
        integer, intent(in) :: nlp
        type(t_config), intent(in) :: config
        ! allocate arrays for radial profiles
        if (allocated(dfer_arr)) deallocate(dfer_arr)
        if (allocated(logxir)) deallocate(logxir)
        if (allocated(gsdr)) deallocate(gsdr)
        if (allocated(logner)) deallocate(logner)
        allocate(dfer_arr(config%xe))
        allocate(logxir (config%xe))
        allocate(gsdr (config%xe))
        allocate(logner (config%xe))

        ! allocate GR arrays
        if (allocated (cosd)) deallocate (cosd)
        if (allocated (dcosdr)) deallocate (dcosdr)
        if (allocated (rlp)) deallocate (rlp)
        if (allocated (tlp)) deallocate (tlp)
        if (allocated (npts)) deallocate (npts)
        if (allocated (gso)) deallocate (gso)
        if (allocated (tauso)) deallocate (tauso)
        if (allocated (cosdelta_obs)) deallocate (cosdelta_obs)
        allocate (cosd(ndelta, nlp))
        allocate (dcosdr(ndelta, nlp))
        allocate (rlp(ndelta, nlp))
        allocate (tlp(ndelta, nlp))
        allocate (npts(nlp))
        allocate (gso(nlp))
        allocate (tauso(nlp))
        allocate (cosdelta_obs(nlp))

        if (allocated(re1)) deallocate(re1)
        if (allocated(taudo1)) deallocate(taudo1)
        if (allocated(pem1)) deallocate(pem1)
        allocate(re1(config%nphi, config%nro))
        allocate(taudo1(config%nphi, config%nro))
        allocate(pem1(config%nphi, config%nro))
    end subroutine setup_global_arrays

    subroutine print_header()
        use reltrans_version_mod
        write(*, *)"----------------------------------------------------"
        write(*, *)"This is RELTRANS " // reltrans_version // ": a transfer function model for"
        write(*, *)"X-ray reverberation mapping."
        write(*, *)"Please cite Ingram et al (2019) MNRAS 488 p324-347, "
        write(*, *)"Mastroserio et al (2021) MNRAS 507 p55-73, and "
        write(*, *)"Lucchini et al (2023) arXiv 230505039L."
        write(*, *)"----------------------------------------------------"
    end subroutine print_header

    ! Calculate the frequency limits and number of frequency bins, depending on
    ! ReIm.
    subroutine config_frequency(config, model_args)
        type(t_config), intent(inout) :: config
        type(t_model_arguments), intent(inout) :: model_args

        logical :: fhizero, flozero

        ! TODO: magic constants

        ! TODO: this is a breaking change compared to the comparison to see if
        ! fhiHz etc are greater than tiny(...)

        ! rework this logic so that low frequencies always result in time
        ! independent spectrum, not just for reim<7
        if (model_args%ReIm .eq. 7) then
            if (model_args%Mass .gt. 1000) then
                model_args%floHz = 1.e-5
                model_args%fhiHz = 5.e-2
            else
                model_args%floHz = 0.07
                model_args%fhiHz = 700.0
            end if
            ! TODO: why 0.01 here, but 0.09 in the other branch?
            ! overwrite the frequency spacing
            config%dlogf = 0.01
        endif
        ! else:
        ! if doing lag-energy spectra, just work out how many frequencies to
        ! average over

        ! Convert frequency bounds from Hz to c/Rg (now being more accurate with
        ! constants)
        config%fhi = dble(model_args%fhiHz) * 4.92695275718945d-06 *           &
            model_args%Mass
        config%flo = dble(model_args%floHz) * 4.92695275718945d-06 *           &
            model_args%Mass
        config%nf = ceiling(log10(model_args%fhiHz / model_args%floHz) /       &
            config%dlogf)
        config%fc = 0.5d0 * (model_args%floHz + model_args%fhiHz)

        ! Check how we did and adjust if necessary
        fhizero = model_args%fhiHz .lt. tiny(model_args%fhiHz)
        flozero = model_args%floHz .lt. tiny(model_args%floHz)
        if (fhizero .or. flozero) then
            model_args%fhiHz = 0.d0
            model_args%floHz = 0.d0
            config%nf = 1
        end if
    end subroutine config_frequency

    subroutine do_convolutions(config, model_args, arrays)
        use conv_mod, only: nex, conv_one_FFTw
        use gr_continuum, only: gso, lens
        use radial_grids, only: logner, gsdr, logxir
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(in) :: model_args
        type(t_arrays), intent(inout) :: arrays

        double precision, parameter :: pi = acos(-1.d0)
        real, parameter :: dyn = 0.0 !1e-7
        real :: Gamma0, Gamma1, Gamma2, Cutoff_0, E, mue, thetae
        real :: Hx_delta(nex), Hx_dlogxi(nex)
        real :: dlogxi1, dlogxi2, logne, logxi0, logxi1, logxi2
        integer :: i, j, m, ionvariation, mubin, rbin
        real :: photarx(nex), photarx_1(nex), photarx_2(nex)
        real :: reline_w0(model_args%nlp, nex), imline_w0(model_args%nlp, nex)
        real :: reline_w1(model_args%nlp, nex), imline_w1(model_args%nlp, nex)
        real :: reline_w2(model_args%nlp, nex), imline_w2(model_args%nlp, nex)
        real :: reline_w3(model_args%nlp, nex), imline_w3(model_args%nlp, nex)
        real :: photarx_dlogxi(nex)
        real :: Hx(nex), photarx_delta(nex)

        ! Initialize arrays for transfer functions
        arrays%ReW0 = 0.0
        arrays%ImW0 = 0.0
        arrays%ReW1 = 0.0
        arrays%ImW1 = 0.0
        arrays%ReW2 = 0.0
        arrays%ImW2 = 0.0
        arrays%ReW3 = 0.0
        arrays%ImW3 = 0.0
        Gamma1 = real(model_args%Gamma) - 0.5*config%DeltaGamma
        Gamma2 = real(model_args%Gamma) + 0.5*config%DeltaGamma
        ! Get logxi values corresponding to Gamma1 and Gamma2
        call xilimits(nex, arrays%earx, model_args%nlp, arrays%contx,          &
             config%DeltaGamma, real(gso), real(lens),                         &
             real(model_args%zcos), dlogxi1, dlogxi2)
        ! Set the ion-variation to 1, there is an if inside the radial loop to
        ! check if either the ionvar is 0 or the logxi is 0 to
        ! set ionvariation to 0  it is important that ionvariation is different
        ! than ionvar because ionvar  is used also later in
        ! the rawS subroutine to calculate the cross-spectrum
        ionvariation = 1
        ! Loop over radius, emission angle and frequency
        do rbin = 1, config%xe !Loop over radial zones
            ! Set parameters with radial dependence
            Gamma0 = real(model_args%Gamma)
            logne = logner(rbin)
            Cutoff_0 = real(gsdr(rbin)) * model_args%Cutoff_s
            logxi0 = real(logxir(rbin))
            if (config%xe .eq. 1)then
                Cutoff_0 = model_args%Cutoff_s
                logne = model_args%lognep
                logxi0 = model_args%logxi
            end if
            ! Avoid negative values of the ionisation parameter
            if (logxi0 .eq. 0.0 .or. config%ionvar .eq. 0) then
                ionvariation = 0.0
            end if
            do mubin = 1, config%me !loop over emission angle zones
                ! Calculate input emission angle
                mue = (real(mubin) - 0.5) / real(config%me)
                thetae = acos(mue) * 180.0 / real(pi)
                if (config%me .eq. 1) thetae = real(model_args%inc)
                ! Call restframe reflection model
                call rest_frame(model_args, arrays, Gamma0, logne,             &
                    Cutoff_0, logxi0, thetae, photarx)
                ! NON LINEAR EFFECTS
                if (config%DC .eq. 0) then
                   ! Gamma variations
                   logxi1 = logxi0 + ionvariation * dlogxi1
                   call rest_frame(model_args, arrays, Gamma1, logne,          &
                       Cutoff_0, logxi1, thetae, photarx_1)
                   logxi2 = logxi0 + ionvariation * dlogxi2
                   call rest_frame(model_args, arrays, Gamma2, logne,          &
                       Cutoff_0, logxi2, thetae, photarx_2)
                   photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
                   ! xi variations
                   call rest_frame(model_args, arrays, Gamma0, logne,          &
                       Cutoff_0, logxi1, thetae, photarx_1)
                   call rest_frame(model_args, arrays, Gamma0, logne,          &
                       Cutoff_0, logxi2, thetae, photarx_2)
                   photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10
                end if
                ! Multiply by E^{Gamma-1} to make less steep
                do i = 1, nex
                   E = 0.5 * (arrays%earx(i) + arrays%earx(i-1))
                   Hx(i) = photarx(i) * E**(Gamma0-1)
                   Hx_delta(i) = photarx_delta(i) * E**(Gamma0-1)
                   Hx_dlogxi(i) = photarx_dlogxi(i) * E**(Gamma0-1)
                end do
                ! Loop through frequencies and lamp posts
                do j = 1, config%nf
                   do i = 1, nex
                      do m = 1, model_args%nlp
                         reline_w0(m, i) = real(arrays%ker_W0(m, i, j,         &
                             mubin, rbin))
                         imline_w0(m, i) = aimag(arrays%ker_W0(m, i, j,        &
                             mubin, rbin))
                         reline_w1(m, i) = real(arrays%ker_W1(m, i, j,         &
                             mubin, rbin))
                         imline_w1(m, i) = aimag(arrays%ker_W1(m, i, j,        &
                             mubin, rbin))
                         reline_w2(m, i) = real(arrays%ker_W2(m, i, j,         &
                             mubin, rbin))
                         imline_w2(m, i) = aimag(arrays%ker_W2(m, i, j,        &
                             mubin, rbin))
                         reline_w3(m, i) = real(arrays%ker_W3(m, i, j,         &
                             mubin, rbin))
                         imline_w3(m, i) = aimag(arrays%ker_W3(m, i, j,        &
                             mubin, rbin))
                      end do
                   end do
                   ! TODO: this test wrapping should be inside the conv functions,
                   ! not at their callsites
                   ! Do the convolution (involves multiplying by E^{1-Gamma})
                   if (config%test) then
                      call conv_one_FFT(dyn, arrays%earx, Gamma0, Hx,          &
                           reline_w0, imline_w0, arrays%ReW0(:, :, j),         &
                           arrays%ImW0(:, :, j), config%DC, model_args%nlp)
                      if (config%DC .eq. 0 .and. config%refvar .eq. 1) then
                         call conv_one_FFT(dyn, arrays%earx, Gamma0, Hx,       &
                              reline_w1, imline_w1, arrays%ReW1(:, :, j),      &
                              arrays%ImW1(:, :, j), config%DC,                 &
                              model_args%nlp)
                         call conv_one_FFT(dyn, arrays%earx, Gamma0,           &
                              Hx_delta, reline_w2, imline_w2, arrays%ReW2(:,   &
                              :, j), arrays%ImW2(:, :, j), config%DC,          &
                              model_args%nlp)
                      end if
                      if (config%DC .eq. 0 .and. config%ionvar .eq. 1) then
                         call conv_one_FFT(dyn, arrays%earx, Gamma0,           &
                              Hx_dlogxi, reline_w3, imline_w3,                 &
                              arrays%ReW3(:, :, j), arrays%ImW3(:, :, j),      &
                              config%DC, model_args%nlp)
                      end if
                   else
                      call conv_one_FFTw(dyn, arrays%earx, Gamma0, Hx,         &
                           reline_w0, imline_w0, arrays%ReW0(:, :, j),         &
                           arrays%ImW0(:, :, j), config%DC, model_args%nlp)
                      if (config%DC .eq. 0 .and. config%refvar .eq. 1) then
                         call conv_one_FFTw(dyn, arrays%earx, Gamma0, Hx,      &
                              reline_w1, imline_w1, arrays%ReW1(:, :, j),      &
                              arrays%ImW1(:, :, j), config%DC,                 &
                              model_args%nlp)
                         call conv_one_FFTw(dyn, arrays%earx, Gamma0,          &
                              Hx_delta, reline_w2, imline_w2, arrays%ReW2(:,   &
                              :, j), arrays%ImW2(:, :, j), config%DC,          &
                              model_args%nlp)
                      end if
                      if (config%DC .eq. 0 .and. config%ionvar .eq. 1) then
                         call conv_one_FFTw(dyn, arrays%earx, Gamma0,          &
                              Hx_dlogxi, reline_w3, imline_w3,                 &
                              arrays%ReW3(:, :, j), arrays%ImW3(:, :, j),      &
                              config%DC, model_args%nlp)
                      end if
                   endif
                end do
            end do
        end do
    end subroutine do_convolutions
end module m_genreltrans

! -----------------------------------------------------------------------
subroutine genreltrans(Cp, dset, nlp, ear, ne, param, ifl, photar)
! All reltrans flavours are calculated in this subroutine.
! Cp and dset are the settings:
! |Cp|=1 means use cut-off power-law, |Cp|=2 means use nthcomp
! Cp>1 means there is a density parameter, Cp<1 means density is hardwired
! dset=0 means ionisation is a parameter, dset=1 means ionization is calculated
! from distance. What to do about ION_ZONES=1 in the distance model?

! The parameter array has 27 parameters. No one model actually has 27
! parameters. In each model, some of these parameters are hardwired, but
! the parameters must be sorted into the param(1:27) array for this subroutine.

! Arg:

! Internal variables:
! constants:
! pi: greek pi
! rnmax: maximum radius to consider GR effects
! nphi, rno: resolution variables, number of pixels on the observer's camera(b
! and phib)
! Emax, Emin: minimum and maximum range of the internal energy grid which is
! different than the xspec one
! dlogf: resolution parameter of the frequency grid
! dyn:   limit to check the saved values
! ionvar: sets the ionisation variation (1 = w/ ion var; 0 = w/o ion var)

    use dyn_gr
    use conv_mod
    use radial_grids
    use gr_continuum
    use m_genreltrans
    use env_variables
    use saved_variables
    use telematrix2
    use rtconstants
    implicit none
    ! Constants
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, dlogf = 0.09 !This is a resolution parameter (base 10)
    ! Args:
    integer, intent(inout) :: ifl
    integer, intent(in) :: Cp, dset, ne, nlp
    real, intent(inout) :: param(32)
    real, intent(out) :: photar(ne)
    ! Variables of the subroutine
    ! initializer
    integer :: m, prev_nf, Cpsave, i, j, Cp_cont
    double precision :: d
    real :: f, fac, dE, ear(0:ne)
    ! relativistic parameters and limit on rin and h
    ! lens needs to be allocatable to save it.
    double precision, allocatable :: frobs(:), frrel(:)
    real :: photerx(nex), absorbx(nex), ReS(ne), ImS(ne)
    double precision :: fhisave, flosave, fcons, contx_temp

    real time_start, time_end !runtime stuff

    ! SAVE
    real, save :: paramsave(32)

    data Cpsave/2/
    data prev_nf /-1/
    ! Save the first call variables
    save d, fhisave, flosave, prev_nf, frobs, frrel, Cpsave

    type(t_config), pointer :: config
    type(t_arrays), save :: arrays
    type(t_model_arguments) :: model_args
    ! make arrays static so its values are kept between function calls

    config => global_config
    call unwrap_arguments(model_args, nlp, dset, param, Cp)
    call config_frequency(config, model_args)
    call arguments_check(config, model_args)

    ! TODO: check to make sure nlp hasn't changed, else many arrays need to be
    ! freed and re-allocated

    if (config%firstcall) then
        call init_fftw_allconv()
        ! initialise environment and allocate all arrays
        call read_environment_variables(config)
        call setup_global_arrays(config, model_args%nlp)
        call setup_arrays(config, arrays, model_args%nlp)

        config%firstcall = .false.
        config%needtrans = .true.
        config%needconv = .true.
        prev_nf = 0 !this is needed to reallocate arrays with realloc_arrays, if firstcall is set to true externally

        ! set sensible distance for observer from the BH
        d = max(1.0d4, 2.0d2 * config%rnmax**2)

        spinsav = -2.d0 !this is needed to force the run of the GRtrace routine

        ! finally, let the people know what they are witnessing!
        call print_header()
    end if

    ! reallocated frequency dependent arrays
    call realloc_arrays(config, model_args, arrays, prev_nf)

    ifl = 1

    ! Note: the two different calls are because for the double lP we set the
    ! temperature from the coronal frame(s), but for the single
    ! LP we use the temperature in the observer frame

    ! Decide if this is the DC component/time averaged spectrum or not
    if (config%flo .lt. tiny(config%flo) .or. config%fhi .lt. tiny(config%fhi))then
        config%DC = 1
        model_args%g = 0.0
        model_args%DelAB = 0.0
        model_args%DelA = 0.0
        model_args%ReIm = MODE_CROSS_SPEC_REAL_REF_FOLDED
        model_args%eta = model_args%eta_0
        ! this is an ugly hack for the double LP model to calculate the time-
        ! averaged spectrum
        model_args%beta_p = 1.
    else
        config%DC = 0
        model_args%boost = abs(model_args%boost)
    end if

    ! Determine if I need to calculate the kernel
    call need_check(model_args%Cp, Cpsave, param, paramsave, config%fhi,       &
         config%flo, fhisave, flosave, config%nf, prev_nf, config%needtrans,   &
         config%needconv)
    if (config%verbose .gt. 2) call CPU_TIME (time_start)
    if (config%needtrans)then
       ! allocate lensing/reflection fraction arrays if necessary
       if (allocated(lens)) deallocate(lens)
       allocate (lens(nlp))
       if (allocated(frobs)) deallocate(frobs)
       allocate (frobs(nlp))
       if (allocated(frrel)) deallocate(frrel)
       allocate (frrel(nlp))
       ! Calculate the Kernel for the given parameters
       call rtrans(config, model_args, arrays, dset, d, nex, frobs, frrel)
       ! print *, 'gso ', gso(1)
    end if
    if (config%verbose .gt. 2) then
       call CPU_TIME (time_end)
       print *, 'Transfer function runtime: ', time_end - time_start, ' seconds'
    end if

    ! set up the continuum spectrum plus relative quantities (cutoff
    ! energies, lensing/gfactors, luminosity, etc)
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    if (dset .eq. 0) then
       call radfunctions_dens(config, model_args, arrays)
    else
       call radfuncs_dist(config, model_args, fcons)
     end if

     ! do this for each lamp post, then find some sort of weird average?
    if (config%verbose .gt. 0) write(*,*)"Observer's reflection fraction for each source:",model_args%boost*frobs
    if (config%verbose .gt. 0) write(*,*)"Relxill reflection fraction for each source:", frrel

    if (config%verbose .gt. 2) call CPU_TIME (time_start)
    if (config%needconv)then
        call do_convolutions(config, model_args, arrays)
    end if
    if (config%verbose .gt. 2) then
        call CPU_TIME (time_end)
        print *, 'Convolutions runtime: ', time_end - time_start, ' seconds'
    endif

    ! Calculate absorption
    call tbabs(arrays%earx, nex, model_args%nh, Ifl, absorbx, photerx)

    ! TBD coherence check - if zero coherence between lamp posts, call a
    ! different subroutine
    if (model_args%ReIm .eq. 7) then
        ! tbd - implement zero cohernece in lag_freq
        if (nlp .gt. 1 .and. model_args%beta_p .eq. 0.) then
            call lag_freq_nocoh(nex, arrays%earx, config%nf, arrays%fix,       &
                 real(config%flo), real(config%fhi), config%Emin,              &
                 config%Emax, nlp, arrays%contx, absorbx, real(tauso),         &
                 real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1,             &
                 arrays%ImW1, arrays%ReW2, arrays%ImW2, arrays%ReW3,           &
                 arrays%ImW3, real(model_args%h), real(model_args%zcos),       &
                 real(model_args%Gamma), real(model_args%eta),                 &
                 model_args%boost, model_args%g, model_args%DelAB,             &
                 config%ionvar, arrays%ReGbar, arrays%ImGbar)
        else
            call lag_freq(nex, arrays%earx, config%nf, arrays%fix,             &
                 real(config%flo), real(config%fhi), config%Emin,              &
                 config%Emax, nlp, arrays%contx, absorbx, real(tauso),         &
                 real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1,             &
                 arrays%ImW1, arrays%ReW2, arrays%ImW2, arrays%ReW3,           &
                 arrays%ImW3, real(model_args%h), real(model_args%zcos),       &
                 real(model_args%Gamma), real(model_args%eta),                 &
                 model_args%beta_p, model_args%boost, model_args%g,            &
                 model_args%DelAB, config%ionvar, arrays%ReGbar,               &
                 arrays%ImGbar)
        end if
    else if (nlp .gt. 1 .and. model_args%beta_p .eq. 0.) then
        call rawG(nex, arrays%earx, config%nf, real(config%flo),               &
             real(config%fhi), nlp, arrays%contx, absorbx, real(tauso),        &
             real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1, arrays%ImW1,    &
             arrays%ReW2, arrays%ImW2, arrays%ReW3, arrays%ImW3,               &
             real(model_args%h), real(model_args%zcos),                        &
             real(model_args%Gamma), real(model_args%eta), model_args%boost,   &
             model_args%ReIm, model_args%g, model_args%DelAB, config%ionvar,   &
             config%DC, model_args%resp_matr, arrays%ReGrawa,                  &
             arrays%ImGrawa)
    else
        ! Calculate raw FT of the full spectrum without absorption
        call rawS(nex, arrays%earx, config%nf, real(config%flo),               &
             real(config%fhi), nlp, arrays%contx, real(tauso), real(gso),      &
             arrays%ReW0, arrays%ImW0, arrays%ReW1, arrays%ImW1,               &
             arrays%ReW2, arrays%ImW2, arrays%ReW3, arrays%ImW3,               &
             real(model_args%h), real(model_args%zcos),                        &
             real(model_args%Gamma), real(model_args%eta),                     &
             model_args%beta_p, model_args%boost, model_args%g,                &
             model_args%DelAB, config%ionvar, config%DC, arrays%ReSraw,        &
             arrays%ImSraw)

        ! Include absorption in the model
        do j = 1, config%nf
            do i = 1, nex
                arrays%ReSrawa(i, j) = arrays%ReSraw(i, j) * absorbx(i)
                arrays%ImSrawa(i, j) = arrays%ImSraw(i, j) * absorbx(i)
            end do
        end do
    end if

    if (config%DC .eq. 1)then
        ! Norm is applied internally for DC/time averaged spectrum component of
        ! dset=1
        ! No need for the immaginary part in DC
        do i = 1, nex
            arrays%ReGbar(i) = (model_args%Anorm / real(1. +                   &
                model_args%eta)) * arrays%ReSrawa(i, 1)
        end do
    else if (model_args%ReIm == MODE_LAG_FREQ) then
        ! if calculating the lag-frequency spectrum, just rebin the arrays
        call rebinE(arrays%fix, arrays%ReGbar, config%nf, ear, ReS, ne)
        call rebinE(arrays%fix, arrays%ImGbar, config%nf, ear, ImS, ne)
    else
        ! In this case, calculate the lag-energy spectrum
        ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band
        ! parameters
        ! note: this must be done by rawG for two incoherent lamp posts, hence
        ! the skip below
        if (nlp == 1 .or. model_args%beta_p .ne. 0.) then
            if (is_folded(model_args%ReIm)) then
                call propercross(nex, config%nf, arrays%earx,                  &
                     arrays%ReSrawa, arrays%ImSrawa, arrays%ReGrawa,           &
                     arrays%ImGrawa, model_args%resp_matr)
            else
                call propercross_NOmatrix(nex, config%nf, arrays%earx,         &
                     arrays%ReSrawa, arrays%ImSrawa, arrays%ReGrawa,           &
                     arrays%ImGrawa)
            endif
        end if
        ! Apply phase correction parameter to the cross-spectral model (for bad
        ! calibration)
        ! this is where coherence = 0 or = 1 cases merge back
        do j = 1, config%nf
            do i = 1, nex
                arrays%ReG(i, j) = cos(model_args%DelA) *                      &
                    arrays%ReGrawa(i, j) - sin(model_args%DelA) *              &
                    arrays%ImGrawa(i, j)
                arrays%ImG(i, j) = cos(model_args%DelA) *                      &
                    arrays%ImGrawa(i, j) + sin(model_args%DelA) *              &
                    arrays%ReGrawa(i, j)
            end do
        end do
        arrays%ReGbar = 0.0
        arrays%ImGbar = 0.0
        fac = 2.302585 * config%fc**2 * log10(model_args%fhiHz /               &
            model_args%floHz) / ((model_args%fhiHz - model_args%floHz) *       &
            real(config%nf))
        do j = 1, config%nf
            f = model_args%floHz * (model_args%fhiHz /                         &
                model_args%floHz)**((real(j) - 0.5) / real(config%nf))
            do i = 1, nex
                arrays%ReGbar(i) = arrays%ReGbar(i) + arrays%ReG(i, j) / f
                arrays%ImGbar(i) = arrays%ImGbar(i) + arrays%ImG(i, j) / f
            end do
        end do
        ! This means that norm for the AC components in the dset=1 model is
        ! power in squared fractional rms format
        ! note: the factor eta is to have the same normalization as the single
        ! LP model, it's 100% arbitrary
        arrays%ReGbar = arrays%ReGbar * fac * (model_args%Anorm / real(1.      &
            + model_args%eta))**2
        arrays%ImGbar = arrays%ImGbar * fac * (model_args%Anorm / real(1.      &
            + model_args%eta))**2
    end if


    !------------------------------------------
    ! Write output depending on ReIm parameter
    !------------------------------------------

    !Preparing ReS and ImS: rebin and folding (if needed)
    if (is_both_folded(model_args%reim)) then
       call cfoldandbin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, &
                ear, ReS, ImS, model_args%resp_matr) !S is count rate
    else if (is_folded(model_args%reim)) then
       call crebin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, ear,   &
             ReS, ImS) !S is in photar form
    end if

    if (model_args%reim == MODE_CROSS_SPEC_LAG_TWO_RESPONSES) then
       if( needresp2 ) call initmatrix2
       call cfoldandbin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, &
             ear, ReS, ImS, 2) !folds the spectrum over the second response
    end if

    !Generate the correct output
    if (is_mode(model_args%reim, MODE_CROSS_SPEC_REAL)) then
       photar = ReS
    else if (is_mode(model_args%reim, MODE_CROSS_SPEC_IMAG)) then
       photar = ImS
    else if (is_mode(model_args%reim, MODE_CROSS_SPEC_MODULUS)) then
       photar = sqrt(ReS**2 + ImS**2)
       if (is_mode(model_args%ReIm, MODE_CROSS_SPEC_MODULUS)) then
          write(*, *) "Warning ReIm = 3 should not be used for fitting!"
       end if
    else if (is_mode(model_args%reim, MODE_CROSS_SPEC_LAG)) then
       do i = 1, ne
          dE = ear(i) - ear(i-1)
          photar(i) = atan2(ImS(i), ReS(i)) / (2.0*pi*config%fc) * dE
       end do
       if (is_mode(model_args%ReIm, MODE_CROSS_SPEC_LAG)) then
          write(*, *)"Warning ReIm = 4 should not be used for fitting!"
       end if
    else if (model_args%reim == MODE_LAG_FREQ) then
       do i = 1, ne
          dE = ear(i) - ear(i-1)
          photar(i) = atan2(ImS(i), ReS(i))/(pi*(ear(i) + ear(i-1)))*dE
       end do
    end if


    !--------------------------------------------
    ! If REV_VERB > 1 write components into files
    !--------------------------------------------

    if (config%verbose .gt. 1 .and. abs(model_args%ReIm) .gt. 0 .and. model_args%ReIm .lt. 7) then
        if (config%DC .eq. 0 .and. model_args%beta_p .eq. 0) then
           call write_components(ne, ear, nex, arrays%earx, config%nf,         &
                real(config%flo), real(config%fhi), nlp, arrays%contx,         &
                absorbx, real(tauso), real(gso), arrays%ReW0, arrays%ImW0,     &
                arrays%ReW1, arrays%ImW1, arrays%ReW2, arrays%ImW2,            &
                arrays%ReW3, arrays%ImW3, real(model_args%h),                  &
                real(model_args%zcos), real(model_args%Gamma),                 &
                real(model_args%eta), model_args%beta_p, model_args%boost,     &
                model_args%floHz, model_args%fhiHz, model_args%ReIm,           &
                model_args%DelA, model_args%DelAB, model_args%g,               &
                config%ionvar, model_args%resp_matr)
        ! catch case here for coherence = 0 or 1
        end if
        ! this writes the full model as returned to Xspec
        ! note that xspec gets output in e.g. lags*dE, and we want just the
        ! lags, so a factor dE needs to be included
        ! add writing of components for lag frequency spectrum
        open (unit = 14, file = 'Output/Total.dat', status = 'replace', action = 'write')
        do i = 1, ne
            dE = ear(i) - ear(i-1)
            write (14, *) (ear(i)+ear(i-1))/2., photar(i)/dE
        end do
        close(14)
        ! print continuum for both single and multiple LPs REDO THIS
        open (unit = 24, file = 'Output/Continuum_spec.dat', status = 'replace', action = 'write')
        do i = 1, nex
            dE = arrays%earx(i) - arrays%earx(i-1)
            if (nlp .eq. 1) then
                contx_temp = arrays%contx(i, 1)/dE
            else
                contx_temp = 0.
                do m = 1, nlp
                    contx_temp = contx_temp + arrays%contx(i, m)
                end do
                contx_temp = contx_temp/((1.+model_args%eta)*dE)
            end if
            write (24, *) (arrays%earx(i)+arrays%earx(i-1))/2., contx_temp
        end do
        close(24)
    else if (is_mode(model_args%ReIm, MODE_LAG_FREQ)) then
       open (unit = 14, file = 'Output/Total.dat', status = 'replace', action = 'write')
        do i = 1, ne
            dE = ear(i) - ear(i-1)
            write (14, *) (ear(i)+ear(i-1))/2., photar(i)/dE
        end do
        close(14)
    endif

    fhisave = config%fhi
    flosave = config%flo
    prev_nf = config%nf
    paramsave = param
    Cpsave = model_args%Cp
end subroutine genreltrans
! -----------------------------------------------------------------------
