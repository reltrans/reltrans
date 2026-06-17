! gfortran post_processor.f90 -J ./build/cache -lcfitsio -L"${HEADAS}/lib" -lreltrans -L./build/lib

program post_processor
  implicit none
  real param(15),par(32)
  integer xe,adensity,verbose,Cp,dset,ne,nlp,ifl
  parameter (ne=10)
  real ear(0:ne),photar(ne)
  
  
! Input parameters
  param(1)  = 11.4144          !h
  param(2)  = 0.998            !a
  param(3)  = 42.2816          !inc
  param(4)  = -2.95166         !rin
  param(5)  = 2e4              !rout
  param(6)  = 0.0              !zcos
  param(7)  = 1.70545          !Gamma
  param(8)  = 3.39811          !logxi
  param(9)  = 2.64304          !Afe
  param(10) = 19.0404          !lognep
  param(11) = 359.389          !kTe
  param(12) = 0.0              !Nh
  param(13) = 0.294556         !boost
  param(14) = 14.8223          !Mass
  param(15) = 0.0999336        !Anorm

! Settings
  xe       = 20       !Number of radial zones
  adensity = 1        !1 = zone A ne; 0 = const ne
  verbose  = 0
  
! Pack into general array
  call dummy_reltransDCp(param,par,Cp,dset)

! Call the main subroutine
  ear = 0.0
  nlp = 1
! call genreltrans0(Cp, dset, nlp, xe, adensity, ear, ne, par, ifl, photar)

  call genreltrans1(Cp, dset, nlp, ear, ne, par, ifl, photar)

end program post_processor







! -----------------------------------------------------------------------
subroutine genreltrans1(Cp, dset, nlp, ear, ne, param, ifl, photar)
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
    implicit none
    ! Constants
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, dlogf = 0.09 !This is a resolution parameter (base 10)
    ! Args:
    integer, intent(inout) :: ifl, dset
    integer, intent(in) :: Cp, ne, nlp
    real, intent(inout) :: param(32)
    real, intent(out) :: photar(ne)
    ! Variables of the subroutine
    ! initializer
    integer :: m, prev_nf, Cpsave, i, j, Cp_cont
    double precision :: d, distance
    real :: f, fac, dE, ear(0:ne)
    ! relativistic parameters and limit on rin and h
    ! lens needs to be allocatable to save it.
    double precision, allocatable :: frobs(:), frrel(:)
    real :: photerx(nex), absorbx(nex), ReS(ne), ImS(ne)
    double precision :: fhisave, flosave, fcons, contx_temp

    double precision, allocatable :: logxir_dens(:)

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
        model_args%ReIm = 1
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
       call radfuncs_dist(config%xe, model_args%rin, rnmax,model_args%b1,     &
            model_args%b2, model_args%qboost, fcons,                           &
            & dble(model_args%lognep), model_args%a, model_args%h(1),          &
            model_args%honr, rlp, dcosdr, cosd, ndelta, config%rmin,npts(1),   &
            & logxir, gsdr, logner, pnorm)
    end if

    call radfunctions_dens(config, model_args, arrays)
    allocate( logxir_dens(config%xe) )
    logxir_dens = logxir
    write(123,*)"skip on"
    do i = 1,config%xe
       write(123,*)i,logxir_dens(i)
    end do
    write(123,*)"no no"

    dset = 1
    model_args%Dkpc = 1.0
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    write(*,*)"fcons=",fcons

    call radfuncs_dist(config%xe, model_args%rin, rnmax,model_args%b1,     &
         model_args%b2, model_args%qboost, fcons,                           &
         & dble(model_args%lognep), model_args%a, model_args%h(1),          &
         model_args%honr, rlp, dcosdr, cosd, ndelta, config%rmin,npts(1),   &
         & logxir, gsdr, logner, pnorm)

    do i = 1,config%xe
       write(123,*)i,logxir(i)
       write(124,*)i,logxir_dens(i)-logxir(i)
    end do

    Distance = 10.0**( 0.5 * ( logxir_dens(2)-logxir(2) )  ) / boost
    write(*,*)"Distance (kpc) = ",distance

    !  ! do this for each lamp post, then find some sort of weird average?
    ! if (config%verbose .gt. 0) write(*,*)"Observer's reflection fraction for each source:",model_args%boost*frobs
    ! if (config%verbose .gt. 0) write(*,*)"Relxill reflection fraction for each source:", frrel

    ! if (config%verbose .gt. 2) call CPU_TIME (time_start)
    ! if (config%needconv)then
    !     call do_convolutions(config, model_args, arrays)
    ! end if
    ! if (config%verbose .gt. 2) then
    !     call CPU_TIME (time_end)
    !     print *, 'Convolutions runtime: ', time_end - time_start, ' seconds'
    ! endif

    ! ! Calculate absorption
    ! call tbabs(arrays%earx, nex, model_args%nh, Ifl, absorbx, photerx)

    ! ! TBD coherence check - if zero coherence between lamp posts, call a
    ! ! different subroutine
    ! if (model_args%ReIm .eq. 7) then
    !     ! tbd - implement zero cohernece in lag_freq
    !     if (nlp .gt. 1 .and. model_args%beta_p .eq. 0.) then
    !         call lag_freq_nocoh(nex, arrays%earx, config%nf, arrays%fix,       &
    !              real(config%flo), real(config%fhi), config%Emin,              &
    !              config%Emax, nlp, arrays%contx, absorbx, real(tauso),         &
    !              real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1,             &
    !              arrays%ImW1, arrays%ReW2, arrays%ImW2, arrays%ReW3,           &
    !              arrays%ImW3, real(model_args%h), real(model_args%zcos),       &
    !              real(model_args%Gamma), real(model_args%eta),                 &
    !              model_args%boost, model_args%g, model_args%DelAB,             &
    !              config%ionvar, arrays%ReGbar, arrays%ImGbar)
    !     else
    !         call lag_freq(nex, arrays%earx, config%nf, arrays%fix,             &
    !              real(config%flo), real(config%fhi), config%Emin,              &
    !              config%Emax, nlp, arrays%contx, absorbx, real(tauso),         &
    !              real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1,             &
    !              arrays%ImW1, arrays%ReW2, arrays%ImW2, arrays%ReW3,           &
    !              arrays%ImW3, real(model_args%h), real(model_args%zcos),       &
    !              real(model_args%Gamma), real(model_args%eta),                 &
    !              model_args%beta_p, model_args%boost, model_args%g,            &
    !              model_args%DelAB, config%ionvar, arrays%ReGbar,               &
    !              arrays%ImGbar)
    !     end if
    ! else if (nlp .gt. 1 .and. model_args%beta_p .eq. 0.) then
    !     call rawG(nex, arrays%earx, config%nf, real(config%flo),               &
    !          real(config%fhi), nlp, arrays%contx, absorbx, real(tauso),        &
    !          real(gso), arrays%ReW0, arrays%ImW0, arrays%ReW1, arrays%ImW1,    &
    !          arrays%ReW2, arrays%ImW2, arrays%ReW3, arrays%ImW3,               &
    !          real(model_args%h), real(model_args%zcos),                        &
    !          real(model_args%Gamma), real(model_args%eta), model_args%boost,   &
    !          model_args%ReIm, model_args%g, model_args%DelAB, config%ionvar,   &
    !          config%DC, model_args%resp_matr, arrays%ReGrawa,                  &
    !          arrays%ImGrawa)
    ! else
    !     ! Calculate raw FT of the full spectrum without absorption
    !     call rawS(nex, arrays%earx, config%nf, real(config%flo),               &
    !          real(config%fhi), nlp, arrays%contx, real(tauso), real(gso),      &
    !          arrays%ReW0, arrays%ImW0, arrays%ReW1, arrays%ImW1,               &
    !          arrays%ReW2, arrays%ImW2, arrays%ReW3, arrays%ImW3,               &
    !          real(model_args%h), real(model_args%zcos),                        &
    !          real(model_args%Gamma), real(model_args%eta),                     &
    !          model_args%beta_p, model_args%boost, model_args%g,                &
    !          model_args%DelAB, config%ionvar, config%DC, arrays%ReSraw,        &
    !          arrays%ImSraw)

    !     ! Include absorption in the model
    !     do j = 1, config%nf
    !         do i = 1, nex
    !             arrays%ReSrawa(i, j) = arrays%ReSraw(i, j) * absorbx(i)
    !             arrays%ImSrawa(i, j) = arrays%ImSraw(i, j) * absorbx(i)
    !         end do
    !     end do
    ! end if

    ! if (config%DC .eq. 1)then
    !     ! Norm is applied internally for DC/time averaged spectrum component of
    !     ! dset=1
    !     ! No need for the immaginary part in DC
    !     do i = 1, nex
    !         arrays%ReGbar(i) = (model_args%Anorm / real(1. +                   &
    !             model_args%eta)) * arrays%ReSrawa(i, 1)
    !     end do
    ! else if (model_args%ReIm .eq. 7) then
    !     ! if calculating the lag-frequency spectrum, just rebin the arrays
    !     call rebinE(arrays%fix, arrays%ReGbar, config%nf, ear, ReS, ne)
    !     call rebinE(arrays%fix, arrays%ImGbar, config%nf, ear, ImS, ne)
    ! else
    !     ! In this case, calculate the lag-energy spectrum
    !     ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band
    !     ! parameters
    !     ! note: this must be done by rawG for two incoherent lamp posts, hence
    !     ! the skip below
    !     if (nlp .eq. 1 .or. model_args%beta_p .ne. 0.) then
    !         if (model_args%ReIm .gt. 0.0) then
    !             call propercross(nex, config%nf, arrays%earx,                  &
    !                  arrays%ReSrawa, arrays%ImSrawa, arrays%ReGrawa,           &
    !                  arrays%ImGrawa, model_args%resp_matr)
    !         else
    !             call propercross_NOmatrix(nex, config%nf, arrays%earx,         &
    !                  arrays%ReSrawa, arrays%ImSrawa, arrays%ReGrawa,           &
    !                  arrays%ImGrawa)
    !         endif
    !     end if
    !     ! Apply phase correction parameter to the cross-spectral model (for bad
    !     ! calibration)
    !     ! this is where coherence = 0 or = 1 cases merge back
    !     do j = 1, config%nf
    !         do i = 1, nex
    !             arrays%ReG(i, j) = cos(model_args%DelA) *                      &
    !                 arrays%ReGrawa(i, j) - sin(model_args%DelA) *              &
    !                 arrays%ImGrawa(i, j)
    !             arrays%ImG(i, j) = cos(model_args%DelA) *                      &
    !                 arrays%ImGrawa(i, j) + sin(model_args%DelA) *              &
    !                 arrays%ReGrawa(i, j)
    !         end do
    !     end do
    !     arrays%ReGbar = 0.0
    !     arrays%ImGbar = 0.0
    !     fac = 2.302585 * config%fc**2 * log10(model_args%fhiHz /               &
    !         model_args%floHz) / ((model_args%fhiHz - model_args%floHz) *       &
    !         real(config%nf))
    !     do j = 1, config%nf
    !         f = model_args%floHz * (model_args%fhiHz /                         &
    !             model_args%floHz)**((real(j) - 0.5) / real(config%nf))
    !         do i = 1, nex
    !             arrays%ReGbar(i) = arrays%ReGbar(i) + arrays%ReG(i, j) / f
    !             arrays%ImGbar(i) = arrays%ImGbar(i) + arrays%ImG(i, j) / f
    !         end do
    !     end do
    !     ! This means that norm for the AC components in the dset=1 model is
    !     ! power in squared fractional rms format
    !     ! note: the factor eta is to have the same normalization as the single
    !     ! LP model, it's 100% arbitrary
    !     arrays%ReGbar = arrays%ReGbar * fac * (model_args%Anorm / real(1.      &
    !         + model_args%eta))**2
    !     arrays%ImGbar = arrays%ImGbar * fac * (model_args%Anorm / real(1.      &
    !         + model_args%eta))**2
    ! end if

    
    ! !------------------------------------------
    ! ! Write output depending on ReIm parameter
    ! !------------------------------------------

    ! !Preparing ReS and ImS: rebin and folding (if needed)
    ! select case (abs(model_args%ReIm)) 

    ! case(1:4)
    !    call crebin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, ear,   &
    !          ReS, ImS) !S is in photar form        
    ! case(5:6)
    !    call cfoldandbin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, &
    !             ear, ReS, ImS, model_args%resp_matr) !S is count rate
    ! case(8)
    !    if( needresp2 ) call initmatrix2
    !    call cfoldandbin(nex, arrays%earx, arrays%ReGbar, arrays%ImGbar, ne, &
    !          ear, ReS, ImS, 2) !folds the spectrum over the second response
    ! end select


    ! !Generate the correct output 
    ! select case (abs(model_args%ReIm))

    ! case(1)          !Real part
    !    photar = ReS

    ! case(2)          !Imaginary part
    !    photar = ImS

    ! case(3,5)        !Modulus
    !    photar = sqrt(ReS**2 + ImS**2)
    !    if (model_args%ReIm==3) then 
    !       write(*, *) "Warning ReIm = 3 should not be used for fitting!"
    !    end if

    ! case(4,6)        !Time lag
    !    do i = 1, ne
    !       dE = ear(i) - ear(i-1)
    !       photar(i) = atan2(ImS(i), ReS(i)) / (2.0*pi*config%fc) * dE
    !    end do
    !    if (model_args%ReIm==4) then 
    !       write(*, *)"Warning ReIm = 4 should not be used for fitting!"
    !    end if

    ! case(7)       !Lag Frequency spectrum
    !    do i = 1, ne
    !       dE = ear(i) - ear(i-1)
    !       photar(i) = atan2(ImS(i), ReS(i))/(pi*(ear(i) + ear(i-1)))*dE
    !    end do
    ! end select

    
    ! !--------------------------------------------
    ! ! If REV_VERB > 1 write components into files
    ! !--------------------------------------------

    ! if (config%verbose .gt. 1 .and. abs(model_args%ReIm) .gt. 0 .and. model_args%ReIm .lt. 7) then
    !     if (config%DC .eq. 0 .and. model_args%beta_p .eq. 0) then
    !        call write_components(ne, ear, nex, arrays%earx, config%nf,         &
    !             real(config%flo), real(config%fhi), nlp, arrays%contx,         &
    !             absorbx, real(tauso), real(gso), arrays%ReW0, arrays%ImW0,     &
    !             arrays%ReW1, arrays%ImW1, arrays%ReW2, arrays%ImW2,            &
    !             arrays%ReW3, arrays%ImW3, real(model_args%h),                  &
    !             real(model_args%zcos), real(model_args%Gamma),                 &
    !             real(model_args%eta), model_args%beta_p, model_args%boost,     &
    !             model_args%floHz, model_args%fhiHz, model_args%ReIm,           &
    !             model_args%DelA, model_args%DelAB, model_args%g,               &
    !             config%ionvar, model_args%resp_matr)
    !     ! catch case here for coherence = 0 or 1
    !     end if
    !     ! this writes the full model as returned to Xspec
    !     ! note that xspec gets output in e.g. lags*dE, and we want just the
    !     ! lags, so a factor dE needs to be included
    !     ! add writing of components for lag frequency spectrum
    !     open (unit = 14, file = 'Output/Total.dat', status = 'replace', action = 'write')
    !     do i = 1, ne
    !         dE = ear(i) - ear(i-1)
    !         write (14, *) (ear(i)+ear(i-1))/2., photar(i)/dE
    !     end do
    !     close(14)
    !     ! print continuum for both single and multiple LPs REDO THIS
    !     open (unit = 24, file = 'Output/Continuum_spec.dat', status = 'replace', action = 'write')
    !     do i = 1, nex
    !         dE = arrays%earx(i) - arrays%earx(i-1)
    !         if (nlp .eq. 1) then
    !             contx_temp = arrays%contx(i, 1)/dE
    !         else
    !             contx_temp = 0.
    !             do m = 1, nlp
    !                 contx_temp = contx_temp + arrays%contx(i, m)
    !             end do
    !             contx_temp = contx_temp/((1.+model_args%eta)*dE)
    !         end if
    !         write (24, *) (arrays%earx(i)+arrays%earx(i-1))/2., contx_temp
    !     end do
    !     close(24)
    ! else if (model_args%ReIm .eq. 7) then
    !    open (unit = 14, file = 'Output/Total.dat', status = 'replace', action = 'write')
    !     do i = 1, ne
    !         dE = ear(i) - ear(i-1)
    !         write (14, *) (ear(i)+ear(i-1))/2., photar(i)/dE
    !     end do
    !     close(14)
    ! endif

    fhisave = config%fhi
    flosave = config%flo
    prev_nf = config%nf
    paramsave = param
    Cpsave = model_args%Cp
end subroutine genreltrans1
! -----------------------------------------------------------------------




  
  
! -----------------------------------------------------------------------
subroutine genreltrans0(Cp, dset, nlp, xein, adensityin, ear, ne, param, ifl, photar)
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
    implicit none
    ! Constants
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, dlogf = 0.09 !This is a resolution parameter (base 10)
    ! Args:
    integer, intent(inout) :: ifl
    integer, intent(in) :: Cp, dset, ne, nlp, xein, adensityin
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

    real Dkpc0, get_fcons

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

    !Set up environment variables (passed directly into the subroutine)
    config%xe = xein
    adensity = adensityin

    
        write(*,*) 'RADIAL ZONES', config%xe        
        if (adensity .eq. 0.0) then
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is constant'
        else
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is zone A SS73'
        endif

        
    !Set up arrays
    call setup_global_arrays(config, model_args%nlp)
    call setup_arrays(config, arrays, model_args%nlp)

    ifl = 1

    !Check I have the right parameters
    write(*,*)"logxi=",model_args%logxi
    write(*,*)"h=",model_args%h(1)
    write(*,*)"a=",model_args%a
    write(*,*)"inc=",model_args%inc    
    write(*,*)"rin=",model_args%rin
    write(*,*)"rout=",model_args%rout
    write(*,*)"zcos=",model_args%zcos
    write(*,*)"Gamma=",model_args%Gamma
    write(*,*)"boost=",model_args%boost
    write(*,*)"Anorm=",model_args%Anorm
    
    ! Note: the two different calls are because for the double lP we set the
    ! temperature from the coronal frame(s), but for the single
    ! LP we use the temperature in the observer frame

    ! ! Decide if this is the DC component/time averaged spectrum or not
    ! if (config%flo .lt. tiny(config%flo) .or. config%fhi .lt. tiny(config%fhi))then
    !     config%DC = 1
    !     model_args%g = 0.0
    !     model_args%DelAB = 0.0
    !     model_args%DelA = 0.0
    !     model_args%ReIm = 1
    !     model_args%eta = model_args%eta_0
    !     ! this is an ugly hack for the double LP model to calculate the time-
    !     ! averaged spectrum
    !     model_args%beta_p = 1.
    ! else
    !     config%DC = 0
    !     model_args%boost = abs(model_args%boost)
    ! end if



    ! set up the continuum spectrum plus relative quantities (cutoff
    ! energies, lensing/gfactors, luminosity, etc)

    Dkpc0 = 1.0
    

    
    write(*,*)"Calling init_cont"
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    write(*,*)"Out of init_cont"

    
    ! if (dset .eq. 0) then
    !    call radfunctions_dens(config, model_args, arrays)
    ! else
    !    call radfuncs_dist(config%xe, model_args%rin, rnmax,model_args%b1,     &
    !         model_args%b2, model_args%qboost, fcons,                           &
    !         & dble(model_args%lognep), model_args%a, model_args%h(1),          &
    !         model_args%honr, rlp, dcosdr, cosd, ndelta, config%rmin,npts(1),   &
    !         & logxir, gsdr, logner, pnorm)
    !  end if

end subroutine genreltrans0
! -----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine dummy_reltransDCp(param,par,Cp,dset)
  implicit none
  integer :: Cp, dset
  real    :: param(15),par(32)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = 0.0              !floHz
  par(24) = 0.0              !fhiHz
  par(25) = 1.0              !ReIm
  par(26) = 0.0              !DelA
  par(27) = 0.0              !DelAB
  par(28) = 0.0              !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(15)        !Anorm
  par(32) = 1.0              !telescope response  
  return
end subroutine dummy_reltransDCp
!-----------------------------------------------------------------------


