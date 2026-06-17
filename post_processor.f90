! gfortran post_processor.f90 -J ./build/cache -lcfitsio -L"${HEADAS}/lib" -lreltrans -L./build/lib

program post_processor
  implicit none
  real param(15),par(32)
  integer xe,adensity,Cp,nlp
  double precision distance,Dkpc
  
  
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
  
! Pack into general array
  call pack_reltransDCp(param,par,Cp,nlp)

! Calculate distance
  Dkpc = distance(Cp, nlp, xe, adensity, par)
  write(*,*)"Dkpc=",Dkpc
  
  Dkpc = distance(Cp, nlp, xe, adensity, par)
  write(*,*)"Dkpc=",Dkpc
  
end program post_processor







! -----------------------------------------------------------------------
function distance(Cp, nlp, xe_in, adensity_in, param)
!
! Calculates the distance given a set of model parameters for reltransDCp.
! Can't do this for reltransPL because the density is hardwired to 15.
! In principle could also do for the double lamppost model, but I haven't
! yet implemented that.  
! Inputs:
! Cp        = PL/nthcomp flag. Only works for Cp=2.
! nlp       = No of lampposts. Code doesn't work for nlp=2.  
! param(32) = full genreltrans parameter array
! Outputs:
! distance  = distance in kpc
!
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
    integer, intent(in) :: Cp, nlp, xe_in, adensity_in
    real, intent(inout) :: param(32)
    ! Variables of the subroutine
    ! initializer
    integer :: m, prev_nf, Cpsave, i, j, Cp_cont, dset, get_env_int
    double precision :: d, distance
    ! relativistic parameters and limit on rin and h
    ! lens needs to be allocatable to save it.
    double precision, allocatable :: frobs(:), frrel(:)
    double precision :: fhisave, flosave, fcons, contx_temp

    double precision, allocatable :: logxir_dens(:)

    real time_start, time_end !runtime stuff

    ! SAVE
    real, save :: paramsave(32)

    data Cpsave/2/
    data prev_nf /-1/
    ! Save the first call variables
    save d, frobs, frrel

    type(t_config), pointer :: config
    type(t_arrays), save :: arrays
    type(t_model_arguments) :: model_args
    ! make arrays static so its values are kept between function calls

    ! Setup
    config => global_config
    call unwrap_arguments(model_args, nlp, dset, param, Cp)
    call config_frequency(config, model_args)
    call arguments_check(config, model_args)

    if (config%firstcall) then
        call init_fftw_allconv()
        config%firstcall = .false.
        config%needtrans = .true.
        config%needconv  = .true.
        prev_nf = 0 !this is needed to reallocate arrays with realloc_arrays, if firstcall is set to true externally
        ! set sensible distance for observer from the BH
        d = max(1.0d4, 2.0d2 * config%rnmax**2)
        spinsav = -2.d0 !this is needed to force the run of the GRtrace routine
        config%verbose = get_env_int("REV_VERB", 0)
        ! initialise environment and allocate all arrays
        config%xe      = xe_in
        config%me      = 1
        adensity       = adensity_in        
        call setup_global_arrays(config, model_args%nlp)
    end if

    ! Allocate arrays that need to be re-allocated every time incase xe changed
    config%xe      = xe_in
    adensity       = adensity_in        
    if (allocated(dfer_arr)) deallocate(dfer_arr)
    if (allocated(logxir)) deallocate(logxir)
    if (allocated(gsdr)) deallocate(gsdr)
    if (allocated(logner)) deallocate(logner)
    allocate(dfer_arr(config%xe))
    allocate(logxir (config%xe))
    allocate(gsdr (config%xe))
    allocate(logner (config%xe))    
    call setup_arrays(config, arrays, model_args%nlp)
    
    ! reallocated frequency dependent arrays
    call realloc_arrays(config, model_args, arrays, prev_nf)

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

    ! Calculate the kernel
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

! Calculate ionization profile measured by reltransDCp
    dset = 0
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    call radfunctions_dens(config, model_args, arrays)
    allocate( logxir_dens(config%xe) )
    logxir_dens = logxir

! Calculate ionization profile assuming isotropic corona seen from D=1kpc
    dset            = 1
    model_args%Dkpc = 1.0
    call init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    call radfuncs_dist(config%xe, model_args%rin, rnmax,model_args%b1,         &
         model_args%b2, model_args%qboost, fcons,                              &
         & dble(model_args%lognep), model_args%a, model_args%h(1),             &
         model_args%honr, rlp, dcosdr, cosd, ndelta, config%rmin,npts(1),      &
         & logxir, gsdr, logner, pnorm)
    
! Caculate the re-scaling distance
    Distance = 10.0**( 0.5 * ( logxir_dens(2)-logxir(2) )  )
    Distance = Distance / sqrt(model_args%boost)

! Save parameters
    fhisave = config%fhi
    flosave = config%flo
    prev_nf = config%nf
    paramsave = param
    Cpsave = model_args%Cp
    
  end function distance  
! -----------------------------------------------------------------------




  


!-----------------------------------------------------------------------
subroutine pack_reltransDCp(param,par,Cp,nlp)
  implicit none
  real, intent(in)     :: param(15)
  real, intent(out)    :: par(32)
  integer, intent(out) :: Cp,nlp
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter
  nlp  = 1   !Number of lampposts
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
end subroutine pack_reltransDCp
!-----------------------------------------------------------------------


