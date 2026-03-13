module common_types
    implicit none

    ! this should be moved so that all of the wrappers pass the following
    ! instead of a params array
    type :: t_model_arguments
        ! these are hardcoded with a fixed size of 2, since nlp <= 2. by not
        ! havin them be pointers or dynamically allocated, it makes reasoning
        ! about the code a little easier.
        double precision :: h(2)
        real :: DelAB(2), g(2)
        double precision :: a, inc, rin, rout, zcos, Gamma, honr
        real :: logxi, Afe, lognep, Cutoff_obs, Cutoff_s, Dkpc, Anorm, beta_p
        real :: Nh, boost, Mass, floHz, fhiHz, DelA
        integer :: nlp, ReIm, resp_matr, Cp
        double precision :: qboost, b1, b2, eta, eta_0
    end type t_model_arguments

    type :: t_config
        integer :: verbose = 0
        ! firstcall: is this the first time the model has been called?
        logical :: firstcall = .true., needtrans = .true., needconv = .true.,test = .false.
        integer :: me, xe, nex, m, ionvar, refvar

        ! TODO: are these really constants, or are they hidden variables?
        ! constants
        integer :: nphi = 200, nro = 200
        ! Emin, Emax: min and max of the internal energy grid
        ! (different from output grid)
        real :: Emin = 1e-2, Emax = 3e3, dyn = 0.0

        ! rnmax: max radius for which GR ray-tracing is used
        ! dlogf: a resolution parameter in base 10
        double precision :: rnmax = 300.d0, dlogf = 0.09

        real :: DeltaGamma = 0.01

        ! internal frequency grid
        ! number of frequency bins
        integer :: nf
        real :: f, fac
        double precision :: fc, flo, fhi
        ! internal frequency grid, for when we do lag/frequency spectra
        integer :: fbinx
        real, allocatable :: fix(:)

        ! relativistic parameters and limit on rin and h
        double precision :: rmin, rh
        double precision, allocatable :: height(:)

        ! internal energy grid (nex) and output/xspec (ne) energy grid
        ! dloge: logarithmic resolution of the internal energy grid
        real :: E, dE, dloge

        ! Radial and angle profile
        integer :: mubin, rbin, ibin

        ! variable for non linear effects
        integer :: DC, ionvariation
        real :: dlogxi1, dlogxi2
    end type t_config

    type :: t_arrays
        ! earx: internal energy grid array (0:nex)
        real, allocatable :: earx(:), ear(:), fix(:)
        real, allocatable :: ReGbar(:), ImGbar(:)
        real, allocatable :: contx(:,:)
        double precision, allocatable :: contx_int(:)
        ! lens needs to be allocatable to save it.
        double precision, allocatable :: frobs(:), frrel(:)
        ! TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
        complex, dimension(:,:,:,:,:), allocatable :: ker_W0, ker_W1, ker_W2, ker_W3
        real, dimension(:,:,:), allocatable :: ReW0, ImW0, ReW1, ImW1
        real, dimension(:,:,:), allocatable :: ReW2, ImW2, ReW3, ImW3
        real, dimension(:,:), allocatable :: ReSraw, ImSraw, ReSrawa, ImSrawa, ReGrawa, ImGrawa, ReG, ImG
    end type t_arrays
contains

    ! Unwraps the arguments from a parameter array into `args`.
    subroutine unwrap_arguments(args, nlp, dset, params, cutoff_powerlaw)
        integer, intent(in) :: nlp, dset, cutoff_powerlaw
        real, target, intent(in) :: params(32)
        type(t_model_arguments), intent(out) :: args
        integer :: i
        do i = 1,nlp
            args%DelAB(i) = params(27 + (i - 1) * nlp)
            args%g(i) = params(28 + (i - 1) * nlp)
        end do
        if (dset .eq. 1) then
            args%Dkpc = params(9)
        end if
        args%h(1) = dble(params(1))
        args%h(2) = dble(params(2))
        args%nlp = nlp
        args%a = dble(params(3))
        args%inc = dble(params(4))
        args%rin = dble(params(5))
        args%rout = dble(params(6))
        args%zcos = dble(params(7))
        args%Gamma = dble(params(8))
        args%logxi = params(9)
        args%Afe = params(10)
        args%lognep = params(11)
        args%Cutoff_s = params(12)
        args%Cutoff_obs = params(12)
        args%eta_0 = params(13)
        args%eta = params(14)
        args%beta_p = params(15)
        args%Nh = params(16)
        args%boost = params(17)
        args%qboost = dble(params(18))
        args%Mass = dble(params(19))
        args%honr = dble(params(20))
        args%b1 = dble(params(21))
        args%b2 = dble(params(22))
        args%floHz = params(23)
        args%fhiHz = params(24)
        args%ReIm = int(params(25))
        args%DelA = params(26)
        args%Anorm = params(31)
        args%resp_matr = params(32)
        args%Cp = cutoff_powerlaw
    end subroutine unwrap_arguments

    ! Adjust the model parameters to sane values and set the derived values in
    ! `config`. It performs the following steps:
    ! - Checks `a`, `rin`, `h` are in bounds.
    ! - Sets the inner radius to the ISCO.
    subroutine arguments_check(config, model_args)
        type(t_config), intent(inout) :: config
        type(t_model_arguments), intent(inout) :: model_args
        integer :: i
        double precision :: disco

        ! TODO: should this be printing if it modifies the input parameters?
        ! some kind of warning perhaps?
        if (abs(model_args%a) .gt. 0.999) then
            model_args%a = sign(model_args%a,1.d0) * 0.999
        end if
        config%rmin = disco(model_args%a)
        config%rh = 1.d0+sqrt(1.d0-model_args%a**2)
        if (model_args%rin .lt. 0.d0) then
            model_args%rin = abs(model_args%rin) * config%rmin
        end if
        if (model_args%rin .lt. config%rmin)then
            write(*,*)"Warning! rin<ISCO! Set to ISCO"
            model_args%rin = config%rmin
        end if
        do i=1,model_args%nlp
            if (model_args%h(i) .lt. 0.d0) then
                model_args%h(i) = abs(model_args%h(i)) * config%rh
            end if
            if (model_args%h(i) .lt. 1.5d0*config%rh)then
                write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
                model_args%h(i) = 1.5d0 * config%rh
            end if
        end do
    end subroutine arguments_check

    ! Read in environment variables that configure reltrans
    subroutine read_environment_variables(config)
        use env_variables, only: adensity, idum
        use xillver_tables, only: path_tables, xillver, xillverDCp,pathname_xillver, pathname_xillverDCp
        type(t_config), intent(inout) :: config
        integer :: get_env_int
        character (len=200) :: get_env_char
        config%me = get_env_int("MU_ZONES", 1)
        config%xe = get_env_int("ION_ZONES", 20)
        ! verbose:
        ! 0: XSPEC output only
        ! 1: Print quantities to terminal + 0
        ! 2: Model components, radial scalings, impulse responses written to
        ! file + 1
        config%verbose = get_env_int("REV_VERB", 0)
        ! include pivoting reflection
        config%refvar = get_env_int("REF_VAR", 1)
        ! include ionisation changes
        config%ionvar = get_env_int("ION_VAR", 1)

        ! these are set in `env_variables`
        ! seed for simulation
        idum = get_env_int("SEED_SIM", -2851043)
        ! decide between zone A density profile or constant density profile
        adensity = max(min(get_env_int("A_DENSITY", 0), 1), 0)

        ! this is from xillver_tables, sets the paths where the tables are read
        ! from
        path_tables = get_env_char("RELTRANS_TABLES", './')
        write(pathname_xillver, '(A, A, A)') trim(path_tables), '/', trim(xillver)
        write(pathname_xillverDCp, '(A, A, A)') trim(path_tables), '/', trim(xillverDCp)

        write(*,*)"----------------------------------------------------"
        write(*,*)" *** ENVIRONMENT VARIABLES *** "
        write(*,*)
        write(*,*) 'RADIAL ZONES', config%xe
        write(*,*) 'ANGLE ZONES', config%me
        if (adensity .eq. 0.0) then
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is constant'
        else
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is zone A SS73'
        endif
        write(*,*) 'VERBOSE is ', config%verbose
        write(*,*) 'REFVAR is ', config%refvar
        write(*,*) 'IONVAR is ', config%ionvar 
        write(*,*)"----------------------------------------------------"

      end subroutine read_environment_variables

    ! Initialise all of the configuration fields that can be derived after
    ! `read_environment_variables` has been called, and allocate the arrays
    ! in `arrays`
    subroutine setup_arrays(config, arrays, nlp)
        use conv_mod, only: nex
        type(t_config), intent(inout) :: config
        type(t_arrays), intent(inout) :: arrays
        integer, intent(in) :: nlp
        integer :: i

        allocate(arrays%earx(0:nex))
        allocate(arrays%ReGbar(nex))
        allocate(arrays%ImGbar(nex))

        allocate(arrays%contx(nex,nlp))
        allocate(arrays%contx_int(nlp))

        config%dloge = log10(config%Emax / config%Emin) / float(nex)

        ! populate the energy array
        do i = 0,nex
           arrays%earx(i) = config%Emin                                        &
               * (config%Emax/config%Emin)**(float(i)/float(nex))
        end do

    end subroutine setup_arrays

    ! Reallocate arrays depending on whether they need to be resized
    subroutine realloc_arrays(config, model_args, arrays, prev_nf)
        use conv_mod, only: nex
        type(t_config), intent(in) :: config
        type(t_model_arguments), intent(in) :: model_args
        type(t_arrays), intent(inout) :: arrays
        integer, intent(in) :: prev_nf
        integer :: i
        logical :: needs_allocating

        if (allocated(arrays%fix)) then
            needs_allocating = prev_nf .ne. config%nf
        else
            needs_allocating = .true.
        endif

        if (needs_allocating) then
            if (allocated(arrays%fix)) deallocate(arrays%fix)
            allocate(arrays%fix(0:config%nf))

            ! populate the frequency array
            do i = 0, config%nf
                arrays%fix(i) = model_args%floHz                               &
                    *(model_args%fhiHz                                         &
                    / model_args%floHz)**(real(i) / real(config%nf))
            end do

            ! reallocate the transfer function arrays
            if (allocated(arrays%ker_W0)) deallocate(arrays%ker_W0)
            allocate(arrays%ker_W0(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W1)) deallocate(arrays%ker_W1)
            allocate(arrays%ker_W1(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W2)) deallocate(arrays%ker_W2)
            allocate(arrays%ker_W2(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ker_W3)) deallocate(arrays%ker_W3)
            allocate(arrays%ker_W3(model_args%nlp,nex,config%nf,config%me,config%xe))

            if (allocated(arrays%ReW0)) deallocate(arrays%ReW0)
            allocate(arrays%ReW0(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW0)) deallocate(arrays%ImW0)
            allocate(arrays%ImW0(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW1)) deallocate(arrays%ReW1)
            allocate(arrays%ReW1(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW1)) deallocate(arrays%ImW1)
            allocate(arrays%ImW1(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW2)) deallocate(arrays%ReW2)
            allocate(arrays%ReW2(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW2)) deallocate(arrays%ImW2)
            allocate(arrays%ImW2(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReW3)) deallocate(arrays%ReW3)
            allocate(arrays%ReW3(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ImW3)) deallocate(arrays%ImW3)
            allocate(arrays%ImW3(model_args%nlp,nex,config%nf))

            if (allocated(arrays%ReSraw)) deallocate(arrays%ReSraw)
            allocate(arrays%ReSraw(nex,config%nf))

            if (allocated(arrays%ImSraw)) deallocate(arrays%ImSraw)
            allocate(arrays%ImSraw(nex,config%nf))

            if (allocated(arrays%ReSrawa)) deallocate(arrays%ReSrawa)
            allocate(arrays%ReSrawa(nex,config%nf))

            if (allocated(arrays%ImSrawa)) deallocate(arrays%ImSrawa)
            allocate(arrays%ImSrawa(nex,config%nf))

            if (allocated(arrays%ReGrawa)) deallocate(arrays%ReGrawa)
            allocate(arrays%ReGrawa(nex,config%nf))

            if (allocated(arrays%ImGrawa)) deallocate(arrays%ImGrawa)
            allocate(arrays%ImGrawa(nex,config%nf))

            if (allocated(arrays%ReG)) deallocate(arrays%ReG)
            allocate(arrays%ReG(nex,config%nf))

            if (allocated(arrays%ImG)) deallocate(arrays%ImG)
            allocate(arrays%ImG(nex,config%nf))
        end if
    end subroutine
end module common_types
