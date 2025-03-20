!-----------------------------------------------------------------------
subroutine initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans,me,xe,refvar,ionvar,nlp,verbose, test)
!!!  Initialises the model and writes the header
!!!------------------------------------------------------------------
  !    Args:
  !        firstcall: check if this is the first time the model is called
  !        Emin, Emax: (constant) minimum and maximum range of the internal energy grid which is different than the xspec one 
  !        dloge: logarithmic resolution of the internal energy grid
  !        earx:  internal energy grid array (0:nex) [nex is shared variable in conv_mod]
  !        d, rnmax: distance of the source, max radius for which GR ray tracing is used
  !        needtrans: check if transfer function has to be calculated
  !        me, xe: number of angle and radial zones
  !        verbose: check if the verbose env variable is active
  !        nphi, nro: (constant) resolution variables, number of pixels on the observer's camera(b and phib)

  !    Internal variables:
  !        i: loop index

  !   Last change: Gullo - 2024 Oct
!!!-------------------------------------------------------------------  
  use conv_mod
  use dyn_gr
  use env_variables
  use xillver_tables
  use radial_grids
  use gr_continuum
      implicit none
      integer          , intent(out)   :: xe,me,refvar,ionvar,verbose
      integer          , intent(in)    :: nlp !constant
      real             , intent(in)    :: Emin, Emax ! constant
      real             , intent(out)   :: dloge, earx(0:nex)
      double precision , intent(in)    :: rnmax
      double precision , intent(out)   :: d
      logical          , intent(inout) :: firstcall, needtrans, test
      integer i, env_test
      integer get_env_int
      character (len=200) :: get_env_char
 
      needtrans = .false.     
      if( firstcall )then

!call the initializer of the FFtw convolution
! the function is in amodules.f90 it sets the structure for the FFTw       
        call init_fftw_allconv() 
         
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS v1.0.0: a transfer function model for"
        write(*,*)"X-ray reverberation mapping."
        write(*,*)"Please cite Ingram et al (2019) MNRAS 488 p324-347, "
        write(*,*)"Mastroserio et al (2021) MNRAS 507 p55-73, and "  
        write(*,*)"Lucchini et al (2023) arXiv 230505039L."  
        write(*,*)"----------------------------------------------------"

!Create *logarithmic* working energy grid
!Will need to evaluate xillver on this grid to use the FT convolution code 
        dloge = log10( Emax / Emin ) / float(nex)
        do i = 0, nex
           earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
        end do

        ! Call environment variables
        me      = get_env_int("MU_ZONES"  , 1 )   !Set number of mu_e zones used
        xe      = get_env_int("ION_ZONES" , 20)   !Set number of ionisation zones used
        ! Decide between zone A density profile or constant density profile
        adensity = get_env_int("A_DENSITY",0)
        adensity = min( adensity , 1 )
        adensity = max( adensity , 0 )
        verbose = get_env_int("REV_VERB",0)     !Set verbose level
                                          !0: Xspec output only
                                          !1: Also print quantities to terminal
                                          !2: Also print model components, radial scalings and impulse response function to 
                                          !files in /Output folder
        refvar = get_env_int("REF_VAR",1)         !choose whether to include pivoting reflection
        ionvar = get_env_int("ION_VAR",1)         !choose whether to include ionization changes
        idum = get_env_int("SEED_SIM", -2851043)  !seed for simulations

        write(*,*) 'RADIAL ZONES', xe
        write(*,*) 'ANGLE ZONES', me
        if (adensity .eq. 0.0) then
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is constant'
        else
            write(*,*) 'A_DENSITY:', adensity, 'Density profile is zone A SS73'
        endif
        write(*,*) 'VERBOSE is ', verbose
        write(*,*) 'REFVAR is ', refvar
        write(*,*) 'IONVAR is ', ionvar 

! set if it's a TEST run 
        env_test = get_env_int("TEST_RUN",0)   
        if (env_test .eq. 1) then
           write(*,*) '*********  This is a TEST run  *********'
           test = .true.
        else
           test = .false.
        endif
        if (test) then 
           call FNINIT
        endif

        write(*,*)"----------------------------------------------------"

! Set sensible distance for observer from the BH
        d = max( 1.0d4 , 2.0d2 * rnmax**2 )

! set the table names 
        path_tables = get_env_char("RELTRANS_TABLES"  , './' )   !search for the env variable RELTRANS_TABLES otherwise set the path to ./
        write(pathname_xillver, '(A, A, A)') trim(path_tables), '/', trim(xillver)
        write(pathname_xillverD, '(A, A, A)') trim(path_tables), '/', trim(xillverD)
        write(pathname_xillverDCp, '(A, A, A)') trim(path_tables), '/', trim(xillverDCp)
        write(*,'(A, A)') 'Set the XILLVER table to ', trim(pathname_xillver)
        write(*,'(A, A)') 'Set the high density XILLVER table to ', trim(pathname_xillverD)
        write(*,'(A, A)') 'Set the nthComp, high density XILLVER table to ', trim(pathname_xillverDCp)
        
        firstcall = .false.

        !Allocate some useful arrays

        !Allocate arrays for radial profiles 
        if (.not. allocated (dfer_arr)) allocate (dfer_arr(xe))
        if (.not. allocated (logxir)  ) allocate (logxir  (xe))
        if (.not. allocated (gsdr)    ) allocate (gsdr    (xe))
        if (.not. allocated (logner)  ) allocate (logner  (xe))

        !Allocate GR arrays
        if (.not. allocated (cosd)        ) allocate (cosd  (ndelta,nlp))
        if (.not. allocated (dcosdr)      ) allocate (dcosdr(ndelta,nlp))
        if (.not. allocated (rlp   )      ) allocate (rlp   (ndelta,nlp))
        if (.not. allocated (tlp   )      ) allocate (tlp   (ndelta,nlp))
        if (.not. allocated (npts)        ) allocate (npts  (       nlp))
        if (.not. allocated (gso)         ) allocate (gso   (       nlp))
        if (.not. allocated (tauso)       ) allocate (tauso (       nlp))
        if (.not. allocated (cosdelta_obs)) allocate (cosdelta_obs(nlp))

     end if
     return
    end subroutine initialiser
!-----------------------------------------------------------------------
