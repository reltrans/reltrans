!-----------------------------------------------------------------------
subroutine initialiser(firstcall,ReIm,DC,nlp,nphi,nro,nf,nfsave,nexsave,nlpsave,Anorm,Emin,Emax,dloge,earxi,rnmax,d,needtrans,&
                       me,xe,refvar,ionvar,verbose)
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
!!!-------------------------------------------------------------------  
    use conv_mod
    use dyn_gr
    use dyn_en
    use env_variables
    implicit none
    integer          , intent(out)   :: xe,me,refvar,ionvar,verbose
    ! integer          , intent(in)    :: nphi, nro !constant
    real             , intent(in)    :: Emin, Emax ! constant
    real             , intent(out)   :: dloge, earxi(0:nexi), Anorm!, earx(0:nex)
    double precision , intent(in)    :: rnmax
    double precision , intent(out)   :: d
    logical          , intent(inout) :: firstcall, needtrans
    integer i, ReIm, DC, nlp, nphi, nro, nf, nfsave, nexsave, conv_size, grid_red, nlpsave
    integer myenv

 
    !this controls how much more coarse the energy grid becomes going into lag energy or lag freqeuency modes
    !and the size of the arrays used in computing the Fourier transform. TWEAK AT YOUR OWN RISK
    grid_red = 4   
    conv_size = 2 
    if( DC .eq. 1 ) then
        nex = nexi
    else if (ReIm .lt. 7) then 
        nex = nexi/grid_red
    else 
        nex = nexi/(grid_red**2)
    end if

    nex_conv = conv_size * nex
    nec = nex_conv/2 + 1
    
    !THIS IS JUST A PLACEHOLDER TO AVOID THE NORMALIZATION ISSUES
    Anorm = Anorm*grid_red**2
    needtrans = .false.     
    if( firstcall )then    
        !call the initializer of the FFtw convolution
        ! the function is in amodules.f90 it sets the structure for the FFTw       
        call init_fftw_allconv(nexi,conv_size,grid_red) 
                
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS v1.0.0: a transfer function model for"
        write(*,*)"X-ray reverberation mapping."
        write(*,*)"Please cite Ingram et al (2019) MNRAS 488 p324-347, "
        write(*,*)"Mastroserio et al (2021) MNRAS 507 p55-73, and "  
        write(*,*)"Lucchini et al (2023) arXiv 230505039L."  
        write(*,*)"----------------------------------------------------"

        ! Call environment variables
        me      = myenv("MU_ZONES"  , 1 )   !Set number of mu_e zones used
        xe      = myenv("ION_ZONES" , 20)   !Set number of ionisation zones used
        ! Decide between zone A density profile or constant density profile
        adensity = myenv("A_DENSITY",0)
        adensity = min( adensity , 1 )
        adensity = max( adensity , 0 )
        verbose = myenv("REV_VERB",0)     !Set verbose level
                                          !0: Xspec output only
                                          !1: Also print quantities to terminal
                                          !2: Also print model components, radial scalings and impulse response function to 
                                          !files in /Output folder
        refvar = myenv("REF_VAR",1)         !choose whether to include pivoting reflection
        ionvar = myenv("ION_VAR",1)         !choose whether to include ionization changes
        idum = myenv("SEED_SIM", -2851043)  !seed for simulations

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
        write(*,*)"----------------------------------------------------"
        
        ! Set sensible distance for observer from the BH
        d = max( 1.0d4 , 2.0d2 * rnmax**2 )
        firstcall = .false.
   ! else if (nlp .ne. nlpsave) then
        !call init_fftw_allconv(nexi,conv_size,grid_red) 
    end if

    call allocate_arrays(nlp,nphi,nro,nf,nfsave,nexsave,nlpsave,me,xe)    

    !First internal energy grid: static grid used to evaluate xspec models
    dloge = log10( Emax / Emin ) / float(nexi)
    do i = 0, nexi
        earxi(i) = Emin * (Emax/Emin)**(float(i)/float(nexi))
    end do

   !Second internal energy grid: dynamic grid that changes size depending on mode
    dloge = log10( Emax / Emin ) / float(nex)
    do i = 0, nex
        earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
    end do       

    return
    end subroutine initialiser
!-----------------------------------------------------------------------
