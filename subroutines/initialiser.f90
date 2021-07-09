!-----------------------------------------------------------------------
subroutine initialiser(firstcall, Emin, Emax, dloge, earx, rnmax, d, needtrans, me, xe, verbose)
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

  !   Last change: Gullo - 2020 Jun
!!!-------------------------------------------------------------------  
  use conv_mod
  use dyn_gr
      implicit none
      integer          , intent(out)   :: xe, me, verbose
      ! integer          , intent(in)    :: nphi, nro !constant
      real             , intent(in)    :: Emin, Emax ! constant
      real             , intent(out)   :: dloge, earx(0:nex)
      double precision , intent(in)    :: rnmax
      double precision , intent(out)   :: d
      logical          , intent(inout) :: firstcall, needtrans
      integer i
      integer myenv

      needtrans = .false.     
      if( firstcall )then

!call the initializer of the FFtw convolution
! the function is in amodules.f90 it sets the structure for the FFTw       
        call init_fftw_allconv() 
         
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS v0.8.2: a transfer function model for"
        write(*,*)"X-ray reverberation mapping."
        write(*,*)"Please cite Ingram et al (2019) MNRAS 488 p324-347."
        write(*,*)"----------------------------------------------------"

!Create *logarithmic* working energy grid
!Will need to evaluate xillver on this grid to use the FT convolution code 
        dloge = log10( Emax / Emin ) / float(nex)
        do i = 0, nex
          earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
        end do

        me      = myenv("MU_ZONES"  , 5 )   !Set number of mu_e zones used
        xe      = myenv("ION_ZONES" , 50)   !Set number of ionisation zones used
! Call environment variables
        verbose = myenv("REV_VERB",0)     !Set verbose level

        write(*,*) 'RADIAL ZONES', xe
        write(*,*) 'ANGLE ZONES', me
        write(*,*) 'VERBOSE is ', verbose

! Set sensible distance for observer from the BH
        d = max( 1.0d4 , 2.0d2 * rnmax**2 )
        
           
        firstcall = .false.
     end if
     return
    end subroutine initialiser
!-----------------------------------------------------------------------
