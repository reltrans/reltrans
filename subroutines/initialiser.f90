!-----------------------------------------------------------------------
subroutine initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans&
     ,nphi,nro,me,ge,xe)
! Initialises the model and writes the header
  use conv_mod
  use dyn_gr
      implicit none
      integer i, nphi, nro 
      integer  me,ge,xe,myenv
      logical firstcall,needtrans
      real Emin,Emax,dloge,earx(0:nex)
      double precision :: d,rnmax
      needtrans = .false.
      
      if( firstcall )then

        call init_fftw_allconv() !call the initializer of the fftw convolution
         
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS: a transfer function model for"
        write(*,*)"X-ray reverberation mapping."
        write(*,*)"Please cite Ingram et al (2019) MNRAS 488 p324-347."
        write(*,*)"----------------------------------------------------"

!Create *logarithmic* working energy grid
!Will need to evaluate xillver on this grid to use the FT convolution code 
        Emax  = 1e3
        Emin  = 1e-1
        dloge = log10( Emax / Emin ) / float(nex)
        do i = 0,nex
          earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
        end do

        nro   = 200    !resolution variables - these could be made parameters
        nphi  = 200    !  "
        me      = myenv("MU_ZONES"  , 5 )     !Set number of mu_e zones used
        ge      = myenv("ECUT_ZONES", 5 )   !Set number of Ecut zones used
        xe      = myenv("ION_ZONES" , 50)   !Set number of ionisation zones used

        rnmax= 300.d0
! Set sensible distance for observer from the BH
        d = max( 1.0d4 , 2.0d2 * rnmax**2 )
           
        firstcall = .false.
     end if
     return
    end subroutine initialiser
!-----------------------------------------------------------------------
