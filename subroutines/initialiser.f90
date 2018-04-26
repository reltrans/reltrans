!-----------------------------------------------------------------------
      subroutine initialiser(firstcall,Emin,Emax,nex,dloge,earx,needtrans)
! Initialises the model and writes the header
      implicit none
      integer nex,i
      logical firstcall,needtrans
      real Emin,Emax,dloge,earx(0:nex)
      needtrans = .false.
      if( firstcall )then
        needtrans = .true.
        write(*,*)"----------------------------------------------------"
        write(*,*)"This is RELTRANS: a transfer function model for"
        write(*,*)"X-ray reverberation mapping written by Adam Ingram."
        write(*,*)"Please cite Ingram et al (2018)."
        write(*,*)"----------------------------------------------------"
        !Create *logarithmic* working energy grid
        !Will need to evaluate xillver on this grid to use the FT convolution code 
        Emax  = 1e3
        Emin  = 1e-3
        dloge = log10( Emax / Emin ) / float(nex)
        do i = 0,nex
          earx(i) = Emin * (Emax/Emin)**(float(i)/float(nex))
        end do
        firstcall = .false.
      end if
      return
    end subroutine initialiser
!-----------------------------------------------------------------------
