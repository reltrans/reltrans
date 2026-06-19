
module telematrix
  !Module containing definitions needs to fold around the telescope
  !response matrix
  logical              :: needchans, needresp, arf, needbkg
  integer              :: nenerg, numchn, Ilo, Ihi, needEs
  real                 :: Elo, Ehi
  real,    allocatable :: En(:), resp(:,:), ECHN(:)
  integer, allocatable :: NGRP(:), FCHAN(:,:), LCHAN(:,:), NCHAN(:,:)
  integer, allocatable :: bkgcounts(:)
  real, allocatable    :: bkgrate(:)
  character (len=500) respname, arfname, bkgname
  data needresp/.true./
  data needchans/.true./
  data needbkg/.true./
end module telematrix

module telematrix2
  !Module containing definitions needs to fold around the telescope
  !response matrix
  logical              :: needchans2, needresp2, arf2
  integer              :: nenerg2, numchn2, Ilo2, Ihi2, needEs2
  real                 :: Elo2, Ehi2
  real,    allocatable :: En2(:), resp2(:,:), ECHN2(:)
  integer, allocatable :: NGRP2(:), FCHAN2(:,:), LCHAN2(:,:), NCHAN2(:,:)
  character (len=500) respname2, arfname2

  data needresp2/.true./
  data needchans2/.true./
contains

  subroutine get_needresp2(check_needresp2) bind(C, name="get_needresp2")
  !---------------------------------------------------------------------------
  !Assign a value to check_needresp2 based on needresp2 variable in the module
  !
  !If needresp2 is True, the second telescope response has NOT be initialised 
  !If needresp2 is False, the second relescope response has been initialised 
  !
  !check_needresp2 = 1 when needresp2 = True
  !check_needresp2 = 0 when needresp2 = False
  !
  !---------------------------------------------------------------------------
    use iso_c_binding
    integer(c_int), intent(out) :: check_needresp2
    if (needresp2) then
       check_needresp2 = 1
    else
       check_needresp2 = 0
    endif
  end subroutine get_needresp2
      
end module telematrix2

module env_variables
  implicit none
  integer :: adensity, idum
  save idum
end module env_variables

module saved_variables
  implicit none
  double precision spinsav,musav,routsav,mudsav
  save spinsav,musav,routsav,mudsav
end module saved_variables

MODULE dyn_gr
!---------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of the GR ray tracing arrays
!  re1:    radius where each photons from the camera hit the disk  
!  taudo1: time of arrival from the observe to the disk (distance is subtracted)  
!  pem1:   value of root of equation \mu(p)= mu for p, that defines
!          where the photons is on the geodesic 
!          pem=-1.D0, if the photon goto infinity.
!          pem=-2.D0, if the photon fall into event horizon.       
!---------------------------------------------------------------------
    implicit none
    integer, parameter :: ndelta = 1000
    integer         , dimension(:)  , allocatable :: npts
    double precision, dimension(:,:), allocatable :: re1, taudo1, pem1
    double precision, dimension(:,:), allocatable :: dcosdr, cosd, rlp, tlp

  contains

!the follwing functions are being exposed for unit testing purposes
    subroutine allocate_re1(x,y) bind(C, name="allocate_re")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: x,y 
      if (allocated(re1)) deallocate(re1)
      allocate(re1(x,y))
    end subroutine allocate_re1

    subroutine allocate_taudo1(x,y) bind(C, name="allocate_taudo")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: x,y 
      if (allocated(taudo1)) deallocate(taudo1)
      allocate(taudo1(x,y))
    end subroutine allocate_taudo1

    subroutine allocate_pem1(x,y) bind(C, name="allocate_pem")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: x,y 
      if (allocated(pem1)) deallocate(pem1)
      allocate(pem1(x,y))
    end subroutine allocate_pem1

    subroutine get_re1(out, nro, nphi) bind(C, name="get_re")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: nro, nphi
      real(c_double), intent(out) :: out(nro*nphi)
      out = reshape(re1, [nro*nphi])
    end subroutine get_re1
    
    subroutine get_taudo1(out, nro, nphi) bind(C, name="get_taudo")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: nro, nphi
      real(c_double), intent(out) :: out(nro*nphi)
      out = reshape(taudo1, [nro*nphi])
    end subroutine get_taudo1

    subroutine get_pem1(out, nro, nphi) bind(C, name="get_pem")
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: nro, nphi
      real(c_double), intent(out) :: out(nro*nphi)
      out = reshape(pem1, [nro*nphi])
    end subroutine get_pem1
  
END MODULE dyn_gr

module xillver_tables
    implicit none

    interface
        ! An interface to the xsatbl function in libXSFunctions.
        subroutine xsatbl(ear,ne,params,filename,ifl,photar,photer)          &
            bind(C, name="xsatbl")
            use iso_c_binding, only: c_float, c_int, c_char
            real(c_float), dimension(*), intent(in) :: ear, params
            real(c_float), dimension(*), intent(out) :: photar, photer
            character(kind=c_char), intent(in) :: filename(*)
            integer(c_int), value, intent(in) :: ne, ifl
        end subroutine xsatbl
    end interface

    character (len=50), parameter ::  xillver = 'xillver-a-Ec5_normalised.fits'
    character (len=50), parameter ::  xillverDCp = 'xillverCp_v3.4_normalised.fits'
    character (len=200) ::  path_tables 
    character (len=200) ::  pathname_xillver 
    character (len=200) ::  pathname_xillverD 
    character (len=200) ::  pathname_xillverDCp
    character (len=500) ::  path_name_reflionx_table
end module xillver_tables

module gr_continuum
  implicit none
  double precision, dimension(:), allocatable :: tauso, gso, lens, cosdelta_obs
  save lens
end module gr_continuum

module radial_grids
  implicit none
  double precision , dimension(:), allocatable :: logxir, logner, gsdr, dfer_arr
  double precision :: pnorm
  save logxir, gsdr, logner
end module radial_grids

module conv_mod
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  !include <libfftw3.a>

  integer, parameter :: nex = 2**12, nex_conv = 4 * nex, nec = nex_conv/2 + 1
  ! real   , dimension(2 * nex_conv) :: adata,bdata,cdata
  ! complex, dimension(nex_conv) :: ac,bc,cc
  
  double precision, parameter :: nexm1 = 1. / real(nex_conv, kind(8))

  type(C_ptr) :: plan1, plan2
  real(   c_double), pointer, dimension(:) :: in,  out_conv
  complex(c_double), pointer, dimension(:) :: out, in_conv
  type(C_ptr) :: a1, a2, a3, a4


contains
  
  subroutine init_fftw_allconv(use_estimate)
    implicit none
    integer(c_int) :: flags, i
    integer, external :: omp_get_max_threads
    logical, intent(in) :: use_estimate
    INTEGER FFTW_PATIENT
    PARAMETER (FFTW_PATIENT=32)
    !i = fftw_init_threads()
    !call fftw_plan_with_nthreads(omp_get_max_threads())
    !print*, "Using threads num:", omp_get_max_threads()

    ! allocate
    a1 = fftw_alloc_real(   int(nex_conv, c_size_t))
    a2 = fftw_alloc_real(   int(nex_conv, c_size_t))
    a3 = fftw_alloc_complex(int(nec     , c_size_t))
    a4 = fftw_alloc_complex(int(nec     , c_size_t))

    call c_f_pointer(a1, in      , [nex_conv])
    call c_f_pointer(a2, out_conv, [nex_conv])
    call c_f_pointer(a3, out     , [nec     ])
    call c_f_pointer(a4, in_conv , [nec     ])
    if (use_estimate) then
        ! Saves time on running the test suite.
        flags = FFTW_ESTIMATE
    else
        flags = FFTW_PATIENT
    endif

    ! note: these two are what kill the runtime of this subroutine
    plan1 = fftw_plan_dft_r2c_1d(nex_conv,  in, out, flags)
    plan2 = fftw_plan_dft_c2r_1d(nex_conv, in_conv, out_conv, flags)
  end subroutine init_fftw_allconv

  subroutine conv_one_FFTw(dyn,earx,Gamma,photarx,reline,imline,ReW_conv,ImW_conv,DC,nlp)
    implicit none
    integer, intent(in) :: DC, nlp 
    real                :: dyn
    real, intent(in)    :: photarx(nex), earx(0:nex), Gamma
    real, intent(in)    :: reline(nlp,nex), imline(nlp,nex)
    real, intent(inout) :: ReW_conv(nlp,nex), ImW_conv(nlp,nex)
    complex :: conv(nec),padFT_photarx(nec)
    complex :: padFT_reline(nec),  padFT_imline(nec)            
    integer :: m, i
    real    :: photmax, depad_conv(nex), E
    
    do m=1,nlp  
       if (DC .eq. 1 ) then
          ! call padding4FT(photarx,padFT_photarx)
          call padding4FT_xillver(photarx,padFT_photarx)

          call padding4FT(reline(m,:),padFT_reline)                        

          conv = (padFT_photarx * padFT_reline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
                   
       else 
          call padding4FT(photarx       , padFT_photarx)        
          call padding4FT(reline(m,:),padFT_reline)
          call padding4FT(imline(m,:),padFT_imline)

          conv = (padFT_photarx * padFT_reline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do

          conv = (padFT_photarx * padFT_imline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ImW_conv(m,i) = ImW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
          
       endif
    end do

  end subroutine conv_one_FFTw

  subroutine padding4FT(line, padFT_line)
    implicit none 
    real           , intent(in)  :: line(nex)
    complex        , intent(out) :: padFT_line(nec)
    
    integer :: i 

    ! Fill padded arrays
    in(1) = 0.0
    do i = 1, nex
        in(i+1) = line(i)
    end do
    do i = 2, 3 * nex
        in(i + nex) = 0.
    end do

    ! do i=1, nex_conv
    !    write(82,*) i, in(i)
    ! enddo
    ! write(82,*) 'no no'

    call fftw_execute_dft_r2c(plan1, in, out)

    padFT_line = out


  end subroutine padding4FT

  subroutine padding4FT_xillver(line, padFT_line)
    implicit none 
    real           , intent(in)  :: line(nex)
    complex        , intent(out) :: padFT_line(nec)

    integer :: i
    real    :: m

    ! Fill padded arrays
    ! in(1) = 0.0
    m = 0.1
    do i = 1, 641
       in(i) = m * real(i) + line(642) - (m * 642)
       if (in(i) .lt. 0.0) then
          in(i) = 0.0
       endif
    end do
    do i = 642, nex
        in(i) = line(i)
    end do
    do i = 1, 3 * nex
        in(i + nex) = 0.
    end do

    ! do i=1, nex_conv
    !    write(82,*) i, in(i)
    ! enddo
    ! write(82,*) 'no no'

    call fftw_execute_dft_r2c(plan1, in, out)

    padFT_line = out

  end subroutine padding4FT_xillver

  subroutine de_paddingFT(dyn, padFT_line, out_line)
    implicit none 
    real    , intent(in) :: dyn
    complex , intent(in) :: padFT_line(nec)
    real    , intent(out):: out_line(nex)

    integer :: i 
    real    :: photmax

    in_conv = padFT_line

    call fftw_execute_dft_c2r(plan2, in_conv, out_conv)
    ! do i = 1, nex_conv 
    !    write(80,*) i, out_conv(i)
    ! enddo
    
    ! Populate output array
    photmax = 0.0
    do i = 1, nex
       out_line(i) = out_conv(i + nex/2 + 1)
       ! write(81,*) i, out_line(i)
        photmax = max( photmax , out_line(i) )
    end do

    ! Clean any residual edge effects
    do i = 1, nex
        if( abs(out_line(i)) .lt. abs(dyn * photmax) ) out_line(i) = 0.0
    end do

    return 
   end subroutine de_paddingFT

end module conv_mod
