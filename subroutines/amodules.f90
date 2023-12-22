
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
end module telematrix2

module env_variables
  implicit none
  integer :: adensity, idum
  save idum
end module env_variables
  
MODULE dyn_gr
!---------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!---------------------------------------------------------------------
    implicit none
    logical :: status_re_tau  
    double precision,dimension(:,:),allocatable :: re1,taudo1,pem1
    save status_re_tau
END MODULE dyn_gr

module xillver_tables
    implicit none 
    character (len=50), parameter ::  xillver = 'xillver-a-Ec5.fits'
    character (len=50), parameter ::  xillverD = 'xillverD-5.fits'
    character (len=50), parameter ::  xillverDCp = 'xillverCp_v3.4_normalised.fits'
    character (len=200) ::  path_tables 
    character (len=200) ::  pathname_xillver 
    character (len=200) ::  pathname_xillverD 
    character (len=200) ::  pathname_xillverDCp
end module xillver_tables

module conv_mod
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  integer, parameter :: nex = 2**12, nex_conv = 4 * nex, nec = nex_conv/2 + 1
  ! real   , dimension(2 * nex_conv) :: adata,bdata,cdata
  ! complex, dimension(nex_conv) :: ac,bc,cc
  
  double precision, parameter :: nexm1 = 1. / real(nex_conv, kind(8))

  type(C_ptr) :: plan1, plan2
  real(   c_double), pointer, dimension(:) :: in,  out_conv
  complex(c_double), pointer, dimension(:) :: out, in_conv
  type(C_ptr) :: a1, a2, a3, a4


contains
  
    subroutine init_fftw_allconv()
    implicit none
    integer(c_int) :: flags, i
    integer, external :: omp_get_max_threads

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
 !   flags = 0 + FFTW_ESTIMATE
    flags = 0 + FFTW_PATIENT

! note: these two are what kill the runtime of this subroutine
    plan1 = fftw_plan_dft_r2c_1d(nex_conv,  in, out, flags)
    plan2 = fftw_plan_dft_c2r_1d(nex_conv, in_conv, out_conv, flags)
    end subroutine init_fftw_allconv

    subroutine conv_one_FFTw(dyn,photarx,reline,imline,ReW_conv,ImW_conv,DC,nlp)
    implicit none
    integer, intent(in) :: DC, nlp 
    real                :: dyn
    real, intent(in)    :: photarx(nex)
    real, intent(in)    :: reline(nlp,nex), imline(nlp,nex)
    real, intent(inout) :: ReW_conv(nlp,nex), ImW_conv(nlp,nex)
    complex :: conv(nec),padFT_photarx(nec)
    complex :: padFT_reline(nec),  padFT_imline(nec)            
    integer :: m
    real    :: photmax, depad_conv(nex)
    
    do m=1,nlp  
        if (DC .eq. 1 ) then   
            call padding4FT(photarx,padFT_photarx)         
            call padding4FT(reline(m,:),padFT_reline)                        
            
            conv = (padFT_photarx * padFT_reline) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW_conv(m,:) = ReW_conv(m,:) + depad_conv
        else 
            call padding4FT(photarx       , padFT_photarx)        
            call padding4FT(reline(m,:),padFT_reline)
            call padding4FT(imline(m,:),padFT_imline)
     
            conv = (padFT_photarx * padFT_reline) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW_conv(m,:) = ReW_conv(m,:) + depad_conv
            
            conv = (padFT_photarx * padFT_imline) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ImW_conv(m,:) = ImW_conv(m,:) + depad_conv           
        endif
    end do

    end subroutine conv_one_FFTw
 
    subroutine conv_all_FFTw(dyn,photarx,photarx_delta,photarx_dlogxi,reline_w0,imline_w0,reline_w1,imline_w1,& 
                             reline_w2,imline_w2,reline_w3,imline_w3,ReW0_conv,ImW0_conv,ReW1_conv,ImW1_conv,&
                             ReW2_conv,ImW2_conv,ReW3_conv,ImW3_conv,DC,nlp)

    implicit none
    integer, intent(in) :: DC, nlp 
    real                :: dyn
    real, intent(in)    :: photarx(nex), photarx_delta(nex), photarx_dlogxi(nex) 
    real, intent(in)    :: reline_w0(nlp,nex), imline_w0(nlp,nex)
    real, intent(in)    :: reline_w1(nlp,nex), imline_w1(nlp,nex)
    real, intent(in)    :: reline_w2(nlp,nex), imline_w2(nlp,nex)
    real, intent(in)    :: reline_w3(nlp,nex), imline_w3(nlp,nex)
    real, intent(inout) :: ReW0_conv(nlp,nex), ImW0_conv(nlp,nex)
    real, intent(inout) :: ReW1_conv(nlp,nex), ImW1_conv(nlp,nex)
    real, intent(inout) :: ReW2_conv(nlp,nex), ImW2_conv(nlp,nex)
    real, intent(inout) :: ReW3_conv(nlp,nex), ImW3_conv(nlp,nex)
    complex :: conv(nec),padFT_photarx(nec),padFT_photarx_delta(nec),padFT_photarx_dlogxi(nec) 
    complex :: padFT_reline_w0(nec),  padFT_imline_w0(nec)
    complex :: padFT_reline_w1(nec),  padFT_imline_w1(nec)
    complex :: padFT_reline_w2(nec),  padFT_imline_w2(nec)
    complex :: padFT_reline_w3(nec),  padFT_imline_w3(nec)              
    integer :: m
    real    :: photmax, depad_conv(nex)

    do m=1,nlp  
        if (DC .eq. 1 ) then   
            call padding4FT(photarx,padFT_photarx)         
            call padding4FT(reline_w0(m,:),padFT_reline_w0)                        
            
            conv = (padFT_photarx * padFT_reline_w0) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW0_conv(m,:) = ReW0_conv(m,:) + depad_conv
        else 
            call padding4FT(photarx       , padFT_photarx)
            call padding4FT(photarx_delta , padFT_photarx_delta)
            call padding4FT(photarx_dlogxi, padFT_photarx_dlogxi)  
        
            call padding4FT(reline_w0(m,:),padFT_reline_w0)
            call padding4FT(reline_w1(m,:),padFT_reline_w1)
            call padding4FT(reline_w2(m,:),padFT_reline_w2)
            call padding4FT(reline_w3(m,:),padFT_reline_w3)      
            
            call padding4FT(imline_w0(m,:),padFT_imline_w0)
            call padding4FT(imline_w1(m,:),padFT_imline_w1)
            call padding4FT(imline_w2(m,:),padFT_imline_w2)
            call padding4FT(imline_w3(m,:),padFT_imline_w3)      

            conv = (padFT_photarx * padFT_reline_w0) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW0_conv(m,:) = ReW0_conv(m,:) + depad_conv
            
            conv = (padFT_photarx * padFT_imline_w0) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ImW0_conv(m,:) = ImW0_conv(m,:) + depad_conv
            
            conv = (padFT_photarx * padFT_reline_w1) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW1_conv(m,:) = ReW1_conv(m,:) + depad_conv

            conv = (padFT_photarx * padFT_imline_w1) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ImW1_conv(m,:) = ImW1_conv(m,:) + depad_conv  
            
            conv = (padFT_photarx_delta * padFT_reline_w2) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW2_conv(m,:) = ReW2_conv(m,:) + depad_conv

            conv = (padFT_photarx_delta * padFT_imline_w2) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ImW2_conv(m,:) = ImW2_conv(m,:) + depad_conv  
            
            conv = (padFT_photarx_dlogxi * padFT_reline_w3) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ReW3_conv(m,:) = ReW3_conv(m,:) + depad_conv
            
            conv = (padFT_photarx_dlogxi * padFT_imline_w3) * nexm1
            call de_paddingFT(dyn, conv, depad_conv)
            ImW3_conv(m,:) = ImW3_conv(m,:) + depad_conv  
        endif
    end do
        
    end subroutine conv_all_FFTw

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

    call fftw_execute_dft_r2c(plan1, in, out)

    padFT_line = out
    
    end subroutine padding4FT


    subroutine de_paddingFT(dyn, padFT_line, out_line)
    implicit none 
    real    , intent(in) :: dyn
    complex , intent(in) :: padFT_line(nec)
    real    , intent(out):: out_line(nex)

    integer :: i 
    real    :: photmax

    in_conv = padFT_line

    call fftw_execute_dft_c2r(plan2, in_conv, out_conv)
    
    ! Populate output array
    photmax = 0.0
    do i = 1, nex
        out_line(i) = out_conv(i + nex/2 + 1)
        photmax = max( photmax , out_line(i) )
    end do
    ! Clean any residual edge effects
    do i = 1, nex
        if( abs(out_line(i)) .lt. abs(dyn * photmax) ) out_line(i) = 0.0
    end do

    return 
   end subroutine de_paddingFT

end module conv_mod



!***********************************************************************
!This is code of YNOGK used for ray-traycing with General Relativity  
!***********************************************************************

!************************************************************
      module constants
!********************************************************************
!*    This module defines many constants often uesd in our code.
!*    One can use these constants through a command "use constants" in their
!*    own subroutines or functions. 
!********************************************************************
      implicit none
      Double precision  infinity,pi,dtors,sixteen,twopi,zero,one,two,three,four,six,half,&
            half2,mh,hbar,pho_v,plankc,five,dtor,eight        
      parameter(infinity=1.D40,dtors=asin(1.D0)*2.D0/180.D0, &
            sixteen=16.D0, twopi=4.D0*dasin(1.D0), pi = dasin(1.D0)*2.D0)!3.141592653589793D0
      PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0, six=6.D0, half=0.5D0, half2=0.25D0, &
            mh=1.6726231D-24, hbar = 1.0545887D-27, plankc=6.626178D-27, pho_v=2.99792458D10, five=5.D0,&
               dtor=asin(1.D0)*2.D0/180.D0, eight=8.D0)    
!********************************************************************************************
      end module constants 
!********************************************************************************************



!********************************************************************************************
      module rootsfinding
!********************************************************************
!* This module aim on solve cubic and quartic polynomial equations.
!* One can use these subroutine root3 and root4 to find roots of cubic and 
!* quartic equations respectively. 
!********************************************************************
      use constants
      implicit none

      contains
!*********************************************************************************************
      subroutine root3(b,c,d,r1,r2,r3,del)
!********************************************************************
!* PURPOSE:  This subroutine aim on solving cubic equations: x^3+b*x^2+c*x+d=0.
!* INPUTS:    b, c, d-----they are the coefficients of equation. 
!* OUTPUTS:   r1,r2,r3----roots of the equation, with complex number form.
!*            del---------the number of real roots among r1,r2,r3. 
!* ROUTINES CALLED:  sort    
!* This code comes from internet.
!********************************************************************
      implicit none
      Double precision a,b,c,d,p,q,DD,temp1,temp2,phi,y1,y2,&
             y3,y2r,y2i,u,v
      complex*16 r1,r2,r3
      integer del
        
      a=1.D0        
! Step 1: Calculate p and q --------------------------------------------
      p  = c/a - b*b/a/a/3.D0
      q  = (two*b*b*b/a/a/a - 9.D0*b*c/a/a + 27.D0*d/a) / 27.D0

! Step 2: Calculate DD (discriminant) ----------------------------------
      DD = p*p*p/27.D0 + q*q/4.D0

! Step 3: Branch to different algorithms based on DD -------------------
      if(DD .lt. 0.D0)then
!         Step 3b:
!         3 real unequal roots -- use the trigonometric formulation
          phi = acos(-q/two/sqrt(abs(p*p*p)/27.D0))
          temp1=two*sqrt(abs(p)/3.D0)
          y1 =  temp1*cos(phi/3.D0)
          y2 = -temp1*cos((phi+pi)/3.D0)
          y3 = -temp1*cos((phi-pi)/3.D0)
      else
!         Step 3a:
!         1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
          temp1 = -q/two + sqrt(DD)
          temp2 = -q/two - sqrt(DD)
          u = abs(temp1)**(1.D0/3.D0)
          v = abs(temp2)**(1.D0/3.D0)
          if(temp1 .lt. 0.D0) u=-u
          if(temp2 .lt. 0.D0) v=-v
          y1  = u + v
          y2r = -(u+v)/two
          y2i =  (u-v)*sqrt(3.D0)/two
      endif
! Step 4: Final transformation -----------------------------------------
      temp1 = b/a/3.D0
      y1 = y1-temp1
      y2 = y2-temp1
      y3 = y3-temp1
      y2r=y2r-temp1
! Assign answers -------------------------------------------------------
      if(DD .lt. 0.D0)then
          call sort(y1,y2,y3,y1,y2,y3)
          r1 = dcmplx( y1,  0.D0)
          r2 = dcmplx( y2,  0.D0)
          r3 = dcmplx( y3,  0.D0)
          del=3
      elseif(DD .eq. 0.D0)then
          call sort(y1,y2r,y2r,y1,y2r,y2r)
          r1 = dcmplx( y1,  0.D0)
          r2 = dcmplx(y2r,  0.D0)
          r3 = dcmplx(y2r,  0.D0)
          del=3
      else
          r1 = dcmplx( y1,  0.D0)
          r2 = dcmplx(y2r, y2i)
          r3 = dcmplx(y2r,-y2i)
          del=1
      endif
      return
      end subroutine root3       
!*********************************************************************************************
      subroutine root4(b,c,d,e,r1,r2,r3,r4,reals)
!********************************************************************************************* 
      !* PURPOSE:  This subroutine aim on solving quartic equations: x^4+b*x^3+c*x^2+d*x+e=0.
      !* INPUTS:    b, c, d, e-----they are the coefficients of equation. 
      !* OUTPUTS:   r1,r2,r3,r4----roots of the equation, with complex number form.
      !*           reals------------the number of real roots among r1,r2,r3,r4.  
      !* ROUTINES CALLED:  root3   
      !* AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
      !* DATE WRITTEN:  1 Jan 2012 
      !********************************************************************
      implicit none
      Double precision b,c,d,e,q,r,s,two
      parameter(two=2.D0)
      complex*16 r1,r2,r3,r4,s1,s2,s3,temp(1:4),temp1
      integer i,j,del,reals
                         
      reals=0
      q=c-3.D0*b**2/8.D0
      r=d-b*c/two+b**3/8.D0
      s=e-b*d/4.D0+b**2*c/16.D0-3.D0*b**4/256.D0
      call root3(two*q,q**2-4.D0*s,-r**2,s1,s2,s3,del)

      If(del.eq.3)then
          If(real(s3).ge.0.D0)then
              reals=4
              s1=dcmplx(real(sqrt(s1)),0.D0)                
              s2=dcmplx(real(sqrt(s2)),0.D0)                
              s3=dcmplx(real(sqrt(s3)),0.D0)
          else
              reals=0
              s1=sqrt(s1)                
              s2=sqrt(s2)        
              s3=sqrt(s3)
          endif
      else
          If(real(s1).ge.0.D0)then
              reals=2
              s1=dcmplx(real(sqrt(s1)),0.D0)
              s2=sqrt(s2)
              s3=dcmplx(real(s2),-aimag(s2))                 
          else
              reals=0
              s1=sqrt(s1)                
              s2=sqrt(s2)        
              s3=sqrt(s3)
          endif 
      endif 

      if(real(s1*s2*s3)*(-r) .lt. 0.D0)then
         s1=-s1
      end if
      temp(1)=(s1+s2+s3)/two-b/4.D0
      temp(2)=(s1-s2-s3)/two-b/4.D0
      temp(3)=(s2-s1-s3)/two-b/4.D0
      temp(4)=(s3-s2-s1)/two-b/4.D0

      Do i=1,4
          Do j=1+i,4
              If(real(temp(i)).gt.real(temp(j)))then
                  temp1=temp(i)
                  temp(i)=temp(j)
                  temp(j)=temp1
              endif                
          enddo                
      enddo
      r1=temp(1)
      r2=temp(2)
      r3=temp(3)
      r4=temp(4)
      return
      end subroutine root4
!*********************************************************************************************
      subroutine sort(a1,a2,a3,s1,s2,s3)
!********************************************************************
!* PURPOSE:  This subroutine aim on sorting a1, a2, a3 by decreasing way.
!* INPUTS:    a1,a2,a3----they are the number list required to bo sorted. 
!* OUTPUTS:   s1,s2,s3----sorted number list with decreasing way. 
!*      
!* AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!* DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      implicit none
  
      Double precision s1,s2,s3,temp,arr(1:3),a1,a2,a3
      integer i,j
   
      arr(1)=a1
      arr(2)=a2
      arr(3)=a3
   
      Do i=1,3
          Do j=i+1,3
              If(arr(i)<arr(j))then
                  temp=arr(i)
                  arr(i)=arr(j)
                  arr(j)=temp
              end if     
          end do 
      end do
      s1=arr(1)
      s2=arr(2)
      s3=arr(3)
      end subroutine sort
      end module rootsfinding
!*******************************************************************************


     module ellfunction
!********************************************************************
!*    PURPOSE:  This module includes supporting functions and subroutines to compute 
!*              Weierstrass' and Jacobi's elliptical integrals and functions by Carlson's 
!*              integral method. Those codes mainly come from Press (2007) and geokerr.f of
!*              Dexter & Agol (2009).   
!*    AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      use constants
      use rootsfinding
      implicit none

      contains
!***************************************************************************
      Double precision function weierstrassP(z,g2,g3,r1,del)
!***************************************************************************
!*    PURPOSE:   to compute Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    z----------the independent variable value.
!*               g_2, g_3---two parameters.
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    weierstrassP----the value of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  sncndn
!*    AUTHOR:     Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      implicit none
      Double precisionz,g2,g3,e1,e2,e3,k2,u,sn,alp,bet,sig,lamb,cn,dn,&
             two,four,zero
      parameter(two=2.D0,four=4.D0,zero=0.D0)                
      complex*16 r1(1:3)        
      integer  del
        
      If(z.eq.zero)then
          weierstrassP=infinity
      else                 
          z=abs(z)     
          if(del.eq.3)then
              e1=real(r1(1))
              e2=real(r1(2))
              e3=real(r1(3))
              k2=(e2-e3)/(e1-e3)
              u=z*sqrt(e1-e3)
              call sncndn(u,1.D0-k2,sn,cn,dn)
              weierstrassP=e3+(e1-e3)/sn**2
          else
              alp=-real(r1(1))/two
              bet=abs(aimag(r1(2)))
              sig=(9.D0*alp**2+bet**2)**(one/four)
              lamb=(sig**2-three*alp)/bet
              k2=0.5D0+1.5D0*alp/sig**2
              u=two*sig*z
              call sncndn(u,1.D0-k2,sn,cn,dn)   
              if(cn.gt.zero)then                
                  weierstrassP=two*sig**two*(1.D0+sqrt(1.D0-sn**2))/sn**2-two*alp-sig**2
              else
                  If(abs(sn).gt.1.D-7)then
                      weierstrassP=two*sig**two*(1.D0-sqrt(1.D0-sn**2))/sn**2-two*alp-sig**2  
                  else
                      weierstrassP=-two*alp
                  endif       
              endif       
          end if
      endif        
      return
      end function weierstrassP
!*************************************************************************************************
      Function halfperiodwp(r1,del)
!************************************************************************************************* 
!*    PURPOSE:   to compute the semi period of Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
!*               this function involved are real numbers.   
!*    INPUTS:    g_2, g_3---two parameters. 
!*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
!*               del--------number of real roots among r1(1:3). 
!*    RETURN:    halfperiodwp----the semi period of function \wp(z;g_2,g_3).  
!*    ROUTINES CALLED:  rf
!*    AUTHOR:    Yang, Xiao-lin & Wang, Jian-cheng (2012)
!*    DATE WRITTEN:  1 Jan 2012 
!********************************************************************
      implicit none
      Double precision halfperiodwp,e1,e2,e3,zero,one,two,&
                           three,four,EK,alp,bet,sig,k2
      parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0) 
      complex*16 r1(3)      
        integer  del
 
        if(del.eq.3)then
         e1=real(r1(1))
         e2=real(r1(2))
           e3=real(r1(3))
         k2=(e2-e3)/(e1-e3)
         EK=rf(zero,one-k2,one)
         halfperiodwp=EK/sqrt(e1-e3)
      else
         alp=-real(r1(1))/two
          bet=abs(aimag(r1(2)))
         sig=(9*alp**two+bet**two)**(one/four)
         k2=one/two+three/two*alp/sig**two
         EK=rf(zero,one-k2,one)
         halfperiodwp=EK/sig
      endif
      return
      End Function halfperiodwp
!*************************************************************************************************
!       Function halfperiodwp(g2,g3,r1,del)
! !************************************************************************************************* 
! !*    PURPOSE:   to compute the semi period of Weierstrass' elliptical function \wp(z;g_2,g_3) and all of 
! !*               this function involved are real numbers.   
! !*    INPUTS:    g_2, g_3---two parameters. 
! !*               r1(1:3)----an array which is the roots of equation W(t)=4t^3-g_2t-g_3=0.
! !*               del--------number of real roots among r1(1:3). 
! !*    RETURN:    halfperiodwp----the semi period of function \wp(z;g_2,g_3).  
! !*    ROUTINES CALLED:  rf
! !*    AUTHOR:    Yang, Xiao-lin & Wang, Jian-cheng (2012)
! !*    DATE WRITTEN:  1 Jan 2012 
! !********************************************************************
!       implicit none
!       Double precision halfperiodwp,g2,g3,e1,e2,e3,zero,one,two,&
!                            three,four,EK,alp,bet,sig,k2
!       parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0) 
!       complex*16 r1(3)      
!         integer  del
 
!         if(del.eq.3)then
!          e1=real(r1(1))
!          e2=real(r1(2))
!            e3=real(r1(3))
!          k2=(e2-e3)/(e1-e3)
!          EK=rf(zero,one-k2,one)
!          halfperiodwp=EK/sqrt(e1-e3)
!       else
!          alp=-real(r1(1))/two
!           bet=abs(aimag(r1(2)))
!          sig=(9*alp**two+bet**two)**(one/four)
!          k2=one/two+three/two*alp/sig**two
!          EK=rf(zero,one-k2,one)
!          halfperiodwp=EK/sig
!       endif
!       return
!       End Function halfperiodwp
!*************************************************************************************************
      subroutine sncndn(uu,emmc,sn,cn,dn)
!************************************************************************************************* 
!*    PURPOSE:   Returns the Jacobian elliptic functions sn(u|k^2), cn(u|k^2), 
!*            and dn(u|k^2). Here uu=u, while emmc=1-k^2. 
!*    RETURN:    sn, cn, dn----Jacobian elliptic functions sn(u|k^2), cn(u|k^2), dn(u|k^2). 
!*    AUTHOR:    Press et al. (2007) 
!********************************************************************
      implicit none      
      Double precision uu,emmc,sn,CA
      Double precision ,optional :: cn,dn
      parameter (CA=3.D0-8)
      integer i,ii,l
      Double precisiona,b,c,d,emc,u,em(13),en(13)
      logical bo

      emc=emmc
      u=uu
      if(emc.ne.0.D0)then
          bo=(emc.lt.0.D0)        
          if(bo)then
              d=1.D0-emc   !t'=t*k, u'=k*u, k'^2=1./k^2,  
              emc=-emc/d
              d=sqrt(d)
              u=d*u
          end if 
          a=1.D0
          dn=1.D0
          l1: do i=1,13 
              l=i
              em(i)=a
              emc=sqrt(emc)
              en(i)=emc
              c=0.5D0*(a+emc)
              if(abs(a-emc).le.CA*a) exit
              emc=a*emc
              a=c
          end do l1
          u=c*u
          sn=sin(u)  
          cn=cos(u)                
          if(sn.eq.0.D0)then
              if(bo)then
                  a=dn
                  dn=cn
                  cn=a
                  sn=sn/d
              end if
              return        
          endif 
          a=cn/sn
          c=a*c
          l2: do ii=l,1,-1
              b=em(ii)
              a=c*a
              c=dn*c
              dn=(en(ii)+a)/(b+a)
              a=c/b
          enddo l2
          a=1.D0/sqrt(c**2+1.D0)
          if(sn.lt.0.D0)then
              sn=-a
          else
              sn=a
          endif
              cn=c*sn
          return           
      else
          cn=1.D0/cosh(u)
          dn=cn                   
          sn=tanh(u)
          return
      end if
      end subroutine sncndn
!*************************************************************************************************
      Double precision FUNCTION rf(x,y,z) 
!*************************************************************************************************
!*     PURPOSE: Compute Carlson fundamental integral RF
!*              R_F=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2)
!*     ARGUMENTS: Symmetric arguments x,y,z
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (2007).
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision x,y,z,ERRTOL,TINY1,BIG,THIRD,C1,C2,C3,C4,delta,zero
      parameter (ERRTOL=0.0025D0,TINY1=1.5D-38,BIG=3.D37,THIRD=1.D0/3.D0,C1=1.D0/24.D0,&
                          C2=0.1D0,C3=3.D0/44.D0,C4=1.D0/14.D0,zero=0.D0)
      Double precision alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      !if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY1.or.max(x,y,z).gt.BIG)then
      !      rf=0.D0
      !      return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero
      IF(z.lt.zero)z=zero
      xt=x
      yt=y
      zt=z
      sqrtx=sqrt(xt)
      sqrty=sqrt(yt)
      sqrtz=sqrt(zt)
      alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
      xt=0.25D0*(xt+alamb)
      yt=0.25D0*(yt+alamb)
      zt=0.25D0*(zt+alamb)
      ave=THIRD*(xt+yt+zt)
      delx=(ave-xt)/ave
      dely=(ave-yt)/ave
      delz=(ave-zt)/ave      
      delta=max(abs(delx),abs(dely),abs(delz))
      Do while(delta.gt.ERRTOL)
            sqrtx=sqrt(xt)
            sqrty=sqrt(yt)
            sqrtz=sqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=THIRD*(xt+yt+zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delta=max(abs(delx),abs(dely),abs(delz))
      enddo
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.D0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return      
      end Function rf
!************************************************************************ 
      Double precision FUNCTION rj(x,y,z,p) 
!************************************************************************ 
!*     PURPOSE: Compute Carlson fundamental integral RJ
!*     RJ(x,y,z,p) = 3/2 \int_0^\infty dt
!*                      (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-1/2) (t+p)^(-1)
!*     ARGUMENTS: x,y,z,p
!*     ROUTINES CALLED:  RF, RC.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!*********************************************************************** 
       implicit none
      Double precision x,y,z,p,ERRTOL,TINY1,BIG,C1,C2,C3,C4,C5,C6,C7,C8,delta,zero
      Double precision a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,fac,pt,&
      rcx,rho,sqrtx,sqrty,sqrtz,sum1,tau,xt,yt,zt
      Parameter(ERRTOL=0.0015D0,TINY1=2.5D-13,BIG=9.D11,C1=3.D0/14.D0,&
                     C2=1.D0/3.D0,C3=3.D0/22.D0,C4=3.D0/26.D0,C5=0.75D0*C3,&
                           C6=1.5D0*C4,C7=0.5D0*C2,C8=C3+C3,zero=0.D0)
      
      !if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY1.or.max(x,y,z,abs(p)).gt.BIG)then
      !   rj=0.D0            
      !   write(*,*)'ffsdfa'
      !   return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero
      IF(z.lt.zero)z=zero
      sum1=0.D0
      fac=1.D0
      If(p.gt.0.D0)then
            xt=x
            yt=y
            zt=z
            pt=p
      else
            xt=min(x,y,z)
            zt=max(x,y,z)
            yt=x+y+z-xt-zt
            a=1.D0/(yt-p)
            b=a*(zt-yt)*(yt-xt)
            pt=yt+b
            rho=xt*zt/yt
            tau=p*pt/yt
            rcx=rc(rho,tau)
      endif
            sqrtx=sqrt(xt)
            sqrty=sqrt(yt)
            sqrtz=sqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
            beta=pt*(pt+alamb)**2
            sum1=sum1+fac*rc(alpha,beta)
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            pt=0.25D0*(pt+alamb)
            ave=0.2D0*(xt+yt+zt+pt+pt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delp=(ave-pt)/ave
            delta=max(abs(delx),abs(dely),abs(delz),abs(delp))
      Do while(delta.gt.ERRTOL)
            sqrtx=sqrt(xt)
            sqrty=sqrt(yt)
            sqrtz=sqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
            beta=pt*(pt+alamb)**2
            sum1=sum1+fac*rc(alpha,beta)
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            pt=0.25D0*(pt+alamb)
            ave=0.2D0*(xt+yt+zt+pt+pt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
            delp=(ave-pt)/ave
            delta=max(abs(delx),abs(dely),abs(delz),abs(delp))
      enddo  
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.D0*ec
      ee=eb+2.D0*delp*(ea-ec)
      rj=3.D0*sum1+fac*(1.D0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+&
                     delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-&
      C2*delp*ec)/(ave*sqrt(ave))      
      If(p.le.0.D0)rj=a*(b*rj+3.D0*(rcx-rf(xt,yt,zt)))
      return
      end Function rj
!************************************************************************ 
      FUNCTION rc(x,y)
!************************************************************************
!*     PURPOSE: Compute Carlson degenerate integral RC
!*              R_C(x,y)=1/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1)
!*     ARGUMENTS: x,y
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision rc,x,y,ERRTOL,TINY1,sqrtNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,C3,C4
      PARAMETER(ERRTOL=0.0012D0,TINY1=1.69D-38,sqrtNY=1.3D-19,BIG=3.D37,&
      TNBG=TINY1*BIG,COMP1=2.236D0/sqrtNY,COMP2=TNBG*TNBG/25.D0,&
      THIRD=1.D0/3.D0,C1=0.3D0,C2=1.D0/7.D0,C3=0.375D0,C4=9.D0/22.D0)
      Double precisionalamb,ave,s,w,xt,yt 

      if((x.lt.0.D0).or.(y.eq.0.D0).or.((x+abs(y)).lt.TINY1).or.((x+abs(y)).gt.BIG).or.&
                        ((y.lt.-COMP1).and.(x.gt.0.D0).and.(x.lt.COMP2)))then
          rc=0.D0            
          return
      endif
      if(y.gt.0.D0)then
            xt=x
            yt=y
            w=1.D0
      else
            xt=x-y
            yt=-y
            w=sqrt(x)/sqrt(xt)
      endif
            alamb=2.D0*sqrt(xt)*sqrt(yt)+yt
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            ave=THIRD*(xt+yt+yt)
            s=(yt-ave)/ave
      Do While(abs(s).gt.ERRTOL)
            alamb=2.D0*sqrt(xt)*sqrt(yt)+yt
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            ave=THIRD*(xt+yt+yt)
            s=(yt-ave)/ave
      ENDdo
      rc=w*(1.D0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END FUNCTION rc
!**********************************************************************
      FUNCTION rd(x,y,z)
!**********************************************************************
!*     PURPOSE: Compute Carlson degenerate integral RD
!*              R_D(x,y,z)=3/2 \int_0^\infty dt (t+x)^(-1/2) (t+y)^(-1/2) (t+z)^(-3/2)
!*     ARGUMENTS: x,y,z
!*     ROUTINES CALLED:  None.
!*     ALGORITHM: Due to B.C. Carlson.
!*     ACCURACY:  The parameter ERRTOL sets the desired accuracy.
!*     REMARKS:  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,zero
        PARAMETER (ERRTOL=0.0015D0,TINY=1.D-25,BIG=4.5D21,C1=3.D0/14.D0,C2=1.D0/6.D0,&
                         C3=9.D0/22.D0,C4=3.D0/26.D0,C5=0.25D0*C3,C6=1.5D0*C4,zero=0.D0) 
      Double precision alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt

      !if(min(x,y).lt.0.D0.or.min(x+y,z).lt.TINY.or. max(x,y,z).gt.BIG)then
      !      rd=0.D0      
      !      return
      !endif
      IF(x.lt.zero)x=zero
      IF(y.lt.zero)y=zero      
            xt=x
            yt=y
            zt=z
            sum=0.D0
            fac=1.D0
      
            sqrtx=sqrt(xt)
            sqrty=sqrt(yt)
            sqrtz=sqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            sum=sum+fac/(sqrtz*(zt+alamb))
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=0.2D0*(xt+yt+3.D0*zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
      DO While(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)
            sqrtx=sqrt(xt)
            sqrty=sqrt(yt)
            sqrtz=sqrt(zt)
            alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
            sum=sum+fac/(sqrtz*(zt+alamb))
            fac=0.25D0*fac
            xt=0.25D0*(xt+alamb)
            yt=0.25D0*(yt+alamb)
            zt=0.25D0*(zt+alamb)
            ave=0.2D0*(xt+yt+3.D0*zt)
            delx=(ave-xt)/ave
            dely=(ave-yt)/ave
            delz=(ave-zt)/ave
      End DO
            ea=delx*dely
            eb=delz*delz
            ec=ea-eb
            ed=ea-6.D0*eb
            ee=ed+ec+ec
            rd=3.D0*sum+fac*(1.D0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+&
                                 delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
            return
      END function rd
!******************************************************************* 
      Function EllipticF(t,k2)
!******************************************************************* 
!*     PURPOSE: calculate Legendre's first kind elliptic integral: 
!*              F(t,k2)=\int_0^t dt/sqrt{(1-t^2)*(1-k2*t^2)}.  
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF  
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!***********************************************************************
      implicit none
      Double precision t,k2,EllipticF,x1,y1,z1  
      
      x1=1.D0-t*t
      y1=1.D0-k2*t*t
      z1=1.D0
!Press et al. 2007 (6.12.19)
      EllipticF=t*rf(x1,y1,z1)
      return
      end function EllipticF       
!********************************************************************* 
      Function EllipticE(t,k2)
!********************************************************************* 
!*     PURPOSE: calculate Legendre's second kind elliptic integrals: 
!*              E(t,k2)=\int_0^t sqrt{1-k2*t^2}/sqrt{(1-t^2)}dt.
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RD 
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!********************************************************************** 
      implicit none
      Double precision t,k2,EllipticE,x1,y1,z1 
      
      x1=1.D0-t*T
      y1=1.D0-k2*t*t
      z1=1.D0
!Press et al. 2007 (6.12.20)
      EllipticE=t*rf(x1,y1,z1)-1.D0/3.D0*k2*t**3*rd(x1,y1,z1)
      return
      end function EllipticE
!*********************************************************************** 
      Function EllipticPI(t,n,k2)
!*********************************************************************** 
!*     PURPOSE: calculate Legendre's third kind elliptic integrals: 
!*              PI(t,n,k2)=\int_0^t /(1+nt^2)/sqrt{(1-k2*t^2)(1-t^2)}dt. 
!*     ARGUMENTS: t, k2
!*     ROUTINES CALLED:  RF, RJ
!*     AUTHOR:  Press et al (1992)
!*     DATE WRITTEN:  25 Mar 91.
!*     REVISIONS:
!********************************************************************** 
      implicit none
      Double precision t,k2,EllipticPI,x1,y1,z1,w1,n 
      
      x1=1.D0-t*t
      y1=1.D0-k2*t*t
      z1=1.D0
      w1=1.D0+n*t*t
!Press et al. 2007 (6.12.20)
      EllipticPI=t*rf(x1,y1,z1)-1.D0/3.D0*n*t*t*t*rj(x1,y1,z1,w1)
      return
      end function EllipticPI
!*************************************************************************************************
      subroutine weierstrass_int_J3(y,x,bb,del,a4,b4,p4,rff_p,integ,cases)
!*************************************************************************************************
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b4*t+a4)^(k/2)*(4*t^3-g_2*t-g_3)^(-1/2)dt.
!*              Where integer index k can be 0, -2, -4 and 2. (75) and (76) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.
!*              bb(1:3) -- Roots of equation 4*t^3-g_2*t-g_3=0 solved by routine root3.
!*              del -- Number of real roots in bb(1:3).
!*              p4,rff_p,integ -- p4(1:4) is an array which specifies the value of index k of J_k(h). 
!*                 If p4(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p4(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p4(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p4(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p4(4)=-4, then J_{-4} was computed and sent to integ(4).
!*              cases -- If cases=1, then only J_0 was computed.
!*                       If cases=2, then only J_0 and J_{-2} are computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} are computed.    
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} are computed.            
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  ellcubicreals, ellcubiccomplexs     
!*     ACCURACY:   Machine.
!*     REMARKS: Based on Yang & Wang (2012).
!*     AUTHOR:     Yang & Wang (2012).
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      implicit none
      Double precision y,x,yt,xt,a1,b1,a2,b2,a3,b3,a4,b4,rff_p,integ(4),b,three,two,one,&
                  tempt,f,g,h,a44,b44,sign_h
      parameter  (three=3.0D0,two=2.0D0,one=1.D0)
      integer  del,p4(4),i,cases
      complex*16 bb(1:3)
      logical :: inverse,neg

      xt=x
      yt=y
      a44=a4
      b44=b4
      inverse=.false.
      neg=.false.
      If(abs(xt-yt).le.1.D-9)then
         integ=0.D0
         return      
      endif
      if(yt.gt.xt)then
         tempt=xt
         xt=yt
         yt=tempt 
         inverse=.true.                  
      endif
      b=0.D0
      sign_h=sign(one,b44*xt+a44)
      If(del.eq.3)then
! equation (75) of Yang & Wang (2012).
            a44=sign_h*a44
            b44=sign_h*b44
            a1=-real(bb(1))
            a2=-real(bb(2))
            a3=-real(bb(3))
            b1=1.D0
            b2=1.D0
            b3=1.D0
            call ellcubicreals(p4,a1,b1,a2,b2,a3,b3,a44,b44,yt,xt,rff_p*two,integ,cases)
            if(inverse)then
                  integ(1)=-integ(1)
                  integ(2)=-integ(2)
                  integ(3)=-integ(3)
                  integ(4)=-integ(4)
            endif
            Do i=1,4
                integ(i)=integ(i)/two      
                integ(i)=integ(i)*(sign_h)**(-p4(i)/2)
            Enddo
      else
! equation (76) of Yang & Wang (2012).
            a44=sign_h*a44
            b44=sign_h*b44
            a1=-real(bb(1))
            b1=one
            f=real(bb(2))**2+aimag(bb(2))**2
            g=-two*real(bb(2))
            h=one
               call ellcubiccomplexs(p4,a1,b1,a44,b44,f,g,h,yt,xt,rff_p*two,integ,cases)
            if(inverse)then
                  integ(1)=-integ(1)
                  integ(2)=-integ(2)
                  integ(3)=-integ(3)
                  integ(4)=-integ(4)
            endif
            Do i=1,4
                integ(i)=integ(i)/two      
                integ(i)=integ(i)*(sign_h)**(-p4(i)/2)
            Enddo
      endif 
      return
      end subroutine weierstrass_int_J3
!**********************************************************************************************
      subroutine carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,p5,rff_p,integ,cases)
!**********************************************************************************************    
!*     PURPOSE: Computes integrals: J_k(h)=\int^x_y (b5*r+a5)^(k/2)*[(h1*r^2+g1*r+f1)(h2*r^2+g2*r+f2)]^(-1/2)dr.
!*              Where integer index k can be 0, -2, -4, 2 and 4. (77) of Yang & Wang (2012).   
!*     INPUTS:  x,y -- limits of integral.  
!*              p5,rff_p,integ -- p5(1:5) is an array which specifies the value of index k of J_k(h). 
!*                 If p5(1)=0, then J_0 was computed and sent to integ(1).
!*                 If p5(1)!=0, then J_0 was replaced by parameter p, and rff_p=p.  
!*                 If p5(2)=-2, then J_{-2} was computed and sent to integ(2).                      
!*                 If p5(3)=2, then J_{2} was computed and sent to integ(3).
!*                 If p5(4)=-4, then J_{-4} was computed and sent to integ(4).
!*                 If p5(4)=4, then J_{4} was computed and sent to integ(5).
!*              cases -- If cases=1, then only J_0 will be computed.
!*                       If cases=2, then only J_0 and J_{-2} will be computed.
!*                       If cases=3, then only J_0, J_{-2} and J_{2} will be computed.
!*                       If cases=4, then J_0, J_{-2}, J_{2} and J_{-4} will be computed. 
!*                       If cases=5, then J_0, J_{-2}, J_{2} and J_{4} will be computed.     
!*     OUTPUTS: integ -- is an array saved the results of J_k(h). 
!*     ROUTINES CALLED:  elldoublecomplexs
!*     ACCURACY:   Machine.
!*     REMARKS: Based on Yang & Wang (2012).
!*     AUTHOR:     Yang & Wang (2012).
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      implicit none
      Double precision y,x,xt,yt,f1,g1,h1,f2,g2,h2,a5,b5,integ(1:5),tempt,a55,&
                  b55,sign_h,one,zero,rff_p
      parameter(one=1.D0,zero=0.D0)
      integer  p5(1:5),cases,i
      logical :: inverse,neg
      xt=x
      yt=y
      a55=a5
      b55=b5
      inverse=.false.
      neg=.false.
      If(abs(xt-yt).le.1.D-9)then
         integ=zero
         return      
      endif
      if(yt.gt.xt)then
         tempt=xt
         xt=yt
         yt=tempt 
         inverse=.true.                  
      endif
      sign_h=sign(one,b55*xt+a55)      
      a55=sign_h*a55
      b55=sign_h*b55
! equation (77) of Yang & Wang (2012).
         call elldoublecomplexs(p5,f1,g1,h1,f2,g2,h2,a55,b55,yt,xt,rff_p,integ,cases)
      if(inverse)then
            integ(1)=-integ(1)
            integ(2)=-integ(2)
            integ(3)=-integ(3)
            integ(4)=-integ(4)
            integ(5)=-integ(5)
      endif      
      Do i=1,5 
            integ(i)=integ(i)*(sign_h)**(-p5(i)/2)
      Enddo
      return
      end subroutine carlson_doublecomplex5
!*******************************************************************************
      subroutine ellcubicreals(index_p4,a1,b1,a2,b2,a3,b3,a4,b4,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(b2*t+a2)*(b3*t+a3)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has three real roots. 
!*              Equation (61) of Yang & Wang (2012).   
!*     INPUTS:  Arguments for above integral. If index_p4(1)=0, then J_0 will be computed, 
!*              else J_0 will be replaced by parameter p, and rff_p=p. 
!*     OUTPUTS:  Value of integral J_k(h).
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!***********************************************************************
      implicit NONE
      Double precision zero,one,half,two,three,ellcubic,d12,d13,d14,d24,d34,X1,X2,X3,X4,&
                       Y1,Y2,Y3,Y4,U1c,U32,U22,W22,U12,Q22,P22,I1c,I3c,r12,r13,r24i,r34i,&
                       I2c,J2c,K2c,a1,b1,a2,b2,a3,b3,a4,b4,y,x,rff_p,r14,r24,r34,rff,integ(4)
      integer  index_p4(4),cases
      PARAMETER ( ZERO=0.D0, ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0 )
      ellcubic=0.d0
      !c (2.1) Carlson (1989)
      d12=a1*b2-a2*b1
      d13=a1*b3-a3*b1
      d14=a1*b4-a4*b1
      d24=a2*b4-a4*b2
      d34=a3*b4-a4*b3
      r14=a1/b1-a4/b4
      r24=a2/b2-a4/b4
      r34=a3/b3-a4/b4            
      !c (2.2) Carlson (1989)
      X1=dsqrt(abs(a1+b1*x))
      X2=dsqrt(abs(a2+b2*x))
      X3=dsqrt(abs(a3+b3*x))
      X4=dsqrt(abs(a4+b4*x))
      Y1=dsqrt(abs(a1+b1*y))
      Y2=dsqrt(abs(a2+b2*y))
      Y3=dsqrt(abs(a3+b3*y))
      Y4=dsqrt(abs(a4+b4*y))
      !c! (2.3) Carlson (1989)
      If(x.lt.infinity)then      
          U1c=(X1*Y2*Y3+Y1*X2*X3)/(x-y)
          U12=U1c**2
          U22=((X2*Y1*Y3+Y2*X1*X3)/(x-y))**2
          U32=((X3*Y1*Y2+Y3*X1*X2)/(x-y))**2
      else
          U1c=dsqrt(abs(b2*b3))*Y1
          U12=U1c**2
          U22=b1*b3*Y2**2
          U32=b2*b1*Y3**2      
      endif             
      !c (2.4) Carlson (1989)
      W22=U12
      W22=U12-b4*d12*d13/d14
      ! (2.5) Carlson (1989)
      If(x.lt.infinity)then            
          Q22=(X4*Y4/X1/Y1)**2*W22 
      else
          Q22=b4/b1*(Y4/Y1)**2*W22       
      endif 
      P22=Q22+b4*d24*d34/d14
      !c Now, compute the three integrals we need [-1,-1,-1],[-1,-1,-1,-2], and 
      !c  [-1,-1,-1,-4]:we need to calculate the [-1,-1,-1,2] integral,we add it in this part.
      if(index_p4(1).eq.0) then
        !c (2.21) Carlson (1989)
        rff=rf(U32,U22,U12)
        integ(1)=two*rff
        if(cases.eq.1)return
      else
        rff=rff_p/two
      endif
          !c (2.12) Carlson (1989)
          I1c=rff_p!two*rff
      If(index_p4(3).eq.2) then
            !c (2.13) Carlson (1989)
            I2c=two/three*d12*d13*rd(U32,U22,U12)+two*X1*Y1/U1c
            !  (2.39) Carlson (1989)
            integ(3)=(b4*I2c-d14*I1c)/b1      
            if(cases.eq.3)return
      endif
      If(X1*Y1.ne.zero) then
      !c (2.14) Carlson (1989)
          I3c=two*rc(P22,Q22)-two*b1*d12*d13/three/d14*rj(U32,U22,U12,W22)
      Else
       ! One can read the paragraph between (2.19) and (2.20) of Carlson (1989).
        I3c=-two*b1*d12*d13/three/d14*rj(U32,U22,U12,W22)
      Endif
        if(index_p4(2).eq.-2) then
      !c (2.49) Carlson (1989)
          integ(2)=(b4*I3c-b1*I1c)/d14
        if(cases.eq.2)return
        endif

        If(index_p4(4).eq.-4)then
            !c (2.1)  Carlson (1989)
               r12=a1/b1-a2/b2
                r13=a1/b1-a3/b3
                r24i=b2*b4/(a2*b4-a4*b2)
                r34i=b3*b4/(a3*b4-a4*b3)      
             If(x.lt.infinity)then
            !c (2.17) Carlson (1989)
               J2c=two/three*d12*d13*rd(U32,U22,U12)+two*d13*X2*Y2/X3/Y3/U1c
            !c (2.59) & (2.6) Carlson (1989)
               K2c=b3*J2c-two*d34*(X1*X2/X3/X4**2-Y1*Y2/Y3/Y4**2)
           else
               J2c=two/three*d12*d13*rd(U32,U22,U12)+two*d13*Y2/b3/Y3/Y1
               K2c=b3*J2c+two*d34*Y1*Y2/Y3/Y4**2
           endif      
               !c (2.62) Carlson (1989)
               integ(4)=-I3c*half/d14*(one/r14+one/r24+one/r34)+&
               half*b4/d14/d24/d34*K2c+(b1/d14)**2*(one-half*r12*r13*r24i*r34i)*I1c
        endif  
      return
      end  subroutine ellcubicreals
!*******************************************************************************
      subroutine ellcubiccomplexs(index_p4,a1,b1,a4,b4,f,g,h,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (b4*t+a4)^(k/2)[(b1*t+a1)*(h*t^2+g*t+f)]^{-1/2}. 
!*              It is the case of equation W(t)=4*t^3-g_2*t-g_3=0 has one real root. 
!*              Equation (62) of Yang & Wang (2012).  
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p4(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: 
!*********************************************************************** 
      Double precision a1,b1,a4,b4,f,g,h,y,x,X1,X4,Y1,Y4,d14,beta1,beta4,a11,c44,integ(4),&
             a142,xi,eta,M2,Lp2,Lm2,I1c,U,U2,Wp2,W2,Q2,P2,rho,I3c,I2c,r24xr34,r12xr13,&
             N2c,K2c,ellcubic,zero,one,two,four,three,half,six,rff,rdd,r14,Ay1,rff_p
      integer   index_p4(4),cases
      PARAMETER ( ONE=1.D0, TWO=2.D0, HALF=0.5d0, THREE=3.D0, FOUR=4.d0, SIX=6.d0, ZERO=0.D0 )
      ellcubic=0.d0
      X1=dsqrt(abs(a1+b1*x))
      X4=dsqrt(abs(a4+b4*x))      
      Y1=dsqrt(abs(a1+b1*y))      
      Y4=dsqrt(abs(a4+b4*y))      
      r14=a1/b1-a4/b4      
      d14=a1*b4-a4*b1
      !c (2.2) Carlson (1991)
      beta1=g*b1-two*h*a1
      beta4=g*b4-two*h*a4
      !c (2.3) Carlson (1991)
      a11=sqrt(two*f*b1*b1-two*g*a1*b1+two*h*a1*a1)
      c44=sqrt(two*f*b4*b4-two*g*a4*b4+two*h*a4*a4)
      a142=two*f*b1*b4-g*(a1*b4+a4*b1)+two*h*a1*a4
      !c (2.4) Carlson (1991)
      xi=sqrt(f+g*x+h*x*x)
      eta=sqrt(f+g*y+h*y*y)
      !c (3.1) Carlson (1991):
      if(x.lt.infinity)then
         M2=((X1+Y1)*sqrt((xi+eta)**two-h*(x-y)**two)/(x-y))**two
      else 
         M2=b1*(two*sqrt(h)*eta+g+two*h*y)
      endif            
      !c (3.2) Carlson (1991):
      Lp2=M2-beta1+sqrt(two*h)*a11
      Lm2=M2-beta1-sqrt(two*h)*a11

      if(index_p4(1).eq.0) then
      !c (1.2)   Carlson (1991)
        rff=rf(M2,Lm2,Lp2)
        integ(1)=four*rff
        if(cases.eq.1)return
      else
        rff=rff_p/four  
      endif
       !c (3.8)  1991
         I1c=rff_p!four*rff
       !c (3.3) 1991
          if(x.lt.infinity)then      
              U=(X1*eta+Y1*xi)/(x-y)
          else
              U=sqrt(h)*Y1
              !If(Y1.eq.zero)U=1.D-9
          endif 
          
        !c (3.5) 1991
          rho=sqrt(two*h)*a11-beta1
          rdd=rd(M2,Lm2,Lp2)
        If(index_p4(3).eq.2)  then
           !  (3.9) Carlson (1991)
           I2c=a11*sqrt(two/h)/three*(four*rho*rdd-six*rff+three/U)+two*X1*Y1/U
           !  (2.39) Carlson (1989)
           integ(3)=(b4*I2c-d14*I1c)/b1      
           if(cases.eq.3)return            
        endif

         U2=U*U
         Wp2=M2-b1*(a142+a11*c44)/d14
         W2=U2-a11**two*b4/two/d14
       !c (3.4) 1991
          if(x.lt.infinity)then        
                Q2=(X4*Y4/X1/Y1)**two*W2
          else
                Q2=(b4/b1)*(Y4/Y1)**two*W2
         endif       
        P2=Q2+c44**two*b4/two/d14
      !c!! (3.9) 1991
      If(X1*Y1.ne.0.D0) then
        I3c=(two*a11/three/c44)*((-four*b1/d14)*(a142+a11*c44)*rj(M2,Lm2,Lp2,Wp2)&
                        -six*rff+three*rc(U2,W2))+two*rc(P2,Q2)
      Else
        I3c=(two*a11/three/c44)*((-four*b1/d14)*(a142+a11*c44)*rj(M2,Lm2,Lp2,Wp2)&
                        -six*rff+three*rc(U2,W2))
      Endif

        if(index_p4(2).eq.-2) then
        !c (2.49) Carlson (1989)  
          integ(2)=(b4*I3c-b1*I1c)/d14
        if(cases.eq.2)return
        endif

        If(index_p4(4).eq.-4)  then
           ! (2.19) Carlson (1991)
             r24Xr34=half*c44**two/h/b4**two
             r12Xr13=half*a11**two/h/b1**two
           !c (3.11) Carlson (1991)                  
             N2c=two/three*sqrt(two*h)/a11*(four*rho*rdd-six*rff)!+three/U)!+two/X1/Y1/U
           If(Y1.eq.zero)then
             Ay1=-two*b1*xi/X1+sqrt(two*h)*a11/U
             If(x.ge.infinity)Ay1=zero!two*d14*eta/Y4**two/U+sqrt(two*h)*a11/U 
           else
             Ay1=(two*d14*eta/Y4**two+a11**two/X1/U)/Y1+sqrt(two*h)*a11/U
           endif      
           !c (2.5) & (3.12) Carlson (1991)                  
             K2c=half*a11**two*N2c-two*d14*(xi/X1/X4**two)+Ay1 
           !c (2.62) Carlson (1989)
              integ(4)=-I3c*half/d14*(one/r14+two*b4*beta4/c44**two)+&
              half/d14/(h*b4*r24Xr34)*K2c+(b1/d14)**two*(one-half*r12Xr13/r24Xr34)*I1c   
        endif 
      return
      end  subroutine ellcubiccomplexs
!********************************************************************************************
      subroutine  elldoublecomplexs(index_p5,f1,g1,h1,f2,g2,h2,a5,b5,y,x,rff_p,integ,cases)
!*******************************************************************************
!*     PURPOSE: Computes J_k(h)=\int_y^x dt (f_1+g_1t+h_1t^2)^{p_1/2} 
!*                       (f_2+g_2t+h_2t^2)^{p_2/2} (a_5+b_5t)^{p_5/2}. 
!*     INPUTS:  Arguments for above integral.
!*     OUTPUTS:  Value of integral. If index_p5(1)=0, then J_0 will be computed, 
!*               else J_0 will be replaced by parameter p, and rff_p=p. 
!*     ROUTINES CALLED: RF,RJ,RC,RD
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Dexter & Agol (2009)
!*     MODIFIED:   Yang & Wang (2012)
!*     DATE WRITTEN:  4 Mar 2009
!*     REVISIONS: [-1,-1,-1,-1,0]=integ(1),[-1,-1,-1,-1,-2]=integ(2),
!*                [-1,-1,-1,-1,-4]=integ(4),[-1,-1,-1,-1,2]=integ(3),
!*                [-1,-1,-1,-1,4]=integ(5)
!***********************************************************************
      implicit NONE 
      Double precision one,half,two,three,four,six,f1,g1,h1,f2,g2,h2,a5,b5,y,x,xi1,xi2,eta1,eta2,integ(5),&
                       theta1,theta2,zeta1,zeta2,M,M2,delta122,delta112,delta222,delta,deltap,lp2,lm2,&
                       deltam,rff,ellquartic,U,U2,alpha15,beta15,alpha25,beta25,lambda,omega2,psi,xi5,eta5,&
                       gamma1,gamma2,Am111m1,A1111m4,XX,S,mu,T,V2,b2,a2,H,A1111m2,xi1p,B,G,Sigma,&
                       H0,S2,T2,eta1p,psi2,A1111,rff_p
      integer  index_p5(5),cases
      one=1.D0
      half=0.5D0
      two=2.D0
      three=3.D0
      four=4.D0
      six=6.D0
      !c (2.1) Carlson (1992)
      If(x.lt.infinity)then
          xi1=sqrt(f1+g1*x+h1*x**two)
          xi2=sqrt(f2+g2*x+h2*x**two)
      else
          xi1=x*sqrt(h1)
          xi2=x*sqrt(h2)      
      endif            
      eta1=sqrt(f1+g1*y+h1*y*y)
      eta2=sqrt(f2+g2*y+h2*y*y)
      !c (2.4) Carlson (1992)
      If(x.lt.infinity)then
          theta1=two*f1+g1*(x+y)+two*h1*x*y
          theta2=two*f2+g2*(x+y)+two*h2*x*y
      else
          theta1=(g1+two*h1*y)*x
          theta2=(g2+two*h2*y)*x
      endif      
      !c (2.5) Carlson (1992)
      zeta1=sqrt(two*xi1*eta1+theta1)
      zeta2=sqrt(two*xi2*eta2+theta2)
      !c (2.6) Carlson (1992)
      If(x.lt.infinity)then
          M=zeta1*zeta2/(x-y)
          M2=M*M
      else
         M2=(two*sqrt(h1)*eta1+g1+two*h1*y)*(two*sqrt(h2)*eta2+g2+two*h2*y)
      endif

      !c (2.7) Carlson (1992)
      delta122=two*f1*h2+two*f2*h1-g1*g2
      delta112=four*f1*h1-g1*g1
      delta222=four*f2*h2-g2*g2
      Delta=sqrt(delta122*delta122-delta112*delta222)
      !c (2.8) Carlson (1992)
      Deltap=delta122+Delta
      Deltam=delta122-Delta
      Lp2=M2+Deltap
      Lm2=M2+Deltam
 
      if(index_p5(1).eq.0) then
        rff=rf(M2,Lm2,Lp2)
        !c (2.36) Carlson (1992)
        integ(1)=four*rff
        if(cases.eq.1)return
      else
        rff=rff_p/four 
      endif
      !c (2.6) Carlson (1992)
          If(x.lt.infinity)then
             U=(xi1*eta2+eta1*xi2)/(x-y)
             U2=U*U
         else
             U=sqrt(h1)*eta2+sqrt(h2)*eta1
             U2=U*U      
         endif
        
      !c (2.11) Carlson (1992)
        alpha15=two*f1*b5-g1*a5 
        alpha25=two*f2*b5-g2*a5
        beta15=g1*b5-two*h1*a5 
        beta25=g2*b5-two*h2*a5
      !c (2.12) Carlson (1992)
        gamma1=half*(alpha15*b5-beta15*a5)
        gamma2=half*(alpha25*b5-beta25*a5)
      !c (2.13) Carlson (1992)
        Lambda=delta112*gamma2/gamma1
        Omega2=M2+Lambda
        psi=half*(alpha15*beta25-alpha25*beta15)
        psi2=psi*psi
      !c (2.15) Carlson (1992)
        xi5=a5+b5*x
        eta5=a5+b5*y
      !c (2.16) Carlson (1992)
          If(x.lt.infinity)then
              Am111m1=one/xi1*xi2-one/eta1*eta2
              A1111m4=xi1*xi2/xi5**two-eta1*eta2/eta5**two
              A1111m2=xi1*xi2/xi5-eta1*eta2/eta5      
              XX = (xi5*eta5*theta1*half*Am111m1-xi1*xi2*eta5**two+&
                                  eta1*eta2*xi5**two)/(x-y)**two
              mu=gamma1*xi5*eta5/xi1/eta1
          else
              Am111m1=sqrt(h2/h1)-one/eta1*eta2
              A1111m4=sqrt(h1*h2)/b5**two-eta1*eta2/eta5**two      
              A1111m2=sqrt(h1*h2)*x/b5-eta1*eta2/eta5
              XX=b5*eta5*h1*y*Am111m1-eta5**two*sqrt(h1*h2)+eta1*eta2*b5**two
              mu=gamma1*b5*eta5/sqrt(h1)/eta1
         endif

      !c (2.17) Carlson (1992)
    
      !c (2.18) Carlson (1992)
        S=half*(M2+delta122)-U2
        S2=S*S
      !c (2.19) Carlson (1992)
        T=mu*S+two*gamma1*gamma2
        T2=T*T
        V2=mu**two*(S2+Lambda*U2)
      !c (2.20) Carlson (1992)
        b2=Omega2**two*(S2/U2+Lambda)
        a2=b2+Lambda**two*psi2/gamma1/gamma2
      !c (2.22) Carlson (1992)
        H=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+half*rc(a2,b2))/gamma1**two-XX*rc(T2,V2)
      If(index_p5(3).eq.2 .or. index_p5(5).eq.4)then
            !(2.23)--(2.29) Carlson (1992)
            psi=g1*h2-g2*h1
            Lambda=delta112*h2/h1
            Omega2=M2+Lambda
            Am111m1=one/xi1*xi2-one/eta1*eta2
            A1111=xi1*xi2-eta1*eta2
            XX=(theta1*half*Am111m1-A1111)/(x-y)**two
            b2=Omega2**two*(S2/U2+Lambda)
            a2=b2+Lambda**two*psi**two/h1/h2
            mu=h1/xi1/eta1
            T=mu*S+two*h1*h2
            T2=T*T
            V2=mu**two*(S2+Lambda*U2)
            H0=delta112*psi*(rj(M2,Lm2,Lp2,Omega2)/three+&
                          half*rc(a2,b2))/h1**two-XX*rc(T2,V2)
          If(index_p5(3).eq.2)then
            !(2.42) Carlson (1992)
            integ(3)=two*b5*H0-two*beta15*rff/h1
            if(cases.eq.3)return
          endif
          If(index_p5(5).eq.4)then      
              !c (2.2) Carlson (1992)
                If(x.lt.infinity)then
                 xi1p=half*(g1+two*h1*x)/xi1
                else
                 xi1p=sqrt(h1)
                endif
                eta1p=half*(g1+two*h1*y)/eta1
              !c (2.3) Carlson (1992)
                B=xi1p*xi2-eta1p*eta2
              !c (2.9) Carlson (1992)
                If(x.lt.infinity)then
                   G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                       (delta122*theta1-delta112*theta2)/four/xi1/eta1/U  
                else
                   G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                  (delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four/sqrt(h1)/eta1/U  
                endif
              !c (2.10) Carlson (1992)  
                  Sigma=G-Deltap*rff+B
              !(2.44) Carlson (1992)
                  integ(5)=-b5*(beta15/h1+beta25/h2)*H0+b5**two*&
                            Sigma/h1/h2+beta15**two*rff/h1**two
               if(cases.eq.5)return
          endif             
      endif
            if (index_p5(2).eq.-2) then
              !c (2.39) Carlson (1992)
                integ(2)=-two*(b5*H+beta15*rff/gamma1)
            if(cases.eq.2)return
            endif
          If(index_p5(4).eq.-4)then
              !c (2.2) Carlson (1992)
                If(x.lt.infinity)then
                 xi1p=half*(g1+two*h1*x)/xi1
                else
                 xi1p=sqrt(h1)
                endif
                eta1p=half*(g1+two*h1*y)/eta1
              !c (2.3) Carlson (1992)
                B=xi1p*xi2-eta1p*eta2
              !c (2.9) Carlson (1992)
                If(x.lt.infinity)then
                G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                         (delta122*theta1-delta112*theta2)/four/xi1/eta1/U  
                else
                 G=two/three*Delta*Deltap*rd(M2,Lm2,Lp2)+half*Delta/U+&
                  (delta122*(g1+two*h1*y)-delta112*(g2+two*h2*y))/four/sqrt(h1)/eta1/U  
                endif
              !c (2.10) Carlson (1992)  
                Sigma=G-Deltap*rff+B
              !c (2.41) Carlson (1992)
                integ(4)=b5*(beta15/gamma1+beta25/gamma2)*H+beta15**two*rff/gamma1**two+&
                                   b5**two*(Sigma-b5*A1111m2)/gamma1/gamma2
          endif
      integ=ellquartic
      return
      end  subroutine  elldoublecomplexs
!*************************************************************************** 
      end module ellfunction



 


!*******************************************************************************
      module blcoordinate
!*******************************************************************************
!*     PURPOSE: This module aims on computing 4 Boyer-Lindquist coordinates (r,\theta,\phi,t)
!*              and affine parameter \sigam.    
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012 
!***********************************************************************
      use constants
      use rootsfinding
      use ellfunction      
      implicit none

      contains
!******************************************************************************************** 
      SUBROUTINE YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        radi,mu,phi,time,sigma) 
!********************************************************************************************
!*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter 
!*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
!*               and \sigma(p). Cf. discussions in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
!*                              momentum of a photon measured under the LNRF frame. This array 
!*                              can be computed by subroutine lambdaq(...), see below.     
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or initialposition of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radi-----------value of function r(p). 
!*               mu-------------value of function \mu(p). 
!*               phi------------value of function \phi(p).
!*               time-----------value of function t(p).
!*               sigma----------value of function \sigma(p).
!*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.            
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        IMPLICIT NONE 
        DOUBLE PRECISION f1234(4),lambda,q,sinobs,muobs,a_spin,robs,scal,&
                zero,one,two,three,four,phi_r,time_r,aff_r,phi_t,time_t,&
                mu_cos,r_coord,radi,mu,time,phi,sigma,p,Rab
        PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
        LOGICAL rotate,err
        INTEGER tm1,tm2,tr1,tr2

        !************************************************************************************
! call integrat_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p) (here is r_coord).
        call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                        scal,phi_r,time_r,aff_r,r_coord,tr1,tr2) 
! call integrat_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p) (here is mu_cos).
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,&
                        scal,phi_t,time_t,mu_cos,tm1,tm2) 
        radi=r_coord
        mu=mu_cos
!time coordinate value, equation (74) of Yang & Wang (2012).
        time=time_r+time_t

!time   write(*,*)time, time_r, time_t
!affine parameter value, equation (74) of Yang & Wang (2012).
!write(*,*)'ynogk=',aff_r,time_t,time_r,time_t
        sigma=aff_r+time_t
!phi coordinate value.
        rotate=.false.
        err=.false.
!write(*,*)'phi2=',p,phi_r,phi_t,tm1,tm2,time_r,time_t,lambda,f1234(3)
        IF(ABS(muobs).NE.ONE)THEN
! equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_t)
            IF(f1234(3).EQ.zero)THEN 
                phi=phi+(tm1+tm2)*PI
            ENDIF
            phi=DMOD(phi,twopi)
            IF(phi.LT.zero)THEN
                phi=phi+twopi
            ENDIF
        ELSE 
! equation (74) of Yang & Wang (2012).
            phi=-(phi_t+phi_r+(tm1+tm2)*PI)

            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            IF(phi.NE.zero)THEN
                rotate=.TRUE.
            ENDIF
            IF(Rab.NE.zero)THEN
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=muobs*phi+asin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=muobs*phi+twopi+asin(f1234(2)/Rab)
                endif
            ELSE
                phi=zero
            ENDIF
            IF(rotate)THEN
                phi=Mod(phi,twopi)
                IF(phi.LT.zero)THEN
                   phi=phi+twopi
                ENDIF
            ENDIF
        ENDIF        
        RETURN
      END SUBROUTINE YNOGK

!============================================================================================
      Function mucos(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,scal)
!============================================================================================
!*     PURPOSE:  Computes function \mu(p) defined by equation (32) in Yang & Wang (2012). That is
!*               \mu(p)=b0/(4*\wp(p+PI0;g_2,g_3)-b1)+\mu_tp1. \wp(p+PI0;g_2,g_3) is the Weierstrass'
!*               elliptic function.  
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_2, \theta component of four momentum of a photon measured under a LNRF.
!*               f12343---------p_3, \phi component of four momentum of a photon measured under a LNRF..
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  mucos----------\mu coordinate of photon corresponding to a given p. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision mucos,f12343,f12342,p,sinobs,muobs,a_spin,q,lambda,mu_tp,&
                        zero,b0,b1,b2,b3,g2,g3,tinf,fzero,a4,b4,AA,BB,two,scal,four,&
                        mu_tp2,one,three,integ4(4),rff_p
      Double precision f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1
      complex*16 dd(3)
      integer ::  reals,p4,index_p4(4),del,cases,count_num=1
      logical :: mobseqmtp
      save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1,&
                mu_tp,b0,b1,b2,b3,g2,g3,dd,fzero,count_num,AA,BB,del
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0,three=3.D0)

10    continue
      If(count_num.eq.1)then
          f12343_1=f12343
          f12342_1=f12342
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs
          a_spin_1=a_spin 
          scal_1=scal
   !*****************************************************************************************************
          If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
              mucos=muobs     !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
              return          !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
          endif               !so Theta_mu=q+a^2mu^2-X^2mu^4=-a^2(1-mu^2) 
! spin is zero.
          If(a_spin.eq.zero)then
              if(q.gt.zero)then        
                  AA=sqrt((lambda**two+q)/q)
                  BB=sqrt(q)
                  If(f12342.lt.zero)then
                      mucos=sin(asin(muobs*AA)+p*BB*AA)/AA        
                  else
                      If(f12342.eq.zero)then
                          mucos=cos(p*AA*BB)*muobs
                      else                              
                          mucos=sin(asin(muobs*AA)-p*AA*BB)/AA        
                      endif        
                  endif
              else
                  mucos=muobs
              endif
          else
! Equatorial plane motion.
              If(muobs.eq.zero.and.q.eq.zero)then
                  mucos=zero
                  return                        
              endif
              call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
              a4=zero
              b4=one
              p4=0
! equations (26)-(29) in Yang & Wang (2012).
              b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
              b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
              b2=-four/three*a_spin**2*mu_tp
              b3=-a_spin**2
! equation (31) in Yang & Wang (2012).
              g2=three/four*(b1**2-b0*b2)
              g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)

              call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
              index_p4(1)=0        
              cases=1
! equation (33) in Yang & Wang (2012).
              If(muobs.ne.mu_tp)then        
                  tinf=b0/(four*(muobs-mu_tp))+b1/four
                  call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
                  fzero=integ4(1)        
              else
                  fzero=zero                
              endif
              If(f12342.lt.zero)then
                  fzero=-fzero         
              endif
! equation (32) in Yang & Wang (2012).
              mucos=mu_tp+b0/(four*weierstrassP(p+fzero,g2,g3,dd,del)-b1)
              !write(*,*)'mu=',weierstrassP(p+fzero,g2,g3,dd,del),b0,mu_tp,b1!tinf,infinity,g2,g3,a4,b4,p4,fzero
              ! If muobs eq 0,q eq 0,and mu_tp eq 0,so b0 eq 0,
              ! so mucos eq mu_tp eq 0.
              count_num=count_num+1
          endif
 !**************************************************************************       
      else
          If(f12343.eq.f12343_1.and.f12342.eq.f12342_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
              sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.a_spin_1&
              .and.scal.eq.scal_1)then
              !******************************************************************        
              If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                  mucos=muobs                 !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
                  return                 !so,mu must =+1 or -1 for ever. q=-a^2, X=lambda/sin(theta)=0
              endif
              If(a_spin.eq.zero)then
                  if(q.gt.zero)then        
                      If(f12342.lt.zero)then
                          mucos=sin(asin(muobs*AA)+p*BB*AA)/AA        
                      else
                          If(f12342.eq.zero)then
                              mucos=cos(p*AA*BB)*muobs
                          else                              
                              mucos=sin(asin(muobs*AA)-p*AA*BB)/AA        
                          endif        
                      endif
                  else
                      mucos=muobs
                  endif
              else
                  If(muobs.eq.zero.and.q.eq.zero)then
                      mucos=zero
                      return                        
                  endif     
! equation (32) in Yang & Wang (2012).    
                  mucos=mu_tp+b0/(four*weierstrassP(p+fzero,g2,g3,dd,del)-b1) 
                  !write(*,*)'mu=',mucos,weierstrassP(p+fzero,g2,g3,dd,del),mu_tp,b0,b1        
              endif
              !******************************************************************** 
          else
              count_num=1        
              goto 10        
          endif        
      endif
      return
      end Function mucos

!********************************************************************************************
      Function radius(p,f1234r,lambda,q,a_spin,robs,scal)
!============================================================================================
!*     PURPOSE:  Computes function r(p) defined by equation (41) and (49) in Yang & Wang (2012). That is
!*               r(p)=b0/(4*\wp(p+PIr;g_2,g_3)-b1)+r_tp1. \wp(p+PIr;g_2,g_3) is the Weierstrass'
!*               elliptic function; Or r=r_+, r=r_-. 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f1234r---------p_1, r components of four momentum of a photon measured under a LNRF. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.  
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  radius---------radial coordinate of photon corresponding to a given p. 
!*     ROUTINES CALLED: weierstrass_int_J3, mutp, root3.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision radius,p,a_spin,rhorizon,q,lambda,scal,zero,integ4(4),&
             cc,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,robs,cr,dr,integ04(4),&
             u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,a4,b4,one,two,four,sqt3,&
             integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,t_inf,PI0_total,&
             PI0_inf_obs,PI0_obs_hori,PI01,PI0_total_2,rff_p
      Double precision f1234r_1,lambda_1,q_1,a_spin_1,robs_1,scal_1
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer ::  reals,cases_int,del,index_p4(4),cases,count_num=1
      logical :: robs_eq_rtp,indrhorizon
      save  f1234r_1,lambda_1,q_1,a_spin_1,robs_1,scal_1,r_tp1,r_tp2,reals,&
                robs_eq_rtp,indrhorizon,cases,bb,rhorizon,b0,b1,b2,b3,g2,g3,dd,del,cc,tinf,tp2,&
                thorizon,tinf1,PI0,u,w,v,L1,L2,m2,t_inf,pinf,a4,b4,PI0_total,PI0_inf_obs,PI0_obs_hori,&
                PI0_total_2

 20   continue
      If(count_num.eq.1)then
          f1234r_1=f1234r 
          lambda_1=lambda
          q_1=q 
          a_spin_1=a_spin
          robs_1=robs
          scal_1=scal          
    !*********************************************************************************************          
          rhorizon=one+sqrt(one-a_spin**2)
          a4=zero
          b4=one
          cc=a_spin**2-lambda**2-q
          robs_eq_rtp=.false.
          indrhorizon=.false.
          call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                             reals,robs_eq_rtp,indrhorizon,cases,bb)
          If(reals.ne.0)then
! equations (35)-(38) in Yang & Wang (2012).
              b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
              b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
              b2=four/three*r_tp1
              b3=one
              g2=three/four*(b1**2-b0*b2)
              g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
! equation (39) in Yang & Wang (2012).
              If(robs-r_tp1.ne.zero)then          
                  tinf=b0/four/(robs-r_tp1)+b1/four
              else
                  tinf=infinity
              endif 
              If(rhorizon-r_tp1.ne.zero)then
                  thorizon=b1/four+b0/four/(rhorizon-r_tp1)
              else
                  thorizon=infinity           
              endif
              tp2=b0/four/(r_tp2-r_tp1)+b1/four   
              tinf1=b1/four

              call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)          
              index_p4(1)=0          
              cases_int=1
! equation (42) in Yang & Wang (2012).
              call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)          
              PI0=integ04(1)
              select case(cases)
              case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        call weierstrass_int_j3(tinf1,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)        
                        PI0_total=PI0+integ14(1)
                        If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).                             
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                                  
                        else
                            radius=infinity  !Goto infinity, far away.
                        endif        
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif                
                    endif
                else
                    If(f1234r.lt.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                        PI0_obs_hori=integ04(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif
                    else
                        call weierstrass_int_J3(tinf1,tinf,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        PI0_inf_obs=integ04(1)        
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                             radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif        
                    endif                 
                endif
              case(2)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then
                        PI01=-PI0
                    else
                        PI01=PI0        
                    endif
! equation (41) in Yang & Wang (2012). 
                    radius=r_tp1+b0/(four*weierstrassP(p+PI01,g2,g3,dd,del)-b1)                    
                else        
                    If(f1234r.le.zero)then
                        call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)        
                        PI0_obs_hori = integ14(1) 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                        
                    else
                        call weierstrass_int_J3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)                        
                        call weierstrass_int_J3(tp2,tinf,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        PI0_total_2=integ14(1)+integ4(1)
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif        
                    endif                              
                endif
              end select                    
              If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(p.lt.one/rhorizon-one/robs)then
                            radius=robs/(robs*p+one)
                        else
                            radius=rhorizon                  
                        endif
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif                        
                    endif        
                endif
                If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                    else        
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                        PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                        -Log(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif                        
                    endif                
                endif
              endif        
          else
            u=real(bb(4))
            w=abs(aimag(bb(4)))
            v=abs(aimag(bb(2)))
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012). 
                L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012). 
                thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012). 
                m2=(L1-L2)/L1
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                t_inf=sqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012). 
                pinf=EllipticF(tinf,m2)/w/sqrt(L1)
                call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then
                    PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif                    
                else
                    PI0=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                        radius=w*tan(atan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                else
                    if(p.lt.(PI/two-atan(robs/w))/w)then
                        radius=w*tan(atan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif                
                endif
            endif                        
          endif
          count_num=count_num+1
      else
        If(f1234r.eq.f1234r_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
        a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
     !***************************************************************************************************
        If(reals.ne.0)then        
            index_p4(1)=0        
            cases_int=1
            select case(cases)
            case(1)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        If(p.lt.PI0_total)then 
! equation (41) in Yang & Wang (2012).                                
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                                  
                        else
                            radius=infinity  !Goto infinity, far away.
                        endif                
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif                        
                    endif
                else
                    If(f1234r.lt.zero)then 
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                                
                    else
                        If(p.lt.PI0_inf_obs)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=infinity !Goto infinity, far away.
                        endif        
                    endif                 
                endif
            case(2)
                If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then
                        PI01=-PI0
                    else
                        PI01=PI0
                    endif
! equation (41) in Yang & Wang (2012). 
                    radius=r_tp1+b0/(four*weierstrassP(p+PI01,g2,g3,dd,del)-b1)                    
                else        
                    If(f1234r.le.zero)then         
                        If(p.lt.PI0_obs_hori)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p-PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif                        
                    else
                        If(p.lt.PI0_total_2)then
! equation (41) in Yang & Wang (2012). 
                            radius=r_tp1+b0/(four*weierstrassP(p+PI0,g2,g3,dd,del)-b1)                
                        else
                            radius=rhorizon !Fall into black hole.
                        endif        
                    endif                              
                endif
            end select                    
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(p.lt.one/rhorizon-one/robs)then
                            radius=robs/(robs*p+one)
                        else
                            radius=rhorizon                  
                        endif
                    else
                        If(p.lt.one/robs)then
                            radius=robs/(one-robs*p)
                        else
                            radius=infinity          
                        endif                        
                    endif        
                endif
                If(cc.eq.-27.D0)then
                        sqt3=sqrt(three)        
                    If(f1234r.lt.zero)then
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                        If(p.ne.zero)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=robs!infinity
                        endif
                    else        
                        cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                        dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                        PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                        -Log(one+two/sqt3)/three/sqt3
                        If(p.lt.PI0)then        
                            radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                        else
                            radius=infinity
                        endif                        
                    endif                
                endif
            endif        
        else
            If(u.ne.zero)then
                call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                If(f1234r.lt.zero)then        
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r <0, r=r_{+}
                        radius=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))
                    else
                        radius=rhorizon
                    endif                    
                else
                    if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012), and p_r >0, r=r_{-}
                        radius=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**two-(L1-one))        
                    else
                        radius=infinity
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    if(p.lt.(atan(robs/w)-atan(rhorizon/w))/w)then
                        radius=w*tan(atan(robs/w)-p*w)        
                    else
                        radius=rhorizon
                    endif
                else
                    if(p.lt.(PI/two-atan(robs/w))/w)then
                        radius=w*tan(atan(robs/w)+p*w)        
                    else
                        radius=infinity
                    endif                
                endif
            endif                        
        endif
      !***************************************************************************************************
        else
            count_num=1        
            goto  20
        endif
      endif        
      return                 
      End function radius

!********************************************************************************************
      Function phi(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
!********************************************************************************************
!*     PURPOSE:  Computes function \phi(p). 
!*     INPUTS:   p--------------independent variable, which must be nonnegative.  
!*               f1234----------array of p_1, p_2, p_3, p_4, which are the components of four-
!*                              momentum of a photon measured under the LNRF frame. This array 
!*                              can be computed by subroutine lambdaq(...), see below.   
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phi------------value of function \phi(p). 
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ******************************************
        Implicit none
        Double precision p,sinobs,muobs,a_spin,phi,phi_r,phi_c,twopi,f1234(4),timet,&
                         two,Rab,robs,scal,zero,one,lambda,q,&
                         time_r,aff_r,mu_cos,r_coord        
        parameter (zero=0.D0,one=1.D0,two=2.D0)
        logical :: rotate,err
        integer  tm1,tm2,tr1,tr2

        twopi=two*PI
        rotate=.false.
        err=.false.
        IF(a_spin.EQ.ZERO.and.(f1234(2).NE.zero.or.muobs.NE.zero))THEN
!When spin a=0, and the motion of photon is not confined in the equatorial plane, then \phi_r = 0.
            phi_r=zero
        ELSE
            call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                                scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
        ENDIF         
! Call integrat_theta_part to evaluate \phi_\mu.
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,phi_c,timet,mu_cos,tm1,tm2)
        !write(*,*)'ss=', p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,scal
        !write(*,*)'phi2=',phi_r,phi_c,tm1,tm2
        If(abs(muobs).ne.one)then
! When the observer is on the equatorial plane, and p_\theta (f1234(2)) = 0, then the photon is
! confined in the equatorial plane.
            If(muobs.eq.zero.and.f1234(2).eq.zero)then
                phi_c=zero
            endif
!Equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_c)        
            If(f1234(3).eq.zero)then
                phi=phi+(tm1+tm2)*PI
            endif        
            phi=dMod(phi,twopi)
            If(phi.lt.zero)then
                phi=phi+twopi
            Endif
        else       
!Equation (74) of Yang & Wang (2012). 
            phi=-(phi_c+phi_r+(tm1+tm2)*PI)
            Rab=sqrt(f1234(3)**two+f1234(2)**two)
            If(phi.ne.zero)then
                rotate=.true.
            endif
            If(Rab.ne.zero)then
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=phi+dasin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=phi+PI-dasin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=phi+twopi+dasin(f1234(2)/Rab)
                endif
            else
                phi=zero
            endif
            If(rotate)then
                phi=Mod(phi,twopi)
                If(phi.lt.zero)then
                    phi=phi+twopi
                Endif
            Endif        
        endif
        return
      End Function phi

!********************************************************************************************
     SUBROUTINE GEOKERR(p_int,rp,mup,varble,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,tm1,tm2,radi,mu,time,phi,sigma) 
!********************************************************************************************
!*     PURPOSE:  Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine parameter 
!*               \sigma as functions of parameter p, i.e. functions r(p), \mu(p), \phi(p), t(p)
!*               and \sigma(p). Cf. discussions in Yang & Wang (2012).    
!*     INPUTS:   p_int----------this parameter will be taken as independent variable, if 
!*                              varble='p', which must be nonnegative.
!*               rp-------------this parameter will be taken as independent variable, if 
!*                              varble='r'.
!*               mup------------this parameter will be taken as independent variable, if 
!*                              varble='mu'.
!*               varble---------Tell the routine which parameter to be as independent variable,
!*                              r, mu or p.
!*               f12342---------array of f_1, f_2, f_3, f_4, which was defined by equation (102)-(105) 
!*                              in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or initialposition of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. If varble='mu', these two parameter must be provided. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively. If varble='r', these two parameter must be provided.        
!*     OUTPUTS:  radi-----------value of function r(p). 
!*               mu-------------value of function \mu(p). 
!*               phi------------value of function \phi(p).
!*               time-----------value of function t(p).
!*               sigma----------value of function \sigma(p).          
!*               tm1,tm2--------number of times of the photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively. 
!*               tr1,tr2--------number of times of the photon meets turning points r_tp1 and r_tp2
!*                              respectively.
!*     ROUTINES CALLED: INTRPART, INTTPART.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        IMPLICIT NONE 
        DOUBLE PRECISION f1234(4),lambda,q,sinobs,muobs,a_spin,robs,scal,radi,mu,time,phi,sigma,&
                zero,one,two,three,four,phi_r,time_r,aff_r,phi_t,time_t,p,Rab,&
                rp,mup,p_int,mu_cos,r_coord
        CHARACTER varble
        PARAMETER(zero=0.D0, one=1.D0, two=2.D0, three=3.D0, four=4.D0)
        LOGICAL rotate,err
        INTEGER tm1,tm2,tr1,tr2
 
        SELECT CASE(varble)
        CASE('r')
            radi=rp
            p=r2p(f1234(1),rp,lambda,q,a_spin,robs,scal,tr1,tr2)
            mu=mucos(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal) 
        CASE('mu')
            mu=mup
            p=mu2p(f1234(3),f1234(2),lambda,q,mup,sinobs,muobs,a_spin,tm1,tm2,scal)
            radi=radius(p,f1234(1),lambda,q,a_spin,robs,scal)
        CASE('p')
            p=p_int 
        END SELECT 
        !************************************************************************************

! Call integrate_r_part to evaluate t_r,\phi_r,\sigma_r, and function r(p)=r_coord.
        call INTRPART(p,f1234(1),f1234(2),lambda,q,sinobs,muobs,a_spin,robs,&
                              scal,phi_r,time_r,aff_r,r_coord,tr1,tr2)
! Call integrate_theta_part to evaluate t_\mu,\phi_\mu,\sigma_\mu, and function \mu(p)=mu_cos.
        call INTTPART(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal,&
                              phi_t,time_t,mu_cos,tm1,tm2) 
        radi=r_coord
        mu=mu_cos
!time coordinate value **************************************************************
        time=time_r+time_t
!affine parameter value *************************************************************
        sigma=aff_r+time_t  
!phi coordinate value ***************************************************************
        rotate=.false.
        err=.false.
        IF(ABS(muobs).NE.ONE)THEN
! Equation (74) of Yang & Wang (2012).
            phi=-(phi_r+phi_t)
            IF(f1234(3).EQ.zero)THEN 
                phi=phi+(tm1+tm2)*PI
            ENDIF
            phi=DMOD(phi,twopi)
            IF(phi.LT.zero)THEN
                phi=phi+twopi
            ENDIF
        ELSE 
! Equation (74) of Yang & Wang (2012).
            phi=-(phi_t+phi_r+(tm1+tm2)*PI)
            Rab=dsqrt(f1234(3)**two+f1234(2)**two)
            IF(phi.NE.zero)THEN
                rotate=.TRUE.
            ENDIF
            IF(Rab.NE.zero)THEN
! a muobs was multiplied to control the rotate direction
                if((f1234(3).ge.zero).and.(f1234(2).gt.zero))then
                    phi=muobs*phi+asin(f1234(2)/Rab)  
                endif
                if((f1234(3).lt.zero).and.(f1234(2).ge.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).le.zero).and.(f1234(2).lt.zero))then
                    phi=muobs*phi+PI-asin(f1234(2)/Rab)
                endif
                if((f1234(3).gt.zero).and.(f1234(2).le.zero))then
                    phi=muobs*phi+twopi+asin(f1234(2)/Rab)
                endif
            ELSE
                phi=zero
            ENDIF
            IF(rotate)THEN
                phi=Mod(phi,twopi)
                IF(phi.LT.zero)THEN
                    phi=phi+twopi
                ENDIF
            ENDIF
        ENDIF      
        RETURN
      END SUBROUTINE GEOKERR

!**************************************************
      Function rms(a_spin)                       
!**************************************************
!*     PURPOSE: Computes inner most stable circular orbit r_{ms}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of inner most stable circular orbit: r_{ms}
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ******************************************************* 
      implicit none
      Double precision rms,a_spin,b,c,d,e
      complex*16 rt(1:4)
      integer  reals,i
      If(a_spin.eq.0.D0)then
            rms=6.D0
            return
      endif
      b=0.D0
      c=-6.D0
      d=8.D0*a_spin
      e=-3.D0*a_spin**2
        ! Bardeen et al. (1972) 
      call root4(b,c,d,e,rt(1),rt(2),rt(3),rt(4),reals)
      Do i=4,1,-1
         If(aimag(rt(i)).eq.0.D0)then
            rms=real(rt(i))**2
            return
         endif             
      enddo
      end function rms

!*********************************************************** 
      Function rph(a_spin)
!***********************************************************
!*     PURPOSE: Computes photon orbit of circluar orbits: r_{ph}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of photon orbit: r_{ph}
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision rph,a_spin
        ! Bardeen et al. (1972) 
      rph=2.D0*(1.D0+cos(2.D0/3.D0*acos(-a_spin)))
      End function  rph      

!************************************************************* 
      Function rmb(a_spin)
!*************************************************************
!*     PURPOSE: Computes marginally bound orbit of circluar orbits: r_{mb}. 
!*     INPUTS:   a_spin ---- Spin of black hole, on interval [-1,1].
!*     OUTPUTS:  radius of marginally bound orbit: r_{mb}
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision rmb,a_spin
        ! Bardeen et al. (1972)  
      rmb=2.D0-a_spin+2.D0*sqrt(1.D0-a_spin)
      End function  rmb      

!********************************************************************************************
      subroutine mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)
!********************************************************************************************
!*     PURPOSE: Returns the coordinates of turning points \mu_tp1 and \mu_tp2 of poloidal motion, judges
!*                whether the initial poloidal angle \theta_{obs} is one of turning points, if 
!*                it is true then mobseqmtp=.TRUE..  
!*     INPUTS:   f12342---------p_2, the \theta component of four momentum of the photon measured 
!*                              under the LNRF, see equation (84) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  mu_tp1, mu_tp2----the turning points, between which the poloidal motion of 
!*                                 the photon was confined, and mu_tp2 <= mu_tp1. 
!*               reals------number of real roots of equation \Theta_\mu(\mu)=0.
!*               mobseqmtp---If mobseqmtp=.TRUE., then muobs equals to be one of the turning points. 
!*     ROUTINES CALLED: NONE
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
       Double precision sinobs,muobs,a_spin,lambda,q,zero,one,two,four,&
                  mu_tp1,mu_tp2,delta,mutemp,f12342,f12343 
      integer  reals
      logical :: mobseqmtp
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)
 
      mobseqmtp=.false.
        If(a_spin .eq. zero)then 
          If(f12342.ne.zero)then
              mu_tp1=sqrt(q/(lambda**two+q))
               mu_tp2=-mu_tp1      
          else
               mu_tp1=abs(muobs)
               mu_tp2=-mu_tp1
               mobseqmtp=.true.
          endif
          reals=2    
        ELSE
          If(lambda.ne.zero)then
            delta=(a_spin**two-lambda**two-q)**two+four*a_spin**two*q
            mu_tp1=dsqrt( dabs((dsqrt(delta)-(lambda**two+q-a_spin**two))/two) )/dabs(a_spin) 
            If(dsqrt(delta)+(lambda**two+q-a_spin**two).le.zero)then
                mu_tp2=dsqrt(-(dsqrt(delta)+(lambda**two+q-a_spin**two))/two)/dabs(a_spin)
                If(f12342.eq.zero)then
                  If(abs(muobs-mu_tp1).le.1.D-4)then
                      mu_tp1=dabs(muobs)
                  else       
                      mu_tp2=dabs(muobs)
                  endif
                  mobseqmtp=.true.
                endif      
                reals=4
            else
                If(f12342.ne.zero)then      
                  mu_tp2=-mu_tp1
                else
                  mu_tp1=dabs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
                endif      
                reals=2
            endif      
          else 
            If(abs(muobs).ne.one)then
                If(q.le.zero)then
                    If(f12342.ne.zero)then
                        mu_tp2=dsqrt(-q)/dabs(a_spin)
                    else
                        mu_tp2=dabs(muobs)!a=B=zero.
                        mobseqmtp=.true.
                    endif
                    reals=4
                else
                    mu_tp2=-one
                    reals=2
                endif 
                mu_tp1=one
            else
                mu_tp1=one
                If(q.le.zero.and.f12342*f12342+f12343*f12343.ne.zero)then
                    mu_tp2=dsqrt(-q)/dabs(a_spin)
                    reals=4
                else
                    mu_tp2=-one
                    reals=2
                endif
            endif
          endif
        ENDIF
      If(abs(muobs).eq.one)mobseqmtp=.true.
      If(muobs.lt.zero.and.reals.eq.4)then
            mutemp=mu_tp1
            mu_tp1=-mu_tp2
            mu_tp2=-mutemp
      endif
      return
      end subroutine mutp      

!============================================================================================
      Subroutine radiustp(f12341,a_spin,robs,lambda,q,r_tp1,&
                        r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
!********************************************************************************************
!*     PURPOSE: Returns the coordinates of turning points r_tp1 and r_tp2 of radial motion, judges
!*                whether the initial radius robs is one of turning points, if 
!*                it is true then robs_eq_rtp=.TRUE.. And if r_tp1 less or equal r_horizon,
!*                then indrhorizon=.TRUE. Where r_horizon is the radius of the event horizon.  
!*     INPUTS:   f12341---------p_r, the r component of four momentum of the photon measured 
!*                              under the LNRF, see equation (83) in Yang & Wang (2012).
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*     OUTPUTS:  r_tp1, r_tp2----the turning points, between which the radial motion of 
!*                                 the photon was confined, and r_tp2 >= r_tp1.
!*               bb(1:4)----roots of equation R(r)=0.                
!*               reals------number of real roots of equation R(r)=0.
!*               robs_eq_rtp---If robs_eq_rtp=.TRUE., then robs equal to be one of turning points. 
!*               cases-------If r_tp2=infinity, then cases=1, else cases=2.
!*               indrhorizon----if r_tp1 less or equals r_horizon, indrhorizon=.TRUE.. 
!*     ROUTINES CALLED: root4
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
       Double precision f12341,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                  zero,one,two,four,b1,c1,d1,e1,r1(2),rhorizon
      integer  reals,i,j,cases
      logical :: robs_eq_rtp,indrhorizon
      complex*16 bb(4)
      parameter (zero=0.D0,two=2.0D0,four=4.D0,one=1.D0)       

      rhorizon=one+sqrt(one-a_spin**two)
      robs_eq_rtp=.false.
      indrhorizon=.false.
      b1=zero
      c1=a_spin**two-lambda**two-q
      d1=two*(q+(lambda-a_spin)**two)
      e1=-q*a_spin**two
      call root4(b1,c1,d1,e1,bb(1),bb(2),bb(3),bb(4),reals)      
          SELECT CASE(reals) 
          CASE(4)   
              IF(f12341.eq.zero)THEN
                  IF(dabs(robs-real(bb(4))) .LE. 1.D-4)THEN
                      r_tp1=robs
                      r_tp2=infinity
                      cases=1
                  ENDIF 
                  IF(dabs(robs-real(bb(2))) .LE. 1.D-4)THEN
                      r_tp1=robs  
                      r_tp2=real(bb(3))
                      cases=2
                  ENDIF
                  IF(dabs(robs-real(bb(3))) .LE. 1.D-4)THEN
                      r_tp1=real(bb(2))       
                      r_tp2=robs
                      cases=2
                  ENDIF
                  IF(dabs(robs-real(bb(1))) .LE. 1.D-4)THEN
                      r_tp1=-infinity
                      r_tp2=robs 
                      cases=3
                      write(*,*)'radiustp(): wrong! 4 roots, cases = 3'
                      stop  
                  ENDIF  
                  robs_eq_rtp = .TRUE. 
              ELSE 
                  If( robs.ge.real(bb(4)) )then      
                      r_tp1=real(bb(4))      
                      r_tp2=infinity
                      cases=1 
                  else
                      If( (robs.ge.real(bb(2)) .and. robs.le.real(bb(3))) )then      
                          r_tp1=real(bb(2))      
                          r_tp2=real(bb(3))
                          cases=2
                      else
                          IF( real(bb(1)) .GT. rhorizon .AND.  robs .LE. real(bb(1))  )THEN
                              write(*,*)'radiustp(): wrong! 4 roots,cases = 3'
                              stop 
                              r_tp2=real(bb(1))  
                              r_tp1=-infinity 
                          ELSE   
                              write(*,*)'radiustp(): wrong! 4 roots',robs,bb
                              stop
                          ENDIF
                      endif
                  endif
              ENDIF          
          CASE(2)
              j=1
              Do  i=1,4
                  If (aimag(bb(i)).eq.zero) then
                      r1(j)=real(bb(i))
                      j=j+1      
                  endif
              Enddo
              If( robs.ge.r1(2) )then      
                  r_tp1=r1(2)      
                  r_tp2=infinity
                  cases=1
              else  
                  If( r1(1).ge.rhorizon .and. robs.le.r1(1) )then
                      write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                      stop      
                  endif
              endif
              IF(f12341.eq.zero)THEN 
                  IF(dabs(robs-r1(2)) .LE. 1.D-4)THEN
                      r_tp1=robs 
                      r_tp2=infinity
                  ENDIF
                  IF(dabs(robs-r1(1)) .LE. 1.D-4)THEN
                      r_tp1=-infinity 
                      r_tp2=robs 
                      write(*,*)'radiustp(): wrong! 2 roots, cases = 3'
                      stop
                  ENDIF
                  robs_eq_rtp=.TRUE. 
              ENDIF
          CASE(0)
              r_tp1=zero
              r_tp2=infinity
              cases=1       
          END SELECT  
 
      IF(rhorizon.ge.r_tp1 .and. rhorizon.le.r_tp2)then
          indrhorizon=.true.
      Endif
      End Subroutine radiustp

!********************************************************************************************  
      Function mu2p(f12343,f12342,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)
!********************************************************************************************
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  value of \mu part of integral of (24). 
!*     ROUTINES CALLED: mu2p_schwartz, mutp, root3, weierstrass_int_J3
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision mu2p,f12342,f12343,mu,sinobs,muobs,a_spin,lambda,q,mu_tp,tposition,tp2,four,&
                   b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,a4,b4,two,mu_tp2,&
                   scal,zero,one,integ4(4),three,rff_p
      parameter (zero=0.D0,two=2.D0,four=4.D0,one=1.D0,three=3.D0)
      integer  t1,t2,reals,p4,index_p4(4),del,cases
      complex*16 dd(3)
      logical :: mobseqmtp  

      If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
          mu2p=zero!-one
          return            
      endif
      If(a_spin.eq.zero)then
            call mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
            return      
      endif
      
      a4=zero
      b4=one
      p4=0
      mobseqmtp=.false.
      call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp)
! equatorial plane motion.
      If(mu_tp.eq.zero)then
            mu2p=zero
            return
      endif
! equations (26)-(29) in Yang & Wang (2012).
      b0=-four*a_spin**2*mu_tp**3+two*mu_tp*(a_spin**2-lambda**2-q)
      b1=-two*a_spin**2*mu_tp**2+one/three*(a_spin**2-lambda**2-q)
      b2=-four/three*a_spin**2*mu_tp
      b3=-a_spin**2
      g2=three/four*(b1**2-b0*b2)
      g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
! equation (30) in Yang & Wang (2012).
      If(abs(mu-mu_tp).ne.zero)then
           tposition=b0/(four*(mu-mu_tp))+b1/four
      else
           tposition=infinity      
      endif
      If(muobs.ne.mu_tp)then      
           tinf=b0/four/(muobs-mu_tp)+b1/four
      else
           tinf=infinity
      endif      

      call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
      index_p4(1)=0
      cases=1 

        If(mu.gt.mu_tp.or.mu.lt.mu_tp2)then
              mu2p=-one
              return
        endif      
! equation (30) in Yang & Wang (2012).
              tp2=b0/four/(mu_tp2-mu_tp)+b1/four
            If(t1.eq.0)then
               p1=zero
            else
! equation (53) in Yang & Wang (2012).
               call weierstrass_int_J3(tposition,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases)
               p1=integ4(1)
            endif
            If(t2.eq.0)then
               p2=zero
            else
               call weierstrass_int_J3(tp2,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)      
               p2=integ4(1)
            endif
               call weierstrass_int_J3(tinf,tposition,dd,del,a4,b4,index_p4,rff_p,integ4,cases)            
               pp=integ4(1) 

! equation (54) in Yang & Wang (2012).
      If(mobseqmtp)then
          If(muobs.eq.mu_tp)then  
              mu2p=-pp+two*(t1*p1+t2*p2)            
          else
              mu2p=pp+two*(t1*p1+t2*p2)            
          endif       
      else
          If(f12342.lt.zero)then
            mu2p=pp+two*(t1*p1+t2*p2)
          endif
          If(f12342.gt.zero)then            
            mu2p=-pp+two*(t1*p1+t2*p2)
          endif      
      endif
      return
      end Function mu2p

!============================================================================================
      subroutine mu2p_schwartz(f12343,f12342,lambda,q,mu,sinobs,muobs,t1,t2,mu2p,scal)
!********************************************************************************************
!*     PURPOSE:  Computes the value of parameter p from \mu coordinate. In other words, to compute 
!*               the \mu part of integral of equation (24), using formula (54) in Yang & Wang (2012).
!*               (54) is: p=-sign(p_\theta)*p_0+2*t1*p_2+2*t2*p_2. where p_\theta is initial \theta
!*               component of 4 momentum of photon.
!*               And black hole spin is zero. 
!*     INPUTS:   f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               t1,t2----------Number of photon meets the turning points \mu_tp1 and \mu_tp2
!*                              respectively.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.     
!*               mu-------------\mu coordinate of photon.     
!*     OUTPUTS:  mu2p-----------value of \mu part of integral of (24). 
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision f12343,f12342,mu,sinobs,muobs,mu2p,pp,p1,p2,BB,two,&
                   lambda,q,scal,zero,one,mu_tp,mu_tp2
      integer  t1,t2      
      parameter(two=2.D0,zero=0.D0,one=1.D0)
      logical :: mobseqmtp
      If(f12343.eq.zero.and.f12342.eq.zero)then !this Theta=q(1-mu^2),so if B=0,then q=0.
          mu2p=-two            !so Theta_mu=0 for ever.But we do not need to 
          return  !consider it,for q=0,so the next part means that it will return
      endif !zero value.
      mobseqmtp=.false.
      If(q.gt.zero)then      
          BB=sqrt(q)
          If(f12342.ne.zero)then
              mu_tp=sqrt(q/(lambda**two+q))
              mu_tp2=-mu_tp      
          else
              mu_tp=muobs
              mu_tp2=-mu_tp
              mobseqmtp=.true.
          endif 
          If(abs(muobs).eq.one)mobseqmtp=.true.            
          pp=(asin(mu/mu_tp)-asin(muobs/mu_tp))*mu_tp/BB      
          If(t1.eq.0)then
              p1=zero
          else      
              p1=(PI/two-asin(mu/mu_tp))*mu_tp/BB      
          endif
          If(t2.eq.0)then
              p2=zero
          else      
              p2=(asin(mu/mu_tp)+PI/two)*mu_tp/BB
            endif
          If(mobseqmtp)then
              If(muobs.eq.mu_tp)then  
                  mu2p=-pp+two*(t1*p1+t2*p2)            
              else
                  mu2p=pp+two*(t1*p1+t2*p2)            
              endif             
          else
              mu2p=sign(one,-f12342)*pp+two*(t1*p1+t2*p2)       
          endif
      else      
          mu2p=zero
      endif            
      return
      end subroutine mu2p_schwartz
!********************************************************************************************
      Function r2p(f1234r,rend,lambda,q,a_spin,robs,scal,t1,t2)
!============================================================================================
!*     PURPOSE:  Computes the value of parameter p from radial coordinate. In other words, to compute 
!*               the r part of integral of equation (24), using formula (58) in Yang & Wang (2012).
!*               (58) is: p=-sign(p_r)*p_0+2*t1*p_2+2*t2*p_2. where p_r is initial radial
!*               component of 4 momentum of photon. 
!*     INPUTS:   f1234r---------f_1, which was defined by equation (106) in Yang & Wang (2012). 
!*               a_spin---------spin of black hole, on interval (-1,1).
!*               robs-----------radial coordinate of the observer. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.
!*               scal-----------a dimentionless parameter to control the size of the images. 
!*                              Which is usually be set to 1.D0.   
!*               t1,t2----------Number of photon meets the turning points r_tp1 and r_tp2
!*                              respectively in radial motion.
!*     OUTPUTS:  r2p------------value of r part of integral (24) in Yang & Wang (2012).
!*     ROUTINES CALLED: radiustp, root3, weierstrass_int_j3, EllipticF.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  4 Jan 2012
!*     REVISIONS: ****************************************** 
      implicit none
      Double precision r2p,a_spin,rhorizon,q,lambda,scal,zero,integ4(4),&
                    cc,b0,b1,b2,b3,g2,g3,tinf,tinf1,PI0,robs,integ04(4),&
                    u,v,w,L1,L2,thorizon,m2,pinf,a4,b4,one,two,four,sqrt3,&
                    integ14(4),three,six,nine,r_tp1,r_tp2,f1234r,tp2,tp,t_inf,&
                    pp,p1,p2,rend,rff_p
      parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0,six=6.D0,nine=9.D0)
      complex*16 bb(1:4),dd(3)
      integer  reals,cases_int,del,index_p4(4),cases,t1,t2
      logical :: robs_eq_rtp,indrhorizon
      
      rhorizon=one+sqrt(one-a_spin**2)
      a4=zero
      b4=one
      cc=a_spin**2-lambda**2-q
      robs_eq_rtp=.false.
      indrhorizon=.false.
      call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,&
                    r_tp2,reals,robs_eq_rtp,indrhorizon,cases,bb)
      If(reals.ne.0)then
          If(rend.lt.r_tp1.or.rend.gt.r_tp2)then
              r2p=-one
              return
          endif
! equations (35)-(38) in Yang & Wang (2012).
          b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
          b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
          b2=four/three*r_tp1
          b3=one
          g2=three/four*(b1**2-b0*b2)
          g3=one/16.D0*(3*b0*b1*b2-2*b1**3-b0**2*b3)
! equation (39) in Yang & Wang (2012).
          If(robs-r_tp1.ne.zero)then      
              tinf=b0/four/(robs-r_tp1)+b1/four
          else
              tinf=infinity
          endif 
          If(rhorizon-r_tp1.ne.zero)then
              thorizon=b1/four+b0/four/(rhorizon-r_tp1)
          else
              thorizon=infinity       
          endif
          If(rend-r_tp1.ne.zero)then
              tp=b1/four+b0/four/(rend-r_tp1)
          else
              tp=infinity       
          endif
          tp2=b0/four/(r_tp2-r_tp1)+b1/four   
          tinf1=b1/four

          call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)      
          index_p4(1)=0      
          cases_int=1
! equation (42) in Yang & Wang (2012).
          call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)      
          PI0=integ04(1)
          select case(cases)
          case(1)
              If(.not.indrhorizon)then
                    If(f1234r.lt.zero)then  
                        call weierstrass_int_j3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)      
                        pp=integ14(1)
                        If(t1.ne.zero)then     
! equation (57) in Yang & Wang (2012).                   
                            call weierstrass_int_j3(tp,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            p1=integ14(1)              
                        else
                            p1=zero  !Goto infinity, far away.
                        endif 
! equation (58) in Yang & Wang (2012).
                        r2p=pp+two*p1*t1
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                        pp=integ04(1)
                        r2p=-pp        
                    endif
                else
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            tp=thorizon
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                            r2p=integ04(1)
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                            r2p=integ04(1)        
                        endif
                    else
                        If(rend.lt.infinity)then
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                            r2p=-pp
                        else
                            call weierstrass_int_J3(tinf,tinf1,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)        
                            r2p=-pp        
                        endif
                    endif                 
                endif
            case(2)
                If(.not.indrhorizon)then   
! equation (57) in Yang & Wang (2012).                 
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t1.eq.zero)then
                            p1=zero
                        else
                            call weierstrass_int_J3(tp,infinity,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            p1=integ4(1)        
                        endif        
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)        
                        endif
                        If(f1234r.ne.zero)then
                            r2p=sign(one,-f1234r)*pp+two*(t1*p1+t2*p2)
                        else
! equation (58) in Yang & Wang (2012).
                            If(robs.eq.r_tp1)then
                                r2p=-pp+two*(t1*p1+t2*p2)
                            else
                                r2p=pp+two*(t1*p1+t2*p2)
                            endif        
                        endif                    
                else        
                    If(f1234r.le.zero)then
                        If(rend.le.rhorizon)then
                            call weierstrass_int_J3(tinf,thorizon,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)                            
                        else
                            call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            pp=integ4(1)
                        endif
                    else
                        call weierstrass_int_J3(tinf,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                        pp=integ4(1)
                        If(t2.eq.zero)then
                            p2=zero
                        else
                            call weierstrass_int_J3(tp2,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            p2=integ4(1)
                        endif
! equation (58) in Yang & Wang (2012).
                        r2p=-pp+two*t2*p2
                    endif                              
                endif
            end select                    
            If(a_spin.eq.zero)then
                If(cc.eq.zero)then
                    If(f1234r.lt.zero)then
                        If(rend.le.rhorizon)then
                            r2p=one/rhorizon-one/robs
                        else
                            r2p=one/rend-one/robs         
                        endif
                    else
                        If(rend.lt.infinity)then
                            r2p=one/robs-one/rend
                        else
                            r2p=one/robs
                        endif 
                    endif 
                endif
                If(cc.eq.-27.D0)then        
                    sqrt3=sqrt(three)        
                    If(f1234r.lt.zero)then
                        If(rend.gt.rhorizon)then
                            r2p=-log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/(sqrt3))/&
                                                         (robs-three)))/(three*sqrt3)+&
                                 log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/&
                                                         (sqrt3))/(rend-three)))/(three*sqrt3)
                        else
                            r2p=-log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/&
                                            (sqrt3))/(robs-three)))/(three*sqrt3)+&
                                 log(abs((sqrt(rhorizon*(rhorizon+6.D0))+(three+two*rhorizon)&
                                            /(sqrt3))/(rhorizon-three)))/(three*sqrt3)
                        endif
                    else        
                        If(rend.lt.infinity)then
                            r2p=-log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/(sqrt3))/&
                                    (rend-three)))/(three*sqrt3)+log(abs((sqrt(robs*(robs+6.D0))+&
                                    (three+two*robs)/(sqrt3))/(robs-three)))/(three*sqrt3)                           
                        else
                            r2p=-log(one+two/sqrt3)/three/sqrt3+&
                                 log(abs((sqrt(rend*(rend+6.D0))+(three+two*rend)/&
                                    (sqrt3))/(rend-three)))/(three*sqrt3)
                        endif                                                
                    endif                
                endif
            endif        
        else
! equation (44) in Yang & Wang (2012).
            u=real(bb(4))
            w=abs(aimag(bb(4)))
            v=abs(aimag(bb(2)))
            If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).
                L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).
                thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
                tp=sqrt((L1-one)/(L1-L2))*(rend-u*(L1+one)/(L1-one))/sqrt((rend-u)**2+w**2)
! equation (48) in Yang & Wang (2012).
                m2=(L1-L2)/L1
                tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                t_inf=sqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).
                pinf=EllipticF(tinf,m2)/w/sqrt(L1) 
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))
                    else
                        r2p=pinf-EllipticF(tp,m2)/(w*sqrt(L1))
                    endif                            
                else
                    If(rend.lt.infinity)then
                        r2p=EllipticF(tp,m2)/(w*sqrt(L1))-pinf
                    else
                        r2p=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                    endif
                endif
            else
                If(f1234r.lt.zero)then
                    If(rend.le.rhorizon)then
                        r2p=(atan(robs/w)-atan(rhorizon/w))/w
                    else
                        r2p=(atan(robs/w)-atan(rend/w))/w
                    endif        
                else
                    if(rend.lt.infinity)then
                        r2p=(atan(rend/w)-atan(robs/w))/w
                    else
                        r2p=(PI/two-atan(robs/w))/w
                    endif                
                endif
            endif                        
      endif        
      return                 
      End function r2p 

!********************************************************************************************
      SUBROUTINE INTTPART(p,f12343,f12342,lambda,q,sinobs,muobs,a_spin,scal,phyt,timet,mucos,t1,t2)    
!********************************************************************************************
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, t and affine parameter \sigma,
!*               expressed by equation (71) and (72) in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f12342---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f12343---------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyt-----------value of integral \phi_\theta expressed by equation (72) in
!*                              Yang & Wang (2012).  
!*               timet----------value of integral \t_\theta expressed by equation (71) in
!*                              Yang & Wang (2012). And \sigma_\theta=time_\theta.
!*               mucos----------value of function \mu(p).
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.            
!*     ROUTINES CALLED: mutp, root3, weierstrass_int_J3 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      USE constants
      IMPLICIT NONE
      Double precision phyt,timet,f12343,f12342,p,sinobs,muobs,a_spin,lambda,q,mu_tp1,tp2,tmu,&
             b0,b1,b2,b3,g2,g3,tinf,p1,p2,pp,Wmup,Wmum,tplus,tminus,&
             a4,b4,mu_tp2,scal,integ4(4),integ(4),rff_p,&
             integ14(4),PI0,integ04(4),mu2p,PI01,h,p1_t,p2_t,pp_t,p1_phi,&
             p2_phi,pp_phi,mucos,p_mt1_mt2,&
             PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p
      Double precision f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1 
      integer ::  t1,t2,i,j,reals,index_p4(4),del,cases_int,count_num=1
      complex*16 dd(3)
      logical :: mobseqmtp
      save  f12343_1,f12342_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,scal_1,a4,b4,mu_tp1,mu_tp2,reals,&
                mobseqmtp,b0,b1,b2,b3,g2,g3,dd,del,PI0,Wmup,Wmum,tplus,tminus,tp2,tinf,h,p_mt1_mt2,&
                PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,PI01

30      continue
        IF(count_num.eq.1)then        
            f12343_1=f12343
            f12342_1=f12342
            lambda_1=lambda
            q_1=q
            muobs_1=muobs
            sinobs_1=sinobs
            a_spin_1=a_spin 
            scal_1=scal
            t1=0
            t2=0
            !************************************************************************
            If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                mucos = sign(one,muobs)
                timet=zero             !this is because that mu==1 for ever
                phyt=zero              !this is because that mu==1 for ever,this 
                count_num=count_num+1
                return        !because that Theta_mu=-a^2(1-mu^2), so,mu must =+1 or -1 for ever.
            endif
            If(muobs.eq.zero.and.(abs(lambda).lt.abs(a_spin)).and.q.eq.zero)then
                timet=zero
                phyt=zero
                mucos=zero 
                count_num=count_num+1
                return                        
            endif
            mobseqmtp=.false.
            call mutp(f12342,f12343,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)        
            If(mu_tp1.eq.zero)then
            !photons are confined in the equatorial plane, so the integrations about \theta are valished.
                timet=zero
                phyt=zero
                mucos=zero
                count_num=count_num+1
                return
            endif
            !**************************************************************************
            If(a_spin.eq.zero)then
                timet=zero
                CALL phyt_schwatz(p,f12343,f12342,lambda,q,sinobs,muobs,scal,phyt,mucos,t1,t2)
                count_num=count_num+1
                return
            endif
            a4=zero
            b4=one 
! equations (26)-(29) in Yang & Wang (2012). 
            b0=-four*a_spin**2*mu_tp1**3+two*mu_tp1*(a_spin**2-lambda**2-q)
            b1=-two*a_spin**2*mu_tp1**2+one/three*(a_spin**2-lambda**2-q)
            b2=-four/three*a_spin**2*mu_tp1
            b3=-a_spin**2
            g2=three/four*(b1**2-b0*b2)
            g3=one/16.D0*(three*b0*b1*b2-two*b1**3-b0**2*b3)
            call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)
! equation (30) in Yang & Wang (2012). 
            If(muobs.ne.mu_tp1)then        
                tinf=b0/four/(muobs-mu_tp1)+b1/four
            else
                tinf=infinity
            endif
            If(mu_tp1-one.ne.zero)then
! equation (64) in Yang & Wang (2012). 
                Wmum=b0/(eight*(-one-mu_tp1)**2)
                Wmup=b0/(eight*(one-mu_tp1)**2) 
                tminus=b0/four/(-one-mu_tp1)+b1/four
                tplus=b0/four/(one-mu_tp1)+b1/four
            endif
            index_p4(1)=0
            cases_int=1
            call weierstrass_int_J3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
            PI0=integ04(1)     
! equation (34) in Yang & Wang (2012).    
            If(f12342.lt.zero)then
                PI01=-PI0        
            else
                PI01=PI0
            endif
            tmu=weierstrassP(p+PI01,g2,g3,dd,del)
! equation (32) in Yang & Wang (2012). 
            mucos = mu_tp1+b0/(four*tmu-b1)
            tp2=b0/four/(mu_tp2-mu_tp1)+b1/four
            h=-b1/four        
            !to get number of turning points of t1 and t2.
            !111111111********************************************************************************** 
                call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int) 
            !write(*,*)tp2,tinf,tmu,mu_tp1,mu_tp2
! equation (51) in Yang & Wang (2012).        
                p_mt1_mt2=integ14(1)
                PI2_p=p_mt1_mt2-PI0
                pp=integ4(1)
                p1=PI0-pp
                p2=p_mt1_mt2-p1
                PI1_phi=zero
                PI2_phi=zero
                PI1_time=zero
                PI2_time=zero 
                Do j=0,10
                    Do i=j,j+1 
                        If(mobseqmtp)then
                            If(muobs.eq.mu_tp1)then
                                t1=j
                                t2=i
! equation (52) in Yang & Wang (2012). 
                                mu2p=-pp+two*(t1*p1+t2*p2)
                            else
                                t1=i
                                t2=j  
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                        else        
                            If(f12342.lt.zero)then        
                                t1=i
                                t2=j
! equation (52) in Yang & Wang (2012). 
                                mu2p=pp+two*(t1*p1+t2*p2)
                            endif
                            If(f12342.gt.zero)then        
                                t1=j
                                t2=i      
! equation (52) in Yang & Wang (2012).    
                                mu2p=-pp+two*(t1*p1+t2*p2)                              
                            endif
                        endif 
                        !write(*,*)p,mu2p,abs(p-mu2p),pp,p1,p2! 
                        If(abs(p-mu2p).lt.1.D-4)goto 400
                    enddo
                enddo
                !11111111*****************************************************************************************
            400 continue
            index_p4(1)=-1
            index_p4(2)=-2
            index_p4(3)=0
            index_p4(4)=-4
           !*****pp part***************************************
            If(lambda.ne.zero)then        
                cases_int=2
                call weierstrass_int_J3(tinf,tmu,dd,del,-tplus,b4,index_p4,abs(pp),integ4,cases_int)
                call weierstrass_int_J3(tinf,tmu,dd,del,-tminus,b4,index_p4,abs(pp),integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                pp_phi=lambda*(pp/(one-mu_tp1*mu_tp1)+integ4(2)*Wmup-integ14(2)*Wmum)
            else 
                pp_phi=zero                 
            endif
            cases_int=4
            call weierstrass_int_J3(tinf,tmu,dd,del,h,b4,index_p4,abs(pp),integ,cases_int)
! equation (71) in Yang & Wang (2012). 
            pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
           !*****p1 part***************************************
            If(t1.eq.0)then        
                p1_phi=zero
                p1_t=zero
            else  
                If(lambda.ne.zero)then  
                    IF(PI1_phi .EQ. zero)THEN
                        cases_int=2        
                        call weierstrass_int_J3(tinf,infinity,dd,del,-tplus,b4,index_p4,PI0,integ4,cases_int)
                        call weierstrass_int_J3(tinf,infinity,dd,del,-tminus,b4,index_p4,PI0,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        PI1_phi=lambda*(PI0/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                    ENDIF 
! equation (51) in Yang & Wang (2012). 
                    p1_phi=PI1_phi-pp_phi 
                else 
                    p1_phi=zero             
                endif 
                IF(PI1_time .EQ. zero)THEN  
                    cases_int=4  
                    call weierstrass_int_J3(tinf,infinity,dd,del,h,b4,index_p4,PI0,integ,cases_int)
! equation (62) in Yang & Wang (2012). 
                    PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen) 
                ENDIF
! equation (51) in Yang & Wang (2012). 
                p1_t=PI1_time-pp_t 
            endif 
          !*****p2 part***************************************
            If(t2.eq.0)then
                p2_phi=zero
                p2_t=zero
            else
                IF(lambda.ne.zero)then  
                    IF(PI2_phi .EQ. zero)THEN  
                        cases_int=2        
                        call weierstrass_int_J3(tp2,tinf,dd,del,-tplus,b4,index_p4,PI2_p,integ4,cases_int)
                        call weierstrass_int_J3(tp2,tinf,dd,del,-tminus,b4,index_p4,PI2_p,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        PI2_phi=lambda*(PI2_p/(one-mu_tp1*mu_tp1)+integ4(2)*Wmup-integ14(2)*Wmum) 
                    ENDIF
! equation (51) in Yang & Wang (2012). 
                    p2_phi=PI2_phi+pp_phi
                ELSE
                    p2_phi=zero                
                ENDIF 

                IF(PI2_time .EQ. zero)THEN  
                    cases_int=4  
                    call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,index_p4,PI2_p,integ,cases_int)
! equation (71) in Yang & Wang (2012).  
                    PI2_time=a_spin**two*(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen) 
                ENDIF   
! equation (51) in Yang & Wang (2012).      
                p2_t=PI2_time+pp_t   
                !write(*,*)'ynogk=',tp2,tinf,h,dd
            endif   
        !**************************************************************
! equation (52) in Yang & Wang (2012). 
            !write(*,*)'ynogk=',pp_t,p1_t,p2_t,t1,t2
            If(mobseqmtp)then
                If(muobs.eq.mu_tp1)then  
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)                
                else
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)                
                endif 
            else
                If(f12342.lt.zero)then
                    phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                    timet=pp_t+two*(t1*p1_t+t2*p2_t)
                endif
                If(f12342.gt.zero)then                
                    phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                    timet=-pp_t+two*(t1*p1_t+t2*p2_t)        
                endif
            endif
            count_num=count_num+1
        ELSE 
            If(f12343_1.eq.f12343.and.f12342_1.eq.f12342.and.lambda_1.eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
            .and.muobs_1.eq.muobs.and.a_spin_1.eq.a_spin.and.scal_1.eq.scal)then  
        !***************************************************************************
                    t1=0
                    t2=0        
                    If(f12343.eq.zero.and.f12342.eq.zero.and.abs(muobs).eq.one)then
                        mucos = sign(one,muobs)
                        timet=zero      !this is because that mu==1 for ever
                        phyt=zero       !this is because that mu==1 for ever,this because that Theta_mu=-a^2(1-mu^2)
                        return          !so,mu must =+1 or -1 for ever.
                    endif
                    If(muobs.eq.zero.and.(abs(lambda).lt.abs(a_spin)).and.q.eq.zero)then
                        timet=zero
                        phyt=zero
                        mucos=zero
                        return                        
                    endif          
                    If(mu_tp1.eq.zero)then
                        !photons are confined in the equatorial plane, so the integrations about \theta are valished.
                        timet=zero
                        phyt=zero
                        mucos=zero
                        return
                    endif
                    If(a_spin.eq.zero)then
                        timet=zero
                        CALL phyt_schwatz(p,f12343,f12342,lambda,q,sinobs,muobs,scal,phyt,mucos,t1,t2)
                        return
                    endif
         
                    tmu=weierstrassP(p+PI01,g2,g3,dd,del)
                    mucos=mu_tp1+b0/(four*tmu-b1)
                    !get numbers of turn points of t1 and t2.
                    !111111111************************************************************  
                        index_p4(1)=0
                        cases_int=1    
                        call weierstrass_int_J3(tinf,tmu,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)                
                        pp=integ4(1)
                        p1=PI0-pp
                        p2=p_mt1_mt2-p1
                        !p1=zero
                        !p2=zero
                        Do j=0,10 
                            Do i=j,j+1
                                If(mobseqmtp)then
                                    If(muobs.eq.mu_tp1)then
                                        t1=j
                                        t2=i
! equation (54) in Yang & Wang (2012). 
                                        mu2p=-pp+two*(t1*p1+t2*p2) 
                                    else
                                        t1=i
                                        t2=j 
                                        mu2p=pp+two*(t1*p1+t2*p2) 
                                    endif
                                else        
                                    If(f12342.lt.zero)then        
                                        t1=i
                                        t2=j
! equation (54) in Yang & Wang (2012). 
                                        mu2p=pp+two*(t1*p1+t2*p2)
                                    endif
                                    If(f12342.gt.zero)then        
                                        t1=j
                                        t2=i     
! equation (54) in Yang & Wang (2012).      
                                        mu2p=-pp+two*(t1*p1+t2*p2)                             
                                    endif
                                endif  
                                !write(*,*)p,mu2p,t1,t2
                                If(abs(p-mu2p).lt.1.D-4)goto 410
                            enddo
                        enddo
                        !11111111********************************************* 
                    410 continue
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4 
                    !*****pp parts************************************
                    If(lambda.ne.zero)then
                        cases_int=2        
                        call weierstrass_int_J3(tinf,tmu,dd,del,-tplus,b4,index_p4,abs(pp),integ4,cases_int)
                        call weierstrass_int_J3(tinf,tmu,dd,del,-tminus,b4,index_p4,abs(pp),integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                        pp_phi=lambda*(pp/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum) 
                    else 
                        pp_phi=zero         
                    endif
                    cases_int=4
                    call weierstrass_int_J3(tinf,tmu,dd,del,h,b4,index_p4,abs(pp),integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                    pp_t=a_spin**two*(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                    !*****p1 parts************************************
                    If(t1.eq.0)then        
                        p1_phi=zero
                        p1_t=zero
                    else  
                        If(lambda.ne.zero)then  
                            IF(PI1_phi .EQ. zero)THEN
                                cases_int=2        
                                call weierstrass_int_J3(tinf,infinity,dd,del,-tplus,b4,index_p4,PI0,integ4,cases_int)
                                call weierstrass_int_J3(tinf,infinity,dd,del,-tminus,b4,index_p4,PI0,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                                PI1_phi=lambda*(PI0/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                            ENDIF
! equation (51) in Yang & Wang (2012). 
                            p1_phi=PI1_phi-pp_phi 
                        else 
                            p1_phi=zero            
                        endif  
                        IF(PI1_time .EQ. zero)THEN
                            cases_int=4  
                            call weierstrass_int_J3(tinf,infinity,dd,del,h,b4,index_p4,PI0,integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                            PI1_time=a_spin**two*(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                        ENDIF
! equation (51) in Yang & Wang (2012). 
                        p1_t=PI1_time-pp_t
                    endif 
                    !*****p2 parts************************************
                    If(t2.eq.0)then
                        p2_phi=zero
                        p2_t=zero
                    else        
                        If(lambda.ne.zero)then
                            IF(PI2_phi .EQ. zero)THEN 
                               cases_int=2
                               call weierstrass_int_J3(tp2,tinf,dd,del,-tplus,b4,index_p4,PI2_p,integ4,cases_int)
                               call weierstrass_int_J3(tp2,tinf,dd,del,-tminus,b4,index_p4,PI2_p,integ14,cases_int)
! equation (72) in Yang & Wang (2012). 
                               PI2_phi=lambda*(PI2_p/(one-mu_tp1**two)+integ4(2)*Wmup-integ14(2)*Wmum)
                            ENDIF
! equation (51) in Yang & Wang (2012). 
                            p2_phi=PI2_phi+pp_phi 
                        else 
                            p2_phi=zero                     
                        endif
                        IF(PI2_time .EQ. zero)THEN 
                            cases_int=4
                            call weierstrass_int_J3(tp2,tinf,dd,del,h,b4,index_p4,PI2_p,integ,cases_int)
! equation (71) in Yang & Wang (2012). 
                            PI2_time=a_spin**two*(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/sixteen)
                        ENDIF
! equation (51) in Yang & Wang (2012). 
                        p2_t=PI2_time+pp_t 
                    endif 
                !**************************************************************
! equation (52) in Yang & Wang (2012). 
                    If(mobseqmtp)then
                        If(muobs.eq.mu_tp1)then  
                            phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timet=-pp_t+two*(t1*p1_t+t2*p2_t)                
                        else
                            phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                            timet=pp_t+two*(t1*p1_t+t2*p2_t)                
                        endif 
                    else
                        If(f12342.lt.zero)then
                            phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)        
                            timet=pp_t+two*(t1*p1_t+t2*p2_t)
                        endif
                        If(f12342.gt.zero)then                
                            phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timet=-pp_t+two*(t1*p1_t+t2*p2_t)         
                        endif
                    endif                
            else
                count_num=1
                goto 30          
            endif                                
        ENDIF
        !write(*,*)'ff1=',phyt,timet,pp_phi,p1_phi,p2_phi,t1,t2
        RETURN
      END SUBROUTINE INTTPART

!********************************************************************************************
      SUBROUTINE phyt_schwatz(p,f3,f2,lambda,q,sinobs,muobs,scal,phyc_schwatz,mucos,t1,t2)
!******************************************************************************************** 
!*     PURPOSE:  Computes \mu part of integrals in coordinates \phi, expressed by equation (72) 
!*               in Yang & Wang (2012) with zero spin of black hole.    
!*     INPUTS:   p--------------independent variable, which must be nonnegative.
!*               f2-------------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               f3-------------p_\phi, which is the \phi component of four momentum of a photon 
!*                              measured under the LNRF, see equation (85) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyc_schwatz-----------value of integral \phi_\theta expressed by equation (71) in
!*                              Yang & Wang (2012).   
!*               mucos----------value of function \mu(p) with zero spin.
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.            
!*     ROUTINES CALLED: schwatz_int
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        USE constants
        implicit none
        Double precision phyc_schwatz,f3,f2,p,sinobs,muobs,pp,p1,p2,mu,&
                        lambda,AA,BB,scal,q,mu_tp1,mu_tp2,&
                        mu2p,mucos,Pt,PI1,PI1_phi,PI2_phi,f3_1,f2_1,lambda_1,q_1,&
                        sinobs_1,muobs_1,scal_1,pp_phi,p1_phi,p2_phi
        !parameter(zero=0.D0,one=1.D0,two=2.D0)
        integer  :: t1,t2,i,j,count_num=1        
        logical :: mobseqmtp
        save :: PI1,PI1_phi,PI2_phi,Pt,f3_1,f2_1,lambda_1,q_1,pp_phi,p1_phi,p2_phi,&
                sinobs_1,muobs_1,scal_1,mobseqmtp,AA,BB,&
                mu_tp1,mu_tp2 

60 continue 
      IF(count_num .EQ. 1)THEN
          f3_1=f3
          f2_1=f2
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs 
          scal_1=scal
          t1=0
          t2=0         
          mobseqmtp=.false.
          If(q.gt.zero)then        
              AA=sqrt((lambda**two+q)/q)
              BB=sqrt(q)
        !*****************************************************
              If(f2.lt.zero)then
                  mu=sin(asin(muobs*AA)+p*BB*AA)/AA        
              else
                  If(f2.eq.zero)then
                      mu=cos(p*AA*BB)*muobs
                  else                              
                      mu=sin(asin(muobs*AA)-p*AA*BB)/AA        
                  endif        
              endif
              mucos = mu  
        !****************************************************
              If(f2.ne.zero)then
                  mu_tp1=sqrt(q/(lambda**two+q))
                  mu_tp2=-mu_tp1        
              else
                  mu_tp1=abs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
              endif
              If(abs(muobs).eq.one)mobseqmtp=.true.        

              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  phyc_schwatz=zero
                  return
              endif

              !***************************************************
              PI1=(PI/two-asin(muobs/mu_tp1))*mu_tp1/BB        
              Pt=PI*mu_tp1/BB        
              pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB        
              p1=PI1-pp
              p2=Pt-p1 
              PI1_phi=zero
              PI2_phi=zero
              Do j=0,100 
                  Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              mu2p=-pp+two*(t1*p1+t2*p2)
                          else
                              t1=i
                              t2=j
                              mu2p=pp+two*(t1*p1+t2*p2)
                          endif
                      else        
                          If(f2.lt.zero)then        
                              t1=i
                              t2=j
                              mu2p=pp+two*(t1*p1+t2*p2)
                          endif
                          If(f2.gt.zero)then        
                              t1=j
                              t2=i    
                              mu2p=-pp+two*(t1*p1+t2*p2)                            
                          endif
                      endif  
                      If(abs(p-mu2p).lt.1.D-4)goto 300
                  enddo
              enddo
              !*************************************************************** 
              300 continue
              If(lambda.eq.zero)then 
                  phyc_schwatz = zero
                  return
              endif
              pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
              If(t1.eq.0)then
                  p1_phi=zero
              else
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                  ENDIF
                  p1_phi=PI1_phi-pp_phi       
              endif
              If(t2.eq.0)then
                  p2_phi=zero
              else
                  IF(PI2_phi .EQ. zero)THEN
                      PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                  ENDIF
                  p2_phi=PI2_phi+pp_phi 
              endif
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then  
                      phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                  else
                      phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                  endif         
              else
                 If(f2.lt.zero)then
                      phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif         
                 If(f2.gt.zero)then
                      phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif          
              endif
          else
              !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
              !                'value, the program should be',&  
              !                'stoped! and q = ',q
              !stop
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyc_schwatz = zero
          endif        
      ELSE
          IF(f3_1.eq.f3.and.f2_1.eq.f2.and.lambda_1.eq.lambda.and.q_1.eq.q.and.sinobs_1.eq.sinobs&
          .and.muobs_1.eq.muobs.and.scal_1.eq.scal)THEN
              If(q.gt.zero)then         
        !*****************************************************
                  If(f2.lt.zero)then
                      mu=sin(asin(muobs*AA)+p*BB*AA)/AA        
                  else
                      If(f2.eq.zero)then
                          mu=cos(p*AA*BB)*muobs
                      else                              
                          mu=sin(asin(muobs*AA)-p*AA*BB)/AA        
                      endif        
                  endif
                  mucos = mu  
        !****************************************************  
                  If(mu_tp1.eq.zero)then
                  !photons are confined in the equatorial plane, 
                  !so the integrations about !\theta are valished.
                      phyc_schwatz=zero
                      return
                  endif

                  !***************************************************  
                  pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB        
                  p1=PI1-pp
                  p2=Pt-p1  
                  Do j=0,100 
                      Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              mu2p=-pp+two*(t1*p1+t2*p2) 
                          else
                              t1=i
                              t2=j
                              mu2p=pp+two*(t1*p1+t2*p2) 
                          endif
                      else        
                          If(f2.lt.zero)then        
                              t1=i
                              t2=j
                              mu2p=pp+two*(t1*p1+t2*p2)
                          endif
                          If(f2.gt.zero)then        
                              t1=j
                              t2=i      
                              mu2p=-pp+two*(t1*p1+t2*p2)                            
                          endif
                      endif   
                      If(abs(p-mu2p).lt.1.D-4)goto 310
                      enddo
                  enddo
                  !************************************************************ 
310 continue
                  If(lambda.eq.zero)then 
                      phyc_schwatz = zero
                      return
                  endif
                  pp_phi=lambda*schwatz_int(muobs,mu,AA)/BB
                  If(t1.eq.0)then
                      p1_phi=zero
                  else
                      IF(PI1_phi .eq. zero)THEN
                          PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB
                      ENDIF
                      p1_phi=PI1_phi-pp_phi       
                  endif
                  If(t2.eq.0)then
                      p2_phi=zero
                  else
                      IF(PI2_phi .EQ. zero)THEN
                          PI2_phi=lambda*schwatz_int(mu_tp2,muobs,AA)/BB
                      ENDIF
                      p2_phi=PI2_phi+pp_phi 
                  endif
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then  
                          phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                      else
                          phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)                
                      endif         
                  else
                      If(f2.lt.zero)then
                          phyc_schwatz=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      endif         
                      If(f2.gt.zero)then
                          phyc_schwatz=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                      endif          
                  endif
              else
                  !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
                  !                'value, the program should be',&  
                  !                'stoped! and q = ',q
                  !stop
                  mucos=muobs 
                  t1 = 0
                  t2 = 0
                  phyc_schwatz = zero
              endif                  
          ELSE
              count_num=1
              goto 60
          ENDIF
      ENDIF                                
      return
      End SUBROUTINE phyt_schwatz 
!************************************************************************* 
      Function schwatz_int(y,x,AA) 
!************************************************************************* 
!*     PURPOSE:  Computes \int^x_y dt/(1-t^2)/sqrt(1-AA^2*t^2) and AA .gt. 1  
!*     INPUTS:   components of above integration.      
!*     OUTPUTS:  valve of integral.             
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      USE constants
      implicit none
      Double precision y,x,yt,xt,AA,schwatz_int,ppx,ppy,A2 

      xt=x
      yt=y        
      If(yt.eq.xt)then
          schwatz_int=0.D0
          return        
      endif
      If(abs(AA).ne.one)then
          A2=AA*AA
          ppx=atan(sqrt(A2-one)*xt/sqrt(abs(one-A2*xt*xt)))
          ppy=atan(sqrt(A2-one)*yt/sqrt(abs(one-A2*yt*yt)))
          schwatz_int=(ppx-ppy)/sqrt(A2-one) 
      ELse
          If(abs(xt).eq.one)then
              schwatz_int=infinity
          Else
              If(abs(yt).eq.one)then
                  schwatz_int=-infinity
              Else
                  ppx=xt/sqrt(abs(one-xt*xt))
                  ppy=yt/sqrt(abs(one-yt*yt))
                  schwatz_int=ppx-ppy 
              endif        
          Endif           
      Endif
      return
      End Function schwatz_int   

!********************************************************************************************
      SUBROUTINE INTRPART(p,f1234r,f1234t,lambda,q,sinobs,muobs,a_spin,&
                            robs,scal,phyr,timer,affr,r_coord,t1,t2)
!******************************************************************************************** 
!*     PURPOSE:  Computes r part of integrals in coordinates \phi, t and affine parameter \sigma,
!*               expressed by equations (62), (63), (65), (67), (68) and (69) in Yang & Wang (2012).    
!*     INPUTS:   p--------------independent variable, which must be nonnegative. 
!*               f1234r---------p_r, which is the r component of four momentum of a photon 
!*                              measured under the LNRF, see equation (83) in Yang & Wang (2012).
!*               f1234t---------p_\theta, which is the \theta component of four momentum of a photon 
!*                              measured under the LNRF, see equation (84) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  phyr-----------value of integral \phi_r expressed by equation (65) or (69) in
!*                              Yang & Wang (2012).  
!*               affr-----------value of integral \sigma_r expressed by equation (62) or (67) in
!*                              Yang & Wang (2012).  
!*               timer----------value of integral t_r expressed by equation (63) or (68) in
!*                              Yang & Wang (2012).  
!*               r_coord--------value of function r(p).
!*               t1,t2----------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.            
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
        USE constants
        IMPLICIT NONE
        DOUBLE PRECISION phyr,radius,p,sinobs,muobs,a_spin,rhorizon,q,lambda,integ4(4),&
                         cc,b0,b1,b2,b3,g2,g3,tobs,tp,pp,p1,p2,PI0,E_add,E_m,&
                         u,v,w,L1,L2,thorizon,m2,pinf,sn,cn,dn,r_add,r_m,B_add,B_m,D_add,D_m,&
                         y,x,f1,g1,h1,f2,h2,a5,b5,a4,b4,robs,&
                         scal,tinf,integ04(4),integ14(4),integ5(5),integ15(5),&
                         r_tp1,r_tp2,t_inf,tp2,f1234r,f1234t,p_temp,PI0_obs_inf,PI0_total,PI0_obs_hori,&
                         PI0_obs_tp2,PI01,timer,affr,r_coord,cr,dr,rff_p,&
                         Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,pp_time,pp_phi,pp_aff,p1_phi,p1_time,p1_aff,&
                         p2_phi,p2_time,p2_aff,time_temp,sqt3,p_tp1_tp2,PI2_p,PI1_p,&
                         PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff        
        DOUBLE PRECISION f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,scal_1
        !PARAMETER(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0)
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER ::  reals,i,j,t1,t2,index_p4(4),index_p5(5),del,cases_int,cases,count_num=1
        LOGICAL :: robs_eq_rtp,indrhorizon
        SAVE :: f1234r_1,f1234t_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,scal_1,rhorizon,r_add,r_m,a4,b4,B_add,&
                B_m,robs_eq_rtp,indrhorizon,r_tp1,r_tp2,reals,cases,bb,b0,b1,b2,b3,g2,g3,tobs,thorizon,&
                tp2,tinf,dd,E_add,E_m,D_add,D_m,PI0,PI0_obs_inf,PI0_total,PI0_obs_hori,PI0_obs_tp2,del,&
                u,v,w,L1,L2,m2,t_inf,pinf,f1,g1,h1,f2,h2,b5,Ap,Am,h,wp,wm,wbarp,wbarm,hm,hp,a5,cc,&
                PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_aff,PI2_aff,PI2_p,PI1_p,p_tp1_tp2,sqt3          

  40 continue        
        If(count_num.eq.1)then
                f1234r_1=f1234r        
                f1234t_1=f1234t
                lambda_1=lambda
                q_1=q
                muobs_1=muobs
                sinobs_1=sinobs
                a_spin_1=a_spin
                robs_1=robs
                scal_1=scal
            !************************************************************************************
                rhorizon=one+sqrt(one-a_spin**two)
! equation (64) in Yang & Wang (2012).
                r_add=rhorizon
                r_m=one-sqrt(one-a_spin**two)      
! equation (64) in Yang & Wang (2012).  
                B_add=(two*r_add-a_spin*lambda)/(r_add-r_m)        
                B_m=(two*r_m-a_spin*lambda)/(r_add-r_m)      
! equation (64) in Yang & Wang (2012).
                Ap=(r_add*(four-a_spin*lambda)-two*a_spin**two)/sqrt(one-a_spin**two)
                Am=(r_m*(four-a_spin*lambda)-two*a_spin**two)/sqrt(one-a_spin**two)
                b4=one
                a4=zero
                cc=a_spin**2-lambda**2-q
                robs_eq_rtp=.false.
                indrhorizon=.false.
                call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                                  reals,robs_eq_rtp,indrhorizon,cases,bb)
! equation (55) in Yang & Wang (2012).
                PI1_phi=zero
                PI2_phi=zero
                PI1_time=zero
                PI2_time=zero
                PI1_aff=zero
                PI2_aff=zero 
!** R(r)=0 has real roots and turning points exists in radial r.
                If(reals.ne.0)then  
! equations (35)-(38) in Yang & Wang (2012).
                    b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
                    b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                    b2=four/three*r_tp1
                    b3=one
                    g2=three/four*(b1**2-b0*b2)
                    g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
! equation (39) in Yang & Wang (2012).
                    If(robs-r_tp1.ne.zero)then        
                        tobs=b0/four/(robs-r_tp1)+b1/four
                    else
                        tobs=infinity
                    endif 
                    If(rhorizon-r_tp1.ne.zero)then
                        thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                    else
                        thorizon=infinity         
                    endif
                    tp2=b0/four/(r_tp2-r_tp1)+b1/four   
                    tinf=b1/four
                    h=-b1/four        
! equation (64), (66) and (70) in Yang & Wang (2012).
                    call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)        
                        E_add=b0/(four*(r_add-r_tp1))+b1/four
                          E_m=b0/(four*(r_m-r_tp1))+b1/four
                        D_add=b0/(four*(r_tp1-r_add)**2)
                          D_m=b0/(four*(r_tp1-r_m)**2)

                               wp=one/(r_tp1-r_add)
                               wm=one/(r_tp1-r_m)
                            wbarp=b0/four/(r_tp1-r_add)**two
                            wbarm=b0/four/(r_tp1-r_m)**two
                               hp=b0/four/(r_add-r_tp1)+b1/four
                               hm=b0/four/(r_m-r_tp1)+b1/four

                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tobs,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int) 
! equation (42) in Yang & Wang (2012).
                    PI0=integ04(1)   
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tinf,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            PI0_obs_inf=integ04(1)
                            If(p.lt.PI0_obs_inf)then    
! equation (41) in Yang & Wang (2012).    
                                tp=weierstrassP(p+PI0,g2,g3,dd,del) 
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=-p                                  
                            else
                                tp=tinf! !Goto infinity, far away. 
                                r_coord = infinity
                                pp=-PI0_obs_inf 
                            endif
                            t1=0
                            t2=0          
                        ELSE 
                            If(.not.indrhorizon)then
                                index_p4(1)=0
                                cases_int=1 
                                call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                                PI0_total=PI0+integ14(1)
                                t2=0
                                If(p.le.PI0)then
                                    t1=0
                                    pp=p  
! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=1
                                    PI1_p=PI0        
                                    If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).     
                                        tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=two*PI0-p
                                        p1=abs(p-PI0)
                                    else        
                                        tp=tinf !Goto infinity, far away.
                                        r_coord = infinity
                                        pp=-PI0_total+two*PI0 
                                        p1=pI0_total-PI0
                                    endif        
                                endif        
                            ELSE     !f1234r<0, photon will fall into black hole unless something encountered.        
                                index_p4(1)=0                
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                                PI0_obs_hori=integ04(1)
                                If(p.lt.PI0_obs_hori)then     
! equation (41) in Yang & Wang (2012).   
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0  
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).
                            tp=weierstrassP(p+PI01,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)        
                            index_p4(1)=0
                            cases_int=1        
                            call weierstrass_int_J3(tobs,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int)
                            call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            pp=integ4(1)
! equation (57) in Yang & Wang (2012).
                            p_tp1_tp2=integ14(1) 
                            PI2_p=p_tp1_tp2-PI0 
                            PI1_p=PI0 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1 
                            !p1=zero
                            !p2=zero
                        !*************************************************************************************
! equation (58) in Yang & Wang (2012).
                            Do j=0,100
                                Do i=j,j+1
                                    If(robs_eq_rtp)then        
                                        If(robs.eq.r_tp1)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        else
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    else
                                        If(f1234r.gt.zero)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        endif
                                        If(f1234r.lt.zero)then
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    endif  
                                    If(abs(p-p_temp).lt.1.D-4)goto 200
                                Enddo
                            Enddo
                        !*************************************************************************************
                        200     continue                    
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then
                                index_p4(1)=0
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                                PI0_obs_hori=integ04(1)
                                If(p.lt.PI0_obs_hori)then   
! equation (41) in Yang & Wang (2012).     
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                            ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.     
                                index_p4(1)=0
                                cases_int=1        
                                call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)  
                                call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                                PI0_obs_tp2=integ04(1)        
                                PI2_p=PI0_obs_tp2
                                PI0_total=integ14(1)+PI0_obs_tp2
                                If(p.le.PI0_obs_tp2)then
                                    t1=0
                                    t2=0
                                    pp=-p 
! equation (41) in Yang & Wang (2012).
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=0
                                    t2=1
                                    If(p.lt.PI0_total)then  
! equation (41) in Yang & Wang (2012).      
                                        tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=p-two*PI0_obs_tp2
                                        p2=p-PI0_obs_tp2
                                    else        
                                        tp=thorizon !Fall into black hole. 
                                        r_coord = rhorizon
                                        pp=PI0_total-two*PI0_obs_tp2
                                        p2=PI0_total-PI0_obs_tp2
                                    endif        
                                endif        
                            ENDIF
                        ENDIF                             
                    END SELECT  
              !****************************************************************** 
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4
                !pp part ***************************************************        
                    cases_int=4
                    call weierstrass_int_J3(tobs,tp,dd,del,h,b4,index_p4,abs(pp),integ4,cases_int)

! equation (62) in Yang & Wang (2012).
                    pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two  
! equation (63) in Yang & Wang (2012).
                    pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                    time_temp=pp*(-Am*wm)           
        
                    cases_int=2        
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_add,b4,index_p4,abs(pp),integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                    pp_time=pp_time-Ap*wbarp*integ4(2)
                    IF(a_spin.NE.zero)THEN
                        call weierstrass_int_J3(tobs,tp,dd,del,-E_m,b4,index_p4,abs(pp),integ14,cases_int)
! equation (63) in Yang & Wang (2012).
                        pp_time=pp_time+Am*wbarm*integ14(2)+time_temp
! equation (65) in Yang & Wang (2012).               
                        pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                    -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                    ELSE
                        pp_phi=zero
                    ENDIF                
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).
                        pp_phi=pp_phi+pp*lambda
                    ENDIF        
                    !p1 part *******************************************************
                    IF(t1 .EQ. 0)THEN
                        p1_phi=ZERO
                        p1_time=ZERO
                        p1_aff=ZERO
                    ELSE
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)
! equation (62) in Yang & Wang (2012).
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI0*r_tp1**two    
! equation (63) in Yang & Wang (2012).     
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff        
                            time_temp=PI0*(-Am*wm)
        
                            cases_int=2        
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                            PI1_time=PI1_time-Ap*wbarp*integ4(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
! equation (63) in Yang & Wang (2012).       
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                              -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).       
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        ENDIF 
! equation (55) in Yang & Wang (2012).
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi 
                    ENDIF
                !p2 part *******************************************************
                    IF(t2.EQ.ZERO)THEN
                        p2_phi=ZERO
                        p2_time=ZERO
                        p2_aff=ZERO
                    ELSE
                        IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p4,PI2_p,integ4,cases_int)
! equation (62) in Yang & Wang (2012).       
                            PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI2_p*r_tp1**two  
! equation (63) in Yang & Wang (2012).       
                            PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff        
                            time_temp=PI2_p*(-Am*wm)

                            cases_int=2        
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,b4,index_p4,PI2_p,integ4,cases_int)
! equation (63) in Yang & Wang (2012).
                            PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,b4,index_p4,PI2_p,integ14,cases_int)
! equation (63) in Yang & Wang (2012).       
                                PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).
                                PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                           -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI2_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).       
                                PI2_phi=PI2_phi+PI2_p*lambda
                            ENDIF
                        ENDIF
! equation (55) in Yang & Wang (2012).       
                        p2_aff=PI2_aff+pp_aff
                        p2_time=PI2_time+pp_time
                        p2_phi=PI2_phi+pp_phi
                    ENDIF
                    !phi, aff,time part *******************************************************
! equation (56) in Yang & Wang (2012).       
                    If(f1234r.ne.zero)then
                        phyr=sign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=sign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=sign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        If(robs.eq.r_tp1)then
                            phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)        
                        else
                            phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        endif        
                    endif 
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon                  
                                endif
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity          
                                endif                        
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                                sqt3=sqrt(three)        
                            If(f1234r.lt.zero)then
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then        
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                            else        
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                                PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -Log(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then        
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif                        
                            endif                
                        endif
                    endif                                   
                ELSE            
! equation (44) in Yang & Wang (2012).   !equation R(r)=0 has no real roots. we use the Legendre elliptic 
                    u=real(bb(4))        !integrations and functions to compute the calculations.
                    w=abs(aimag(bb(4)))
                    v=abs(aimag(bb(2)))
                    If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).       
                        L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                        L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).       
                        thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012).   
                        m2=(L1-L2)/L1
                        tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                        t_inf=sqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).       
                        pinf=EllipticF(tinf,m2)/(w*sqrt(L1))
                        call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn)
                        f1=u**2+w**2
                        g1=-two*u
                        h1=one
                        f2=u**2+v**2
                        g2=-g1
                        h2=one
                        a5=zero
                        b5=one
                        index_p5(1)=-1
                        index_p5(2)=-2
                        index_p5(3)=2
                        index_p5(4)=-4
                        index_p5(5)=4
                        IF(f1234r.lt.zero)THEN
                            PI0=pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))                                
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).       
                                y=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = y 
                                pp=p  
                            else
                                y=rhorizon
                                r_coord = y
                                pp=PI0
                            endif         
                            x=robs 
                        ELSE
                            PI0=EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).       
                                x=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = x
                                pp=p
                            else
                                x=infinity 
                                r_coord = x
                                pp=PI0
                            endif
                            y=robs         
                        ENDIF 
                        !affine parameter part integration **********************************************
                        cases_int=5    
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,index_p5,abs(pp),integ5,cases_int)
! equation (67) in Yang & Wang (2012).   
                        affr=integ5(5)
! equation (68) in Yang & Wang (2012).  
                        timer=two*integ5(3)+four*pp+affr
                        cases_int=2    
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_add,b5,index_p5,abs(pp),integ5,cases_int)
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_m,b5,index_p5,abs(pp),integ15,cases_int)
                        !phy part**************************************************************************
! equation (68) in Yang & Wang (2012).   
                        timer=timer+Ap*integ5(2)-Am*integ15(2)
! equation (69) in Yang & Wang (2012).   
                        phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))        
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).   
                            phyr=phyr+pp*lambda
                        ENDIF 
                    ELSE
                        If(f1234r.lt.zero)then
                            PI0=(atan(robs/w)-atan(rhorizon/w))/w        
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)-p*w)
                                r_coord = radius 
                            else
                                radius=rhorizon
                                r_coord = radius
                            endif
                            !timer part ****************************************    
                            y=radius
                            x=robs
                        ELSE
                            PI0=(PI/two-atan(robs/w))/w        
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)+p*w)
                                r_coord = radius
                            else
                                radius=infinity
                                r_coord = radius
                            endif
                            !timer part ************************************************************************8
                            y=robs
                            x=radius
                        ENDIF
                        pp_time=(x-y)+atan(x/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                                  atan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))
                        pp_time=pp_time+log(x**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two))-&
                                  (log(y**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two)))        
                        timer=pp_time+Ap*log(abs(x-r_add))/(w**two+r_add**two)-Am*log(abs(x-r_m))/(w**two+r_m**two)-&
                                   (Ap*log(abs(y-r_add))/(w**two+r_add**two)-Am*log(abs(y-r_m))/(w**two+r_m**two))  
                        !affine parameter part *****************************************************************
                        affr=(x-y)-w*atan(x/w)+w*atan(y/w)          
                        !phy part ******************************************************************
                        IF(a_spin .NE. zero)THEN
                            phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                                               *(atan(x/w)-atan(y/w))+&
                                 log(abs(x-r_add)/sqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                                 log(abs(y-r_add)/sqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                                 log(abs(x-r_m)/sqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                                 log(abs(y-r_m)/sqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                            phyr=phyr*a_spin 
                        ELSE
                            phyr=zero
                        ENDIF        
                        If(muobs.eq.zero.and.f1234t.eq.zero)then         
                            phyr=phyr+lambda*(atan(x/w)-atan(y/w))/w
                        ENDIF   
                    ENDIF                        
                ENDIF
                count_num=count_num+1
        !*****************************************************************************************
        else 
        !*****************************************************************************************
            If(f1234r.eq.f1234r_1.and.f1234t.eq.f1234t_1.and.lambda.eq.lambda_1.and.q.eq.q_1.and.&
            sinobs.eq.sinobs_1.and.muobs.eq.muobs_1.and.a_spin.eq.a_spin_1.and.robs.eq.robs_1.and.scal.eq.scal_1)then
                !***********************************************************************
                If(reals.ne.0)then  !** R(r)=0 has real roots and turning points exists in radial r. 
!used in the geodesics I'm trying to understand March 6 2017
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity. 
                            If(p.lt.PI0_obs_inf)then  
! equation (41) in Yang & Wang (2012).             
                                tp=weierstrassP(p+PI0,g2,g3,dd,del)        
                                r_coord = r_tp1+b0/(four*tp-b1)
                                pp=-p                                  
                            else
                                tp=tinf! !Goto infinity, far away. 
                                r_coord = infinity
                                pp=-PI0_obs_inf 
                            endif
                            t1=0
                            t2=0          
                        ELSE 
                            If(.not.indrhorizon)then 
                                t2=0
                                If(p.le.PI0)then
                                    t1=0
                                    pp=p  
! equation (41) in Yang & Wang (2012).       
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=1        
                                    If(p.lt.PI0_total)then   
! equation (41) in Yang & Wang (2012).            
                                        tp=weierstrassP(p-PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=two*PI0-p
                                        p1=abs(p-PI0)
                                    else        
                                        tp=tinf !Goto infinity, far away.
                                        r_coord = infinity
                                        pp=-PI0_total+two*PI0 
                                        p1=pI0_total-PI0
                                    endif        
                                endif        
                            ELSE     !f1234r<0, photon will fall into black hole unless something encountered.         
                                If(p.lt.PI0_obs_hori)then  
! equation (41) in Yang & Wang (2012).             
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0  
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).       
                            tp=weierstrassP(p+PI01,g2,g3,dd,del)
                            r_coord = r_tp1+b0/(four*tp-b1)        
                            index_p4(1)=0
                            cases_int=1        
                            call weierstrass_int_J3(tobs,tp,dd,del,a4,b4,index_p4,rff_p,integ4,cases_int) 
! equation (57) in Yang & Wang (2012).       
                            pp=integ4(1) 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1  
                        !*************************************************************************************
! equation (58) in Yang & Wang (2012).       
                            Do j=0,100
                                Do i=j,j+1
                                    If(robs_eq_rtp)then        
                                        If(robs.eq.r_tp1)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        else
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    else
                                        If(f1234r.gt.zero)then
                                            t1=j
                                            t2=i
                                            p_temp=-pp+two*(t1*p1+t2*p2)
                                        endif
                                        If(f1234r.lt.zero)then
                                            t1=i
                                            t2=j
                                            p_temp=pp+two*(t1*p1+t2*p2)
                                        endif
                                    endif  
                                    If(abs(p-p_temp).lt.1.D-4)goto 210
                                Enddo
                            Enddo
                        !*************************************************************************************
                        210        continue                    
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then 
                                If(p.lt.PI0_obs_hori)then    
! equation (41) in Yang & Wang (2012).        
                                    tp=weierstrassP(p-PI0,g2,g3,dd,del)        
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                    pp=p                                  
                                else
                                    tp=thorizon! !Fall into black hole.
                                    r_coord = rhorizon
                                    pp=PI0_obs_hori
                                endif
                                t1=0
                                t2=0
                            ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.  
                                If(p.le.PI0_obs_tp2)then
                                    t1=0
                                    t2=0
                                    pp=-p 
! equation (41) in Yang & Wang (2012).    
                                    tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                    r_coord = r_tp1+b0/(four*tp-b1)
                                else
                                    t1=0
                                    t2=1
                                    If(p.lt.PI0_total)then    
! equation (41) in Yang & Wang (2012).        
                                        tp=weierstrassP(p+PI0,g2,g3,dd,del)
                                        r_coord = r_tp1+b0/(four*tp-b1)
                                        pp=p-two*PI0_obs_tp2
                                        p2=p-PI0_obs_tp2
                                    else        
                                        tp=thorizon !Fall into black hole. 
                                        r_coord = rhorizon
                                        pp=PI0_total-two*PI0_obs_tp2
                                        p2=PI0_total-PI0_obs_tp2
                                    endif        
                                endif        
                            ENDIF
                        ENDIF                             
                    END SELECT 
              !****************************************************************** 
                    index_p4(1)=-1
                    index_p4(2)=-2
                    index_p4(3)=0
                    index_p4(4)=-4
                !pp part *************************************************** 
                   
                    cases_int=4
                    call weierstrass_int_J3(tobs,tp,dd,del,h,b4,index_p4,abs(pp),integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
                    pp_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pp*r_tp1**two  
! equation (63) in Yang & Wang (2012).           
                    pp_time=integ4(2)*b0/two+pp*(two*r_tp1+four+Ap*wp)+pp_aff
                    time_temp=pp*(-Am*wm)         
        
                    cases_int=2        
                    call weierstrass_int_J3(tobs,tp,dd,del,-E_add,b4,index_p4,abs(pp),integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                    pp_time=pp_time-Ap*wbarp*integ4(2)         
                    IF(a_spin.NE.zero)THEN
                        call weierstrass_int_J3(tobs,tp,dd,del,-E_m,b4,index_p4,abs(pp),integ14,cases_int)
! equation (63) in Yang & Wang (2012).    
                        pp_time=pp_time+Am*wbarm*integ14(2)+time_temp     
! equation (65) in Yang & Wang (2012).                   
                        pp_phi=pp*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                    -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                    ELSE
                        pp_phi=zero
                    ENDIF                
                    IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                        pp_phi=pp_phi+pp*lambda
                    ENDIF        
                 !p1 part *******************************************************
                    IF(t1 .EQ. 0)THEN
                        p1_phi=ZERO
                        p1_time=ZERO
                        p1_aff=ZERO

!added by Pieter ************************************************************************
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN !original
!here, PI1_time is calculated the first point after the turning point (t1>0)
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)
! equation (62) in Yang & Wang (2012). 
!********************************************* old ******************************************   
                            !PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+p1*r_tp1**two
!********************************************* old ******************************************   

!********************************************* new ******************************************   
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pI0*r_tp1**two  
!********************************************* new ******************************************    
! equation (63) in Yang & Wang (2012).          
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff        
                            time_temp=PI0*(-Am*wm)
        
                            cases_int=2        
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                            PI1_time=PI1_time-Ap*wbarp*integ4(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
! equation (65) in Yang & Wang (2012).    
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp   
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                              -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                        ENDIF
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi 

!end added by Pieter ************************************************************************


                    ELSE
                        IF(PI1_aff .EQ. zero .AND. PI1_time .EQ. zero)THEN !original
!here, PI1_time is calculated the first point after the turning point (t1>0)
                            cases_int=4
                            call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p4,PI0,integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
!********************************************* old ******************************************   
                            !PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+p1*r_tp1**two
!********************************************* old ******************************************   

!********************************************* new ******************************************   
                            PI1_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+pI0*r_tp1**two  
!********************************************* new ******************************************      
! equation (63) in Yang & Wang (2012).          
                            PI1_time=integ4(2)*b0/two+PI0*(two*r_tp1+four+Ap*wp)+PI1_aff        
                            time_temp=PI0*(-Am*wm)
        
                            cases_int=2        
                            call weierstrass_int_J3(tobs,infinity,dd,del,-E_add,b4,index_p4,PI0,integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                            PI1_time=PI1_time-Ap*wbarp*integ4(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tobs,infinity,dd,del,-E_m,b4,index_p4,PI0,integ14,cases_int)
! equation (65) in Yang & Wang (2012).    
                                PI1_time=PI1_time+Am*wbarm*integ14(2)+time_temp   
                                PI1_phi=PI0*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                              -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI1_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                                PI1_phi=PI1_phi+PI0*lambda
                            ENDIF
                            !write(*,*)'hier:',PI1_time
                        ENDIF 
! equation (55) in Yang & Wang (2012).    
                        p1_aff=PI1_aff-pp_aff
                        p1_time=PI1_time-pp_time
                        P1_phi=PI1_phi-pp_phi 
                    ENDIF
                !p2 part *******************************************************
                    IF(t2.EQ.ZERO)THEN
                        p2_phi=ZERO
                        p2_time=ZERO
                        p2_aff=ZERO
                    ELSE
                        IF(PI2_aff .EQ. zero .AND. PI2_time .EQ. zero)THEN
                            cases_int=4
                            call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p4,PI2_p,integ4,cases_int)
! equation (62) in Yang & Wang (2012).    
                            PI2_aff=integ4(4)*b0**two/sixteen+integ4(2)*b0*r_tp1/two+PI2_p*r_tp1**two   
! equation (63) in Yang & Wang (2012).          
                            PI2_time=integ4(2)*b0/two+PI2_p*(two*r_tp1+four+Ap*wp)+PI2_aff        
                            time_temp=PI2_p*(-Am*wm)

                            cases_int=2        
                            call weierstrass_int_J3(tp2,tobs,dd,del,-E_add,b4,index_p4,PI2_p,integ4,cases_int)
! equation (63) in Yang & Wang (2012).    
                            PI2_time=PI2_time-Ap*wbarp*integ4(2)!+Am*wbarm*integ14(2) 
                            IF(a_spin.NE.zero)THEN
                                call weierstrass_int_J3(tp2,tobs,dd,del,-E_m,b4,index_p4,PI2_p,integ14,cases_int)
! equation (63) in Yang & Wang (2012).    
                                PI2_time=PI2_time+Am*wbarm*integ14(2)+time_temp   
! equation (65) in Yang & Wang (2012).    
                                PI2_phi=PI2_p*a_spin*(B_add/(r_tp1-r_add)-B_m/(r_tp1-r_m))&
                                           -a_spin*B_add*D_add*integ4(2)+a_spin*B_m*D_m*integ14(2)
                            ELSE
                                PI2_phi=zero
                            ENDIF
                            IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                                PI2_phi=PI2_phi+PI2_p*lambda
                            ENDIF
                        ENDIF
! equation (55) in Yang & Wang (2012).    
                        p2_aff=PI2_aff+pp_aff
                        p2_time=PI2_time+pp_time
                        p2_phi=PI2_phi+pp_phi
                    ENDIF                    
                    !phi, aff,time part *******************************************************
! equation (56) in Yang & Wang (2012).    
                    If(f1234r.ne.zero)then
                        phyr=sign(one,-f1234r)*pp_phi+two*(t1*p1_phi+t2*p2_phi)
                        timer=sign(one,-f1234r)*pp_time+two*(t1*p1_time+t2*p2_time)
                        affr=sign(one,-f1234r)*pp_aff+two*(t1*p1_aff+t2*p2_aff)
                    else
                        If(robs.eq.r_tp1)then
                            phyr=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=-pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=-pp_aff+two*(t1*p1_aff+t2*p2_aff)        
                        else
                            phyr=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                            timer=pp_time+two*(t1*p1_time+t2*p2_time)
                            affr=pp_aff+two*(t1*p1_aff+t2*p2_aff)
                        endif        
                    endif 
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                If(p.lt.one/rhorizon-one/robs)then
                                    radius=robs/(robs*p+one)
                                else
                                    radius=rhorizon                  
                                endif
                            else
                                If(p.lt.one/robs)then
                                    radius=robs/(one-robs*p)
                                else
                                    radius=infinity          
                                endif                        
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                                sqt3=sqrt(three)        
                            If(f1234r.lt.zero)then
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(three*sqt3*p)+two/sqt3
                                If(p.ne.zero)then        
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=robs!infinity
                                endif
                            else        
                                cr=-three*abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)-sqt3
                                dr=-abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(three-robs))*exp(-three*sqt3*p)+two/sqt3
                                PI0=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -Log(one+two/sqt3)/three/sqt3
                                If(p.lt.PI0)then        
                                    radius=(three+cr*dr+sqrt(9.D0+6.D0*cr*dr+cr**two))/(dr**two-one)
                                else
                                    radius=infinity
                                endif                        
                            endif                
                        endif
                    endif                                   
                ELSE                        !equation R(r)=0 has no real roots. we use the Legendre elliptic 
                                             !integrations and functions to compute the calculations. 
                    If(u.ne.zero)then 
                        call sncndn(p*w*sqrt(L1)+sign(one,f1234r)*pinf*w*sqrt(L1),one-m2,sn,cn,dn) 
                        index_p5(1)=-1
                        index_p5(2)=-2
                        index_p5(3)=2
                        index_p5(4)=-4
                        index_p5(5)=4
                        IF(f1234r.lt.zero)THEN                                 
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).    
                                y=u+(-two*u+w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = y 
                                pp=p  
                            else
                                y=rhorizon
                                r_coord = y
                                pp=PI0
                            endif         
                            x=robs 
                        ELSE 
                            if(p.lt.PI0)then
! equation (49) in Yang & Wang (2012).    
                                x=u+(-two*u-w*(L1-L2)*sn*abs(cn))/((L1-L2)*sn**2-(L1-one))
                                r_coord = x
                                pp=p
                            else
                                x=infinity 
                                r_coord = x
                                pp=PI0
                            endif
                            y=robs         
                        ENDIF 
                        !affine parameter part integration **********************************************
                        cases_int=5
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,a5,b5,index_p5,abs(pp),integ5,cases_int)
! equation (67) in Yang & Wang (2012).    
                        affr=integ5(5)
! equation (68) in Yang & Wang (2012).    
                        timer=two*integ5(3)+four*pp+affr
                        cases_int=2
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_add,b5,index_p5,abs(pp),integ5,cases_int)
                        call carlson_doublecomplex5(y,x,f1,g1,h1,f2,g2,h2,-r_m,b5,index_p5,abs(pp),integ15,cases_int)
                        !phy part**************************************************************************
! equation (68) in Yang & Wang (2012).    
                        timer=timer+Ap*integ5(2)-Am*integ15(2)
! equation (69) in Yang & Wang (2012).    
                        phyr=a_spin*(B_add*integ5(2)-B_m*integ15(2))        
                        IF(muobs.eq.zero.and.f1234t.eq.zero)THEN
! equation (18) in Yang & Wang (2012).    
                            phyr=phyr+pp*lambda        
                        ENDIF 
                    ELSE
                        If(f1234r.lt.zero)then         
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)-p*w)
                                r_coord = radius 
                            else
                                radius=rhorizon
                                r_coord = radius
                            endif
                            !timer part ******************************************************************        
                            y=radius
                            x=robs
                        ELSE 
                            If(p.lt.PI0)then
                                radius=w*tan(atan(robs/w)+p*w)
                                r_coord = radius
                            else
                                radius=infinity
                                r_coord = radius
                            endif
                            !timer part ************************************************************************8
                            y=robs
                            x=radius
                        ENDIF
                        pp_time=(x-y)+atan(x/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))-&
                                  atan(y/w)*(-w+four/w-r_add*Ap/w/(w**two+r_add**two)+r_m*Am/w/(w**two+r_m**two))
                        pp_time=pp_time+log(x**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two))-&
                                  (log(y**two+w**two)*(one-Ap/two/(w**two+r_add**two)+Am/two/(w**two+r_m**two)))        
                        timer=pp_time+Ap*log(abs(x-r_add))/(w**two+r_add**two)-Am*log(abs(x-r_m))/(w**two+r_m**two)-&
                                   (Ap*log(abs(y-r_add))/(w**two+r_add**two)-Am*log(abs(y-r_m))/(w**two+r_m**two)) 
                        !affine parameter part *****************************************************************
                        affr=(x-y)-w*atan(x/w)+w*atan(y/w)          
                        !phy part ******************************************************************
                        IF(a_spin .NE. zero)THEN
                            phyr=(-B_add*r_add/w/(r_add**two+w**two)+B_m*r_m/w/(r_m**two+w**two))&
                                               *(atan(x/w)-atan(y/w))+&
                                 log(abs(x-r_add)/sqrt(x**two+w**two))*B_add/(r_add*two+w**two)-&
                                 log(abs(y-r_add)/sqrt(y**two+w**two))*B_add/(r_add*two+w**two)-&
                                 log(abs(x-r_m)/sqrt(x**two+w**two))*B_m/(r_m*two+w**two)+&
                                 log(abs(y-r_m)/sqrt(y**two+w**two))*B_m/(r_m*two+w**two)
                            phyr=phyr*a_spin 
                        ELSE
                            phyr=zero
                        ENDIF        
                        If(muobs.eq.zero.and.f1234t.eq.zero)then         
                            phyr=phyr+lambda*(atan(x/w)-atan(y/w))/w
                        ENDIF   
                    ENDIF                        
                ENDIF
            ELSE
                count_num=1
                goto 40
            endif
        endif
           
        RETURN
     END SUBROUTINE INTRPART 

!********************************************************************************************
      Function  Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mu,rout,rin) 
!********************************************************************************************
!*     PURPOSE:  Solves equation \mu(p)=mu, i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.  
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are defined by equation 
!*                              (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon.
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.    
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.       
!*     REMARKS:                 This routine just search the intersection points of geodesic with 
!*                              up surface of disk. Following routine Pemdisk_all will searches 
!*                              intersection points of geodesic with up and down surface of disk.      
!*     ROUTINES CALLED: mutp, mu2p.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision Pemdisk,f1234(4),f3,f2,sinobs,muobs,a_spin,lambda,q,mu_tp,two,&
                         mu,scal,zero,robs,&
                         mu_tp2,four,one,pm,rout,rin,re
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals
        logical :: mobseqmtp 

        f3 = f1234(3)
        f2 = f1234(2) 
        mobseqmtp=.false. 
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp,mu_tp2,reals,mobseqmtp) 
 
        If(reals.eq.2)then
            If(mobseqmtp)then        
                t1=0
                t2=0
            else
                If(muobs.gt.zero)then
                    If(f2.lt.zero)then        
                        t1=1
                        t2=0
                    endif
                    If(f2.gt.zero)then        
                        t1=0
                        t2=0                                
                    endif
                else
                    If(muobs.eq.zero)then
                        If(f2.lt.zero)then        
                            t1=1
                            t2=0
                        endif
                        If(f2.gt.zero)then        
                            t1=0
                            t2=1                                
                        endif                        
                    else
                        If(f2.lt.zero)then        
                            t1=0
                            t2=0
                        endif
                        If(f2.gt.zero)then        
                            t1=0
                            t2=1                                
                        endif
                    endif        
                endif        
            endif
            !write(*,*)'mu=',mu,B,t1,t2
            pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal)  
            re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal)  
            IF(re .le. rout .and. re .ge. rin)THEN
                Pemdisk = pm
                return 
            ENDIF
            IF(re .gt. rout)THEN
                Pemdisk = -one 
                RETURN
            ENDIF 
            IF(re .lt. rin)THEN
                Pemdisk = -two
                return   
            ENDIF
        else
            Pemdisk=-one                 
        endif 
        return
      End Function Pemdisk 
!********************************************************************* 
      Function  Pemdisk_all(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mu,rout,rin) 
!********************************************************************* 
!*     PURPOSE:  Solves equation \mu(p)=mu, where \mu(p)=\mu(p), i.e. to search the value p_{em} of
!*               parameter p, corresponding to the intersection point of geodesic with 
!*               disk, i.e. a surface has a constants inclination angle with respect to 
!*               equatorial plane of black hole.
!* 
!*     INPUTS:   f1234(1:4)-----array of f_1, f_2, f_3, f_0, which are defined by equation 
!*                              (106)-(109) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               mu-------------mu=\cos(\pi/2-\theta_disk), where \theta_disk is the inclination 
!*                              angle of disk surface with respect to the equatorial plane of 
!*                              black hole.  
!*               rin, rout------inner and outer radius of disk.        
!*     OUTPUTS:  pemdisk--------value of root of equation \mu(p)= mu for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.        
!*     REMARKS:                 This routine will searches intersection points of 
!*                              geodesic with double surfaces of disk.      
!*     ROUTINES CALLED: mutp, mu2p, radius.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision Pemdisk_all,f3,f2,sinobs,muobs,a_spin,lambda,q,mu_tp1,two,&
                         mu,scal,robs,zero,&
                         mu_tp2,four,one,pm,f1234(4),rout,rin,re,rhorizon
        parameter(zero=0.D0,two=2.D0,four=4.D0,one=1.D0)
        integer  t1,t2,reals,i,j
        logical :: mobseqmtp 

        f3 = f1234(3)
        f2 = f1234(2) 
        mobseqmtp=.false.
        rhorizon = one + dsqrt(one-a_spin*a_spin) 
        call mutp(f2,f3,sinobs,muobs,a_spin,lambda,q,mu_tp1,mu_tp2,reals,mobseqmtp)         

        If(reals.eq.2)then
            Do j = 0,10
                Do i=j,j+1    
                    If(mobseqmtp)then        
                        If(muobs.eq.mu_tp1)then
                            t1=j
                            t2=i
                        else
                            t1=i
                            t2=j
                        endif
                    else
                        If(muobs.gt.zero)then
                            If(f2.lt.zero)then        
                                t1=i
                                t2=j
                            endif
                            If(f2.gt.zero)then        
                                t1=j
                                t2=i                                
                            endif
                        else
                            If(muobs.eq.zero)then
                                If(f2.lt.zero)then        
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then        
                                    t1=j
                                    t2=i                                
                                endif                        
                            else
                                If(f2.lt.zero)then        
                                    t1=i
                                    t2=j
                                endif
                                If(f2.gt.zero)then        
                                    t1=j
                                    t2=i                                
                                endif
                            endif        
                        endif        
                    endif
                    pm=mu2p(f3,f2,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,scal) 
                    IF(pm .le. zero)cycle
                    re = radius(pm,f1234(1),lambda,q,a_spin,robs,scal) 
                    !write(*,*)'mu=',re,pm,t1,t2
                    IF(re .le. rout .and. re .ge. rin)THEN
                        Pemdisk_all = pm
                        return
                    ELSE
                        IF(re .ge. infinity)THEN
                            Pemdisk_all = -one 
                            RETURN
                        ELSE
                            IF(re .le. rhorizon)THEN
                                Pemdisk_all = -two
                                return  
                            ENDIF
                        ENDIF 
                    ENDIF 
                ENDDO
            ENDDO
        else
            pm=-two                        
        endif
        Pemdisk_all=pm
        return
      End Function Pemdisk_all
!*****************************************************************************************************
      subroutine metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
!*****************************************************************************************************
!*     PURPOSE:  Computes Kerr metric, exp^\nu, exp^\psi, exp^mu1, exp^\mu2, and omiga at position:
!*               r_obs, \theta_{obs}.     
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).         
!*     OUTPUTS:  somiga,expnu,exppsi,expmu1,expmu2------------Kerr metrics under Boyer-Lindquist coordinates.
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision robs,a_spin,Delta,bigA,two,sinobs,muobs,one,sigma
        Double precision ,optional :: somiga,expnu,exppsi,expmu1,expmu2        

        two=2.D0        
        one=1.D0
! equations (1) and (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu=sqrt(sigma*Delta/bigA)
        exppsi=sinobs*sqrt(bigA/sigma)
        expmu1=sqrt(sigma/Delta)
        expmu2=sqrt(sigma)
        return        
      End subroutine metricg        

!********************************************************************************************
      Subroutine lambdaq_old(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using 
!*               formulae (110) and (112) in Yang & Wang (2012).    
!*     INPUTS:   alpha,beta-----Impact parameters.
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of f_1, f_2, f_3, f_4, which was defined by equation 
!*                              (106)-(109) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision f1234(4),robs,sinobs,muobs,a_spin,lambda,q,A1,&
                         Delta,zero,one,two,four,Sigma,at,Bt,&
                         scal,Vr,Vt,Vp,gama,gama_tp,gama_p,Aab,expnu2,&
                         eppsi2,epmu12,epmu22,bigA,three,alpha,beta,&
                         velocity(3),somiga
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)

        if(abs(beta).lt.1.D-7)beta=zero
        if(abs(alpha).lt.1.D-7)alpha=zero
! equations (113), (114) in Yang & Wang (2012).
        at=alpha/scal/robs
        Bt=beta/scal/robs        
        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (92) in Yang & Wang (2012).
        Aab=-one/dsqrt(one+at**two+Bt**two)
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/dsqrt(one-(Vt**two+Vp**two))
        gama_p=one/dsqrt(one-Vp**two)
! equation (106)-(109) in Yang & Wang (2012).
        f1234(1)=(-gama*Vr+gama/gama_tp*Aab)
        f1234(2)=-((-gama*Vt+gama*gama_tp*Vr*Vt*Aab)*robs*scal+gama_tp/gama_p*beta*Aab)
        f1234(3)=-((-gama*Vp+gama*gama_tp*Vr*Vp*Aab)*robs*scal+gama_tp*&
                                  gama_p*Vt*Vp*beta*Aab+gama_p*alpha*Aab)
        f1234(4)=((gama-gama*gama_tp*Vr*Aab)*robs*scal-gama_tp*gama_p*Vt*&
                                  beta*Aab-gama_p*Vp*alpha*Aab)/robs/scal
        If(dabs(f1234(1)).lt.1.D-7)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-7)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-7)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-7)f1234(4)=zero 
! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu2=Sigma*Delta/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma 
! equations (110), (112) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/robs/scal/f1234(4)*&
                    (one-lambda*somiga))**two*bigA/Delta 
        !write(*,*)'kk=',f1234(3),dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal,f1234(3)*somiga*sinobs
       return
       End subroutine lambdaq_old 

!********************************************************************************************
      Subroutine lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from impact parameters alpha and beta by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   alpha,beta-----Impact parameters.
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision f1234(4),robs,sinobs,muobs,a_spin,lambda,q,A1,&
                         Delta,zero,one,two,four,Sigma,at,Bt,&
                         scal,Vr,Vt,Vp,gama,expnu2,&
                         eppsi2,epmu12,epmu22,bigA,three,alpha,beta,&
                         velocity(3),somiga,prt,ptt,ppt
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)

        if(abs(beta).lt.1.D-7)beta=zero
        if(abs(alpha).lt.1.D-7)alpha=zero
! equations (94), (95) in Yang & Wang (2012).
        at=alpha/scal/robs
        Bt=beta/scal/robs        
        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (90) in Yang & Wang (2012). 
        gama=one/dsqrt(one-(Vr**two+Vt**two+Vp**two)) 

! equations (97), (98), (99) in Yang & Wang (2012). 
        prt=-one/dsqrt(one+at**two+Bt**two)
        ptt=Bt*prt
        ppt=at*prt
! equations (89), (90) and (91) in Yang & Wang (2012).
        f1234(1)=( gama*Vr-prt*(one+gama*gama*Vr*Vr/(one+gama))-&
                 ptt*gama*gama*Vr*Vt/(one+gama)-ppt*gama*gama*Vr*Vp/(one+gama) )*robs*scal 
        f1234(2)=( gama*Vt-prt*gama*gama*Vt*Vr/(one+gama)-ptt*(one+&
                 gama*gama*Vt*Vt/(one+gama))-ppt*gama*gama*Vt*Vp/(one+gama) )*robs*scal 
        f1234(3)=( gama*Vp-prt*gama*gama*Vp*Vr/(one+gama)-ptt*gama*gama*Vp*Vt/(one+gama)-&
                 ppt*(one+gama*gama*Vp*Vp/(one+gama)) )*robs*scal
        f1234(4)=gama*(one-prt*Vr-ptt*Vt-ppt*Vp) 

!in the above equations, f1234(1), f1234(2), f1234(3), f1234(4) denote pr, ptheta, pphi, pt
!a transformation was done from the source frame (denoted p_accent in the paper) with velocities vr, vt, vp to the LNRF frame 
!in the paper, see 5th line below eq. 89


! Keep r component p_r of four momentum to be negative, so the photon will go
! to the central black hole.
        f1234(1)=-f1234(1) 
        If(dabs(f1234(1)).lt.1.D-6)f1234(1)=zero
        If(dabs(f1234(2)).lt.1.D-6)f1234(2)=zero
        If(dabs(f1234(3)).lt.1.D-6)f1234(3)=zero
        If(dabs(f1234(4)).lt.1.D-6)f1234(4)=zero 
! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        expnu2=Sigma*Delta/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma 
! equations (86) and (87) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)/robs/scal*&
                    (one-lambda*somiga))**two*bigA/Delta  
        !write(*,*)'kk=',f1234(3)!,dsqrt(Delta)*Sigma/bigA*f1234(4)*robs*scal,f1234(3)*somiga*sinobs
       return
       End subroutine lambdaq

!********************************************************************************************
      Subroutine initialdirection(pr,ptheta,pphi,sinobs,&
                              muobs,a_spin,robs,velocity,lambda,q,f1234)
!********************************************************************************************
!*     PURPOSE:  Computes constants of motion from components of initial 4 momentum 
!*               of photon measured by emitter in its local rest frame, by using 
!*               formulae (86) and (87) in Yang & Wang (2012).    
!*     INPUTS:   pr,ptheta,pphi-----components of initial 4 momentum of photon measured by 
!*               emitter in its local rest frame.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or the initial position of photon.
!*               velocity(1:3)--Array of physical velocities of the observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  f1234(1:4)-----array of p_r, p_theta, p_phi, p_t, which are the components of 
!*                              four momentum of a photon measured under the LNRF frame, and 
!*                              defined by equations (82)-(85) in Yang & Wang (2012). 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  5 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision lambda,q,sinobs,muobs,a_spin,robs,zero,one,two,three,four,&
                        velocity(3),Vr,Vt,Vp,&
                        a1,gama,gama_tp,gama_p,f1234(4),&
                        Delta,Sigma,bigA,eppsi2,epmu12,epmu22,somiga,&
                        pr,ptheta,pphi
        parameter(zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)
        optional pr,ptheta,pphi 

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (92) in Yang & Wang (2012).
        gama=one/sqrt(one-(Vr**two+Vt**two+Vp**two))
        gama_tp=one/sqrt(one-(Vt**two+Vp**two))
        gama_p=one/sqrt(one-Vp**two)
! equation (106)-(109) in Yang & Wang (2012).
        f1234(1)=-(-gama*Vr+gama/gama_tp*pr)
        f1234(2)=-(-gama*Vt+gama*gama_tp*Vr*Vt*pr+gama_tp/gama_p*ptheta)
        f1234(3)=-(-gama*Vp+gama*gama_tp*Vr*Vp*pr+gama_tp*gama_p*Vt*Vp*ptheta+gama_p*pphi)
        f1234(4)=(gama-gama*gama_tp*Vr*pr-gama_tp*gama_p*Vt*ptheta-gama_p*Vp*pphi)

!*************************** added by Pieter ********************************************** codewoord kreeft
! equations (89), (90) and (91) in Yang & Wang (2012).
!        f1234(1)=( gama*Vr-pr*(one+gama*gama*Vr*Vr/(one+gama))-&
!                 pt*gama*gama*Vr*Vt/(one+gama)-pp*gama*gama*Vr*Vp/(one+gama) )
!        f1234(2)=( gama*Vt-pr*gama*gama*Vt*Vr/(one+gama)-pt*(one+&
!                 gama*gama*Vt*Vt/(one+gama))-pp*gama*gama*Vt*Vp/(one+gama) )
!        f1234(3)=( gama*Vp-pr*gama*gama*Vp*Vr/(one+gama)-pt*gama*gama*Vp*Vt/(one+gama)-&
!                 pp*(one+gama*gama*Vp*Vp/(one+gama)) )
!        f1234(4)=gama*(one-pr*Vr-pt*Vt-pp*Vp)
!        write(*,*) f1234(1), f1234(2), f1234(3), f1234(4)
!*************************** End of added by Pieter ***************************************

        If(abs(f1234(1)).lt.1.D-7)f1234(1)=zero
        If(abs(f1234(2)).lt.1.D-7)f1234(2)=zero
        If(abs(f1234(3)).lt.1.D-7)f1234(3)=zero
        If(abs(f1234(4)).lt.1.D-7)f1234(4)=zero

! equations (1), (2) in Yang & Wang (2012).
        Delta=robs**two-two*robs+a_spin**two
        Sigma=robs**two+(a_spin*muobs)**two
        bigA=(robs**two+a_spin**two)**two-(a_spin*sinobs)**two*Delta
        somiga=two*a_spin*robs/bigA
        eppsi2=sinobs**two*bigA/Sigma
        epmu12=Sigma/Delta
        epmu22=Sigma
! equations (86) and (87) in Yang & Wang (2012).
        A1 = f1234(3)/(dsqrt(Delta)*Sigma/bigA*f1234(4)+f1234(3)*somiga*sinobs) 
        lambda = A1*sinobs 
        !write(*,*) lambda, A1, sinobs
        !write(*,*) 'pr, ptheta, pphi', pr, ptheta, pphi
        !write(*,*) 'fr, ftheta, fphi', f1234(1), f1234(2), f1234(3) 
        q=(A1*A1-a_spin*a_spin)*muobs*muobs+(f1234(2)/f1234(4)*&
                   (one-lambda*somiga))**two*bigA/Delta   
      return
      End subroutine initialdirection

!********************************************************************************************        
      Subroutine center_of_image(robs,scal,velocity,alphac,betac)
!********************************************************************************************
!*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (100)  
!*               and (101) in Yang & Wang (2012). alphac, betac are the coordinates of 
!*               center point of images on the screen of observer.    
!*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.  
!*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
!*                              LNRF.        
!*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
!*     ROUTINES CALLED: NONE.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  9 Jan 2012.
!*     REVISIONS: ****************************************** 
        implicit none
        Double precision robs,scal,velocity(3),alphac,betac,zero,one,two,four,&
                        Vr,Vt,Vp,gama,a1,b1,c1,alphap,alpham,betap,betam
        parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

        Vr=velocity(1)
        Vt=velocity(2)
        Vp=velocity(3)
! equation (90) in Yang & Wang (2012).
        gama=one/sqrt(one-(Vr*Vr+Vt*Vt+Vp*Vp))       

        If(Vt.ne.zero)then
            If(Vp.ne.zero)then   
                a1=(one+gama*gama*(Vt*Vt+Vp*Vp)/(one+gama))**two-gama*gama*(Vp*Vp+Vt*Vt)  
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*(Vt*Vt+Vp*Vp))/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt
                betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1                         
                If(betap*Vp.lt.zero)then
                    betac=betap
                Else
                    betac=betam                
                Endif
                alphac=Vp/Vt*betac
            Else
                alphac=zero
                a1=(one+gama*gama*Vt*Vt/(one+gama))**two-gama*gama*Vt*Vt
                b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*Vt*Vt)/(one+gama)**two
                c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt    
                betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
                betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
                If(betap*Vt.lt.zero)then
                    betac=betap
                Else
                    betac=betam                
                Endif
            Endif  
        else
            betac=zero        
            If(Vp.ne.zero)then
                a1=(one+gama*gama*Vp*Vp/(one+gama))**two-gama*gama*Vp*Vp
                b1=two*gama*gama*Vp*Vr*(one+gama+gama*gama*Vp*Vp)/(one+gama)**two
                c1=(gama*gama*Vp*Vr/(one+gama))**two-gama*gama*Vp*Vp   
                alphap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1        
                alpham=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
                If(alphap*Vp.lt.zero)then
                    alphac=alphap        
                Else
                    alphac=alpham        
                Endif
            Else
                alphac=zero
            Endif
        endif
        alphac=alphac*robs*scal
        betac=betac*robs*scal
        return
      End Subroutine center_of_image
!********************************************************************************************        
!       Subroutine center_of_image(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
! !********************************************************************************************
! !*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (100)  
! !*               and (101) in Yang & Wang (2012). alphac, betac are the coordinates of 
! !*               center point of images on the screen of observer.    
! !*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
! !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
! !*                              \theta_{obs} is the inclination angle of the observer.
! !*               a_spin---------spin of black hole, on interval (-1,1).  
! !*               scal-----------a dimentionless parameter to control the size of the images.
! !*                              Which is usually be set to 1.D0.  
! !*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
! !*                              LNRF.        
! !*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
! !*     ROUTINES CALLED: NONE.
! !*     ACCURACY:   Machine.    
! !*     AUTHOR:     Yang & Wang (2012).  
! !*     DATE WRITTEN:  9 Jan 2012.
! !*     REVISIONS: ****************************************** 
!         implicit none
!         Double precision robs,sinobs,muobs,a_spin,scal,velocity(3),alphac,betac,zero,one,two,four,&
!                         Vr,Vt,Vp,gama,a1,b1,c1,alphap,alpham,betap,betam
!         parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

!         Vr=velocity(1)
!         Vt=velocity(2)
!         Vp=velocity(3)
! ! equation (90) in Yang & Wang (2012).
!         gama=one/sqrt(one-(Vr*Vr+Vt*Vt+Vp*Vp))       

!         If(Vt.ne.zero)then
!             If(Vp.ne.zero)then   
!                 a1=(one+gama*gama*(Vt*Vt+Vp*Vp)/(one+gama))**two-gama*gama*(Vp*Vp+Vt*Vt)  
!                 b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*(Vt*Vt+Vp*Vp))/(one+gama)**two
!                 c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt
!                 betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
!                 betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1                         
!                 If(betap*Vp.lt.zero)then
!                     betac=betap
!                 Else
!                     betac=betam                
!                 Endif
!                 alphac=Vp/Vt*betac
!             Else
!                 alphac=zero
!                 a1=(one+gama*gama*Vt*Vt/(one+gama))**two-gama*gama*Vt*Vt
!                 b1=two*gama*gama*Vt*Vr*(one+gama+gama*gama*Vt*Vt)/(one+gama)**two
!                 c1=(gama*gama*Vt*Vr/(one+gama))**two-gama*gama*Vt*Vt    
!                 betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
!                 betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
!                 If(betap*Vt.lt.zero)then
!                     betac=betap
!                 Else
!                     betac=betam                
!                 Endif
!             Endif  
!         else
!             betac=zero        
!             If(Vp.ne.zero)then
!                 a1=(one+gama*gama*Vp*Vp/(one+gama))**two-gama*gama*Vp*Vp
!                 b1=two*gama*gama*Vp*Vr*(one+gama+gama*gama*Vp*Vp)/(one+gama)**two
!                 c1=(gama*gama*Vp*Vr/(one+gama))**two-gama*gama*Vp*Vp   
!                 alphap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1        
!                 alpham=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
!                 If(alphap*Vp.lt.zero)then
!                     alphac=alphap        
!                 Else
!                     alphac=alpham        
!                 Endif
!             Else
!                 alphac=zero
!             Endif
!         endif
!         alphac=alphac*robs*scal
!         betac=betac*robs*scal
!         return
!       End Subroutine center_of_image
!********************************************************************************************        
!********************************************************************************************        
!       Subroutine center_of_image_old(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
! !********************************************************************************************
! !*     PURPOSE:  Solves equations f_3(alphac,betac)=0, f_2(alphac,betac)=0, of (119)  
! !*               and (120) in Yang & Wang (2012). alphac, betac are the coordinates of 
! !*               center point of images on the screen of observer.    
! !*     INPUTS:   robs-----------radial coordinate of observer or the initial position of photon.
! !*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
! !*                              \theta_{obs} is the inclination angle of the observer.
! !*               a_spin---------spin of black hole, on interval (-1,1).  
! !*               scal-----------a dimentionless parameter to control the size of the images.
! !*                              Which is usually be set to 1.D0.  
! !*               velocity(1:3)--Array of physical velocity of observer or emitter with respect to
! !*                              LNRF.        
! !*     OUTPUTS:  alphac,betac---coordinates of center point of images on the screen of observer.            
! !*     ROUTINES CALLED: NONE.
! !*     ACCURACY:   Machine.    
! !*     AUTHOR:     Yang & Wang (2012).  
! !*     DATE WRITTEN:  9 Jan 2012.
! !*     REVISIONS: ****************************************** 
!         implicit none
!         Double precision robs,sinobs,muobs,a_spin,scal,velocity(3),alphac,betac,zero,one,two,four,&
!                         Vr,Vt,Vp,gama,gama_tp,gama_p,a1,b1,c1,alphap,alpham,betap,betam
!         parameter(zero=0.D0,one=1.D0,two=2.D0,four=4.D0)

!         Vr=velocity(1)
!         Vt=velocity(2)
!         Vp=velocity(3)
! ! equation (92) in Yang & Wang (2012).
!         gama=one/sqrt(one-(Vr**two+Vt**two+Vp**two))
!         gama_tp=one/sqrt(one-(Vt**two+Vp**two))
!         gama_p=one/sqrt(one-Vp**two)        

!         If(Vt.ne.zero)then
!           If(Vp.ne.zero)then         
!             a1=(gama_tp/gama_p)**two-((Vp/gama_tp)**two+Vt**two)*gama**two
!             b1=two*gama*gama_tp*Vr*Vp/gama_p
!             c1=-Vp**two
!             alphap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
!             alpham=(-b1-sqrt(b1**two-four*a1*c1))/two/a1                         
!             If(alphap*Vp.lt.zero)then
!                 alphac=alphap
!             Else
!                 alphac=alpham                
!             Endif
!             betac=alphac*Vt*gama_tp/Vp
!           Else
!             alphac=zero
!             a1=one-(gama*Vt/gama_tp)**two
!             b1=two*gama*Vr*Vt
!             c1=(gama*Vt)**two*(Vr**two-one)        
!             betap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1
!             betam=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
!             If(betap*Vt.lt.zero)then
!                 betac=betap
!             Else
!                 betac=betam                
!             Endif
!           Endif  
!         else
!             betac=zero        
!             If(Vp.ne.zero)then
!                 a1=gama_p**two-(gama*Vp)**two
!                 b1=two*gama*gama_tp*gama_p*Vr*Vp
!                 c1=-(gama_tp*Vp)**two
!                 alphap=(-b1+sqrt(b1**two-four*a1*c1))/two/a1        
!                 alpham=(-b1-sqrt(b1**two-four*a1*c1))/two/a1        
!                 If(alphap*Vp.lt.zero)then
!                     alphac=alphap        
!                 Else
!                     alphac=alpham        
!                 Endif
!             Else
!                 alphac=zero
!             Endif
!         endif
!         alphac=alphac*robs*scal
!         betac=betac*robs*scal
!         return
!       End Subroutine center_of_image_old

!********************************************************************************************
      FUNCTION p_total(f1234r,lambda,q,sinobs,muobs,a_spin,robs,scal)
!******************************************************************************************** 
!*     PURPOSE:  Computes the integral value of \int^r dr (R)^{-1/2}, from the starting position to
!*               the termination----either the infinity or the event horizon.   
!*     INPUTS:   f1234r---------p_r, the r component of four momentum of the photon measured 
!*                              under the LNRF, see equation (83) in Yang & Wang (2012).
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0.        
!*     OUTPUTS:  p_total--------which is the value of integrals \int^r dr (R)^{-1/2}, along a 
!*                              whole geodesic, that is from the starting position to either go to
!*                              infinity or fall in to black hole.          
!*     ROUTINES CALLED: root3, weierstrass_int_J3, radiustp, weierstrassP, EllipticF, carlson_doublecomplex5 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************   
        IMPLICIT NONE
        DOUBLE PRECISION sinobs,muobs,a_spin,rhorizon,q,lambda,&
                         cc,b0,b1,b2,b3,g2,g3,tobs,pp,p1,p2,PI0,&
                         u,v,w,L1,L2,thorizon,m2,pinf,r_add,r_m,&
                         a4,b4,robs,&
                         scal,tinf,integ04(4),integ14(4),&
                         r_tp1,r_tp2,t_inf,tp2,f1234r,&
                         PI01,rff_p,p_total,p_tp1_tp2,PI2_p,PI1_p,sqt3         
        !PARAMETER(zero=0.D0,one=1.D0,two=2.D0,four=4.D0,three=3.D0)
        COMPLEX*16 bb(1:4),dd(3)
        INTEGER ::  reals,t1,t2,index_p4(4),del,cases_int,cases
        LOGICAL :: robs_eq_rtp,indrhorizon                

                rhorizon=one+sqrt(one-a_spin**two)
! equation (64) in Yang & Wang (2012).
                r_add=rhorizon
                r_m=one-sqrt(one-a_spin**two)      
! equation (64) in Yang & Wang (2012).   
                b4=one
                a4=zero
                cc=a_spin**2-lambda**2-q
                robs_eq_rtp=.false.
                indrhorizon=.false.
                call radiustp(f1234r,a_spin,robs,lambda,q,r_tp1,r_tp2,&
                                  reals,robs_eq_rtp,indrhorizon,cases,bb) 
!** R(r)=0 has real roots and turning points exists in radial r.
                If(reals.ne.0)then  
! equations (35)-(38) in Yang & Wang (2012).
                    b0=four*r_tp1**3+two*(a_spin**2-lambda**2-q)*r_tp1+two*(q+(lambda-a_spin)**2)
                    b1=two*r_tp1**2+one/three*(a_spin**2-lambda**2-q)
                    b2=four/three*r_tp1
                    b3=one
                    g2=three/four*(b1**2-b0*b2)
                    g3=one/16.D0*(three*b0*b1*b2-two*b1**three-b0**two*b3)
! equation (39) in Yang & Wang (2012).
                    If(robs-r_tp1.ne.zero)then        
                        tobs=b0/four/(robs-r_tp1)+b1/four
                    else
                        tobs=infinity
                    endif 
                    If(rhorizon-r_tp1.ne.zero)then
                        thorizon=b1/four+b0/four/(rhorizon-r_tp1)
                    else
                        thorizon=infinity         
                    endif
                    tp2=b0/four/(r_tp2-r_tp1)+b1/four   
                    tinf=b1/four       
! equation (64), (66) and (70) in Yang & Wang (2012).
                    call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)       

                    index_p4(1)=0
                    cases_int=1
                    call weierstrass_int_J3(tobs,infinity,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int) 
! equation (42) in Yang & Wang (2012).
                    PI0=integ04(1)   
                    select case(cases)
                    CASE(1)
                        If(f1234r .ge. zero)then !**photon will goto infinity.
                            index_p4(1)=0
                            cases_int=1
                            call weierstrass_int_J3(tinf,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                            p_total = integ04(1)      
                        ELSE 
                            If(.not.indrhorizon)then
                                index_p4(1)=0
                                cases_int=1 
                                call weierstrass_int_j3(tinf,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                                p_total = PI0+integ14(1)      
                            ELSE     !f1234r<0, photon will fall into black hole unless something encountered.        
                                index_p4(1)=0                
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                                p_total = integ04(1)  
                            ENDIF
                        ENDIF         
                    CASE(2)
                        If(.not.indrhorizon)then
                            write(*,*)'we come here!'  
                            If(f1234r.lt.zero)then
                                PI01=-PI0
                            else
                                PI01=PI0        
                            endif
! equation (41) in Yang & Wang (2012).    
                            index_p4(1)=0
                            cases_int=1         
                            call weierstrass_int_J3(tp2,infinity,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int)
                            pp=zero
! equation (57) in Yang & Wang (2012).
                            p_tp1_tp2=integ14(1) 
                            PI2_p=p_tp1_tp2-PI0 
                            PI1_p=PI0 
                            p1=PI0-pp
                            p2=p_tp1_tp2-p1   
! equation (58) in Yang & Wang (2012). 
                            t1 = 2
                            t2 = 2
                            If(robs_eq_rtp)then         
                                p_total=abs(pp)+two*(t1*p1+t2*p2)
                            else
                                If(f1234r.gt.zero)then 
                                    p_total=-pp+two*(t1*p1+t2*p2)
                                endif
                                If(f1234r.lt.zero)then 
                                    p_total=pp+two*(t1*p1+t2*p2)
                                endif
                            endif   
                        !*************************************************************************************
                        200     continue
                        else  !photon has probability to fall into black hole.
                            If(f1234r.le.zero)then
                                index_p4(1)=0
                                cases_int=1
                                call weierstrass_int_J3(tobs,thorizon,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)
                                p_total = integ04(1) 
                            ELSE  !p_r>0, photon will meet the r_tp2 turning point and turn around then goto vevnt horizon.     
                                index_p4(1)=0
                                cases_int=1        
                                call weierstrass_int_J3(tp2,tobs,dd,del,a4,b4,index_p4,rff_p,integ04,cases_int)  
                                call weierstrass_int_j3(tp2,thorizon,dd,del,a4,b4,index_p4,rff_p,integ14,cases_int) 
                                p_total = integ14(1)+integ04(1)        
                            ENDIF
                        ENDIF                             
                    END SELECT       
!************************************************************************************************ 
                    If(a_spin.eq.zero)then
                        If(cc.eq.zero)then
                            If(f1234r.lt.zero)then
                                p_total = one/rhorizon-one/robs
                            else
                                p_total = one/robs                       
                            endif        
                        endif
                        If(cc.eq.-27.D0)then
                            sqt3=sqrt(three)        
                            If(f1234r.lt.zero)then
                                p_total=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -Log(one+two/sqt3)/three/sqt3 
                            else        
                                p_total=Log(abs((sqrt(robs*(robs+6.D0))+(three+two*robs)/sqt3)/(robs-three)))/three/sqt3&
                                                -Log(one+two/sqt3)/three/sqt3                     
                            endif                
                        endif
                    endif                                   
                ELSE            
! equation (44) in Yang & Wang (2012).   !equation R(r)=0 has no real roots. we use the Legendre elliptic 
                    u=real(bb(4))        !integrations and functions to compute the calculations.
                    w=abs(aimag(bb(4)))
                    v=abs(aimag(bb(2)))
                    If(u.ne.zero)then
! equation (45) in Yang & Wang (2012).       
                        L1=(four*u**2+w**2+v**2+sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
                        L2=(four*u**2+w**2+v**2-sqrt((four*u**2+w**2+v**2)**2-four*w**2*v**2))/(two*w**2)
! equation (46) in Yang & Wang (2012).       
                        thorizon=sqrt((L1-one)/(L1-L2))*(rhorizon-u*(L1+one)/(L1-one))/sqrt((rhorizon-u)**2+w**2)
! equation (48) in Yang & Wang (2012).   
                        m2=(L1-L2)/L1
                        tinf=sqrt((L1-one)/(L1-L2))*(robs-u*(L1+one)/(L1-one))/sqrt((robs-u)**two+w**two)
                        t_inf=sqrt((L1-one)/(L1-L2))
! equation (50) in Yang & Wang (2012).       
                        pinf=EllipticF(tinf,m2)/(w*sqrt(L1))
                        IF(f1234r.lt.zero)THEN
                            p_total = pinf-EllipticF(thorizon,m2)/(w*sqrt(L1))             
                        ELSE
                            p_total = EllipticF(t_inf,m2)/(w*sqrt(L1))-pinf        
                        ENDIF  
                    ELSE
                        If(f1234r.lt.zero)then
                            p_total = (atan(robs/w)-atan(rhorizon/w))/w         
                        ELSE
                            p_total = (PI/two-atan(robs/w))/w         
                        ENDIF   
                    ENDIF                        
                ENDIF             
        RETURN
     END FUNCTION p_total

!********************************************************************************************        
      end module BLcoordinate 


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      module sphericalmotion
!********************************************************************
!* This module aim on the calculation of the geodesic trajectory of a photon
!* doing spherical motion.
!********************************************************************
      USE blcoordinate
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      contains 
!********************************************************************************************
      SUBROUTINE SPHERICALMOTION_BL(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                               robs,thetamax,phyt,timet,sigmat,mucos,t1,t2)    
!******************************************************************************************** 
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{ini}), muobs=cos(\theta_{ini}), where 
!*                              \theta_{ini} is the initial theta angle of the photon.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               thetamax-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      IMPLICIT NONE
      Double precision phyt,timet,kp,kt,p,sinobs,muobs,a_spin,lambda,q,mu_tp1,tp2,tmu,&
             b0,b1,b2,b3,g2,g3,tobs,p1,p2,pp,c_add,c_m,a_add,a_m,come,&
             a4,b4,delta,mu_tp2,robs,integ5(5),integ(5),rff_p,tp1,&
             integ15(5),PI0,integ05(5),p_mu,PI01,h,p1_t,p2_t,pp_t,p1_phi,&
             p2_phi,pp_phi,mucos,p_mt1_mt2,PI1_sig,PI2_sig,&
             PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,p1_sig,p2_sig,pp_sig,sigmat 
      Double precision kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,&
             c_phi,c_time,thetamax,thetamax_1,sinmax,mumax,Omega,somiga,&
             expnu,exppsi,expmu1,expmu2,c_tau 
      integer ::  t1,t2,i,j,index_p5(5),del,cases_int,count_num=1
      complex*16 dd(3)
      logical :: mobseqmtp
      save  kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,a4,b4,mu_tp1,mu_tp2,&
            mobseqmtp,b0,b1,b2,b3,g2,g3,dd,del,PI0,c_m,c_add,a_m,a_add,tp2,tobs,h,p_mt1_mt2,&
            PI1_phi,PI2_phi,PI1_time,PI2_time,PI2_p,come,tp1,Delta,c_phi,c_time,&
            PI01,sinmax,mumax,thetamax_1,Omega,c_tau

30    continue
      IF(count_num.eq.1)then
          kp_1=kp
          kt_1=kt
          lambda_1=lambda
          q_1=q
          muobs_1=muobs
          sinobs_1=sinobs
          a_spin_1=a_spin
          robs_1=robs
          thetamax_1 = thetamax
          t1=0
          t2=0  
          IF(thetamax.ne.90.D0 .and. thetamax.ne.180.D0)THEN
              sinmax = dsin(thetamax*dtor)
              mumax = dcos(thetamax*dtor)
          Else
              IF(thetamax.eq.90.D0)THEN
                  sinmax = one
                  mumax = zero
              ENDIF
              IF(thetamax.eq.180.D0)THEN
                  sinmax = zero
                  mumax = one
              ENDIF
          ENDIF
          mu_tp1 = abs(mumax)
          mu_tp2 = -mu_tp2
     !***********************************************************************    
          If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, so the integrations about \theta are valished.
              Omega = one/(robs**(three/two)+a_spin)
              mucos=zero
              phyt=p
              timet=phyt/Omega
              call metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
              c_tau = dsqrt(-expnu*expnu+somiga*somiga*exppsi*exppsi-two*exppsi*somiga*Omega+&
                            exppsi*exppsi*Omega*Omega)
              sigmat=c_tau*timet
              count_num=count_num+1 
              return
          endif
     !**************************************************************************
          If(a_spin.eq.zero)then
              timet=zero
              CALL SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                           robs,thetamax,phyt,timet,sigmat,mucos,t1,t2) 
              count_num=count_num+1
              return
          endif
          a4=zero
          b4=one
!      equations (26)-(29) of Yang and Wang (2013).
          come = -one  
          b0=four*a_spin**2*mu_tp1**3*come-two*mu_tp1*(a_spin**2*come+lambda**2+q)
          b1=two*a_spin**2*mu_tp1**2*come-(a_spin**2*come+lambda**2+q)/three
          b2=four/three*a_spin**2*mu_tp1*come
          b3=a_spin**2*come
!      equations (31) of Yang and Wang (2013).
          g2=three/four*(b1**2-b0*b2)
          g3=(three*b0*b1*b2-two*b1**3-b0**2*b3)/16.D0  
          call root3(zero,-g2/four,-g3/four,dd(1),dd(2),dd(3),del)

!      equations (30) of Yang and Wang (2013).
          If(muobs.ne.mu_tp1)then
              tobs=b0/four/(muobs-mu_tp1)+b1/four
          else
              tobs=infinity
          endif
          tp1=infinity
          tp2=b0/four/(mu_tp2-mu_tp1)+b1/four
!      equations (72)-(73) of Yang and Wang (2013).
          If(mu_tp1-one.ne.zero)then
               c_m=b0/(four*(-one-mu_tp1)**2)
             c_add=b0/(four*(one-mu_tp1)**2) 
               a_m=b0/four/(-one-mu_tp1)+b1/four
             a_add=b0/four/(one-mu_tp1)+b1/four
          endif
          index_p5(1)=0
          cases_int=1
!      equations (53) of Yang and Wang (2013).
          call weierstrass_int_J3(tobs,tp1,dd,del,a4,b4,index_p5,rff_p,integ05,cases_int)
          PI0=integ05(1)
          If(kt.lt.zero)then
              PI01=-PI0
          else
              PI01=PI0
          endif
          tmu=weierstrassP(p+PI01,g2,g3,dd,del)
!      equations (32) of Yang and Wang (2013).
          mucos = mu_tp1+b0/(four*tmu-b1)
          h=-b1/four
          !to get number of turn points of t1 and t2.
          !111111111*****************************************************************************************
          !mu=mu_tp+b0/(four*tmu-b1)
          call weierstrass_int_J3(tp2,tp1,dd,del,a4,b4,index_p5,rff_p,integ15,cases_int)
          call weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int)
!      equations (53) of Yang and Wang (2013).	
          p_mt1_mt2=integ15(1)
          PI2_p=p_mt1_mt2-PI0
          pp=integ5(1)
          p1=PI0-pp
!      equations (53) of Yang and Wang (2013).
          p2=p_mt1_mt2-p1
          PI1_phi=zero
          PI2_phi=zero
          PI1_sig=zero
          PI2_sig=zero
          PI1_time=zero
          PI2_time=zero 
!      equations (54) of Yang and Wang (2013).
          Do j=0,10**6
              Do i=j,j+1
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      else
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                  else
                      If(kt.lt.zero)then
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                      If(kt.gt.zero)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      endif
                  endif    
                  If(dabs(p-p_mu).lt.1.D-3)goto 400
              enddo
          enddo
          !11111111*************************************************************** 
400       continue
!      equations (71)-(73) of Yang and Wang (2013).
          Delta=robs*robs+a_spin*a_spin-two*robs
          c_phi = a_spin*(robs*(two)-(lambda*a_spin))/Delta
          c_time = ( (two)*robs**three+(two*a_spin*(a_spin-lambda))*robs )/Delta
          index_p5(1)=-1
          index_p5(2)=-2
          index_p5(3)=0
          index_p5(4)=-4
          index_p5(5)=0
          !*****pp part***************************************
          If(lambda.ne.zero)then
              cases_int=2
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5,abs(pp),integ5,cases_int)
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5,abs(pp),integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
              pp_phi=(pp/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*pp 
          else 
              pp_phi=c_phi*pp
          endif
          cases_int=4
          call weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5,abs(pp),integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
          pp_sig=(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                  sixteen)*a_spin*a_spin+robs*robs*pp
!      equations (22) of Yang and Wang (2013).
          pp_t=pp_sig+c_time*pp  
          !*****p1 part***************************************
          If(t1.eq.0)then
              p1_phi=zero
              p1_sig=zero 
              p1_t=zero
          else  
              If(lambda.ne.zero)then  
                  IF(PI1_phi .EQ. zero)THEN
                      cases_int=2
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_add,b4,index_p5,PI0,integ5,cases_int)
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_m,b4,index_p5,PI0,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI1_phi=(PI0/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI0
                  ENDIF 
                  p1_phi=PI1_phi-pp_phi 
              else 
                  IF(PI1_phi.eq.zero)PI1_phi=c_phi*PI0
                  p1_phi=PI1_phi-pp_phi
              endif 
              IF(PI1_time .EQ. zero .or. PI1_sig.eq.zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI1_sig=(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                           sixteen)*a_spin*a_spin+robs*robs*PI0
!      equations (22) of Yang and Wang (2013).
                  PI1_time=PI1_sig+c_time*PI0  
              ENDIF
              p1_sig=PI1_sig-pp_sig
              p1_t=PI1_time-pp_t 
          endif 
          !*****p2 part***************************************
          If(t2.eq.0)then
              p2_phi=zero
              p2_t=zero
          else
              IF(lambda.ne.zero)then  
                  IF(PI2_phi .EQ. zero)THEN  
                      cases_int=2
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,index_p5,PI2_p,integ5,cases_int)
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,index_p5,PI2_p,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI2_phi=(PI2_p/(one-mu_tp1*mu_tp1)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI2_p
                  ENDIF
                  p2_phi=PI2_phi+pp_phi
              ELSE
                  IF(PI2_phi.eq.zero)PI2_phi=c_phi*PI2_p  
                  p2_phi=PI2_phi+pp_phi             
              ENDIF 
                
              IF(PI2_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI2_sig=(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI2_p
!      equations (22) of Yang and Wang (2013).
                  PI2_time=PI2_sig+c_time*PI2_p  
              ENDIF
              p2_sig=PI2_sig+pp_sig
              p2_t=PI2_time+pp_t   
          endif   
          !write(*,*)pp_phi,p1_phi,p2_phi,t1,t2
	!**************************************************************
          If(mobseqmtp)then 
              If(muobs .eq. mu_tp1)then
!      equations (52) of Yang and Wang (2013).
                  phyt= -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet= -pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat= -pp_sig+two*(t1*p1_sig+t2*p2_sig) 
              else
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)  
              endif
          else
              If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)
              endif
              If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).		
                  phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=-pp_sig+two*(t1*p1_sig+t2*p2_sig)
               endif
          endif
          IF(mu_tp1.eq.one)phyt = phyt+(t1+t2)*PI
          !phyt = mod(phyt,twopi)
          !If(phyt .lt. zero)phyt=phyt+twopi
          count_num=count_num+1
      ELSE  
          If(kp_1.eq.kp.and.kt_1.eq.kt.and.lambda_1.eq.lambda.and.q_1.eq.q.and.&
          sinobs_1.eq.sinobs.and.muobs_1.eq.muobs.and.a_spin_1.eq.a_spin.and.robs_1.eq.robs&
          .and.thetamax_1.eq.thetamax)then   
	!***************************************************************************
          t1=0
          t2=0  
        !***********************************************************************   	
          If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, so the integrations about \theta are valished. 
              mucos=zero
              phyt=p
              timet=phyt/Omega 
              sigmat=c_tau*timet 
              return
          endif
          !**************************************************************************
          If(a_spin.eq.zero)then 
              CALL SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,&
                           robs,thetamax,phyt,timet,sigmat,mucos,t1,t2) 
              return
          endif   
          tmu=weierstrassP(p+PI01,g2,g3,dd,del)
          mucos = mu_tp1+b0/(four*tmu-b1)
          !to get number of turn points of t1 and t2.
          !111111111*****************************************************************************************
          !mu=mu_tp+b0/(four*tmu-b1) 
          index_p5(1)=0
          cases_int=1
          call weierstrass_int_J3(tobs,tmu,dd,del,a4,b4,index_p5,rff_p,integ5,cases_int)  
          pp=integ5(1)
          p1=PI0-pp
          p2=p_mt1_mt2-p1 
          Do j=0,10**6
              Do i=j,j+1
                  If(mobseqmtp)then
                      If(muobs.eq.mu_tp1)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      else
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                  else
                      If(kt.lt.zero)then
                          t1=i
                          t2=j
                          p_mu=pp+two*(t1*p1+t2*p2)
                      endif
                      If(kt.gt.zero)then
                          t1=j
                          t2=i
                          p_mu=-pp+two*(t1*p1+t2*p2)
                      endif
                  endif  
                  !write(*,*)p,p_mu,abs(p-p_mu),t1,t2!(a,B,lambda,q,mu,sinobs,muobs,a_spin,t1,t2,robs),t1,t2
                  If(abs(p-p_mu).lt.1.D-3)goto 410
              enddo
          enddo
          !11111111*************************************************************** 
410       continue
          index_p5(1)=-1
          index_p5(2)=-2
          index_p5(3)=0
          index_p5(4)=-4
          index_p5(5)=0
          !*****pp part***************************************
          If(lambda.ne.zero)then
              cases_int=2
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_add,b4,index_p5,abs(pp),integ5,cases_int)
              call weierstrass_int_J3(tobs,tmu,dd,del,-a_m,b4,index_p5,abs(pp),integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
              pp_phi=(pp/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*pp 
          else 
              pp_phi=c_phi*pp
          endif
          cases_int=4
          call weierstrass_int_J3(tobs,tmu,dd,del,h,b4,index_p5,abs(pp),integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
          pp_sig=(pp*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                    sixteen)*a_spin*a_spin+robs*robs*pp
!      equations (22) of Yang and Wang (2013).
          pp_t=pp_sig+c_time*pp
          !*****p1 part***************************************
          If(t1.eq.0)then
              p1_phi=zero
              p1_t=zero
          else  
              If(lambda.ne.zero)then  
                  IF(PI1_phi .EQ. zero)THEN
                      cases_int=2
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_add,b4,index_p5,PI0,integ5,cases_int)
                      call weierstrass_int_J3(tobs,infinity,dd,del,-a_m,b4,index_p5,PI0,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI1_phi=(PI0/(one-mu_tp1**two)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI0
                  ENDIF 
                  p1_phi=PI1_phi-pp_phi 
              else 
                  IF(PI1_phi.eq.zero)PI1_phi=c_phi*PI0
                  p1_phi=PI1_phi-pp_phi     
              endif 
              IF(PI1_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tobs,infinity,dd,del,h,b4,index_p5,PI0,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI1_sig=(PI0*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI0
!      equations (22) of Yang and Wang (2013).
                  PI1_time=PI1_sig+c_time*PI0  
              ENDIF
              p1_sig=PI1_sig-pp_sig
              p1_t=PI1_time-pp_t 
          endif 
          !*****p2 part***************************************
          If(t2.eq.0)then
              p2_phi=zero
              p2_t=zero
          else
              IF(lambda.ne.zero)then  
                  IF(PI2_phi .EQ. zero)THEN  
                      cases_int=2
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_add,b4,index_p5,PI2_p,integ5,cases_int)
                      call weierstrass_int_J3(tp2,tobs,dd,del,-a_m,b4,index_p5,PI2_p,integ15,cases_int)
!      equations (21) (72) of Yang and Wang (2013).
                      PI2_phi=(PI2_p/(one-mu_tp1*mu_tp1)+(integ5(2)*c_add-integ15(2)*c_m)/two)*lambda+c_phi*PI2_p 
                  ENDIF
                  p2_phi=PI2_phi+pp_phi
              ELSE
                  IF(PI2_phi.eq.zero)PI2_phi=c_phi*PI2_p 
                  p2_phi=PI2_phi+pp_phi           
              ENDIF 
                
              IF(PI2_time .EQ. zero)THEN  
                  cases_int=4  
                  call weierstrass_int_J3(tp2,tobs,dd,del,h,b4,index_p5,PI2_p,integ,cases_int)
!      equations (20) (71) of Yang and Wang (2013).
                  PI2_sig=(PI2_p*mu_tp1**two+integ(2)*mu_tp1*b0/two+integ(4)*b0**two/&
                                  sixteen)*a_spin*a_spin+robs*robs*PI2_p
!      equations (22) of Yang and Wang (2013).
                  PI2_time=PI2_sig+c_time*PI2_p   
              ENDIF
              p2_sig=PI2_sig+pp_sig 
              p2_t=PI2_time+pp_t   
          endif   
          !write(*,*)'kkk=',pp_phi,p1_phi,p2_phi,t1,t2
	!**************************************************************
          If(mobseqmtp)then 
              If(muobs .eq. mu_tp1)then
!      equations (52) of Yang and Wang (2013).
                  phyt= -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet= -pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat= -pp_sig+two*(t1*p1_sig+t2*p2_sig) 
              else
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)  
              endif
          else
              If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).
                  phyt=pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=pp_sig+two*(t1*p1_sig+t2*p2_sig)
              endif
              If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).		
                  phyt=-pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  timet=-pp_t+two*(t1*p1_t+t2*p2_t)
                  sigmat=-pp_sig+two*(t1*p1_sig+t2*p2_sig)
              endif
          endif 
          IF(mu_tp1.eq.one)phyt = phyt+(t1+t2)*PI 
       !***************************************************** 
          else
              count_num=1
              goto 30
          endif
      ENDIF 
      RETURN
      END SUBROUTINE SPHERICALMOTION_BL
!********************************************************************************************
      SUBROUTINE SPINZERO(p,kp,kt,lambda,q,sinobs,muobs,a_spin,robs,&
                         thetamax,phyt,timet,sigmat,mucos,t1,t2) 
!******************************************************************************************** 
!*     PURPOSE:   This routine computs the four B_L coordiants r,\mu,\phi,t and affine 
!*                parameter \sigma of the spherical motion when black hole spin a is zero. 
!*     INPUTS:   p--------------the independent variable.
!*               kt-------------p_\theta, the \theta component of four momentum of the photon.
!*               kp-------------p_\phi,  the \phi component of four momentum of the photon.
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1). 
!*               robs-----------radial coordinate of observer or the initial position of photon. 
!*               thetamax-------the maximum or minimum value of the \theta coordinate of the geodesic.
!*                              which also the \theta turning point of the spherical motion.       
!*     OUTPUTS:  phyt-----------the \phi coordinat of the photon.
!*               timet----------the t coordinats of the photon.
!*               sigmat---------the affine parameter \sigma.
!*               mucos----------the \theta coordinates of the photon, and mucos=cos(\theta).  
!*               t1,t2----------number of times of photon meets turning points \mu_tp1 and \mu_tp2
!*                              respectively.      
!*     ROUTINES CALLED: metricg, circ_mb, weierstrass_int_J3, weierstrassP
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ****************************************** 
      use constants 
      IMPLICIT NONE
      DOUBLE PRECISION kp,kt,lambda,q,p,sinobs,muobs,a_spin,robs,phyt,timet,mucos,&
             kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,sigmat,AA,BB,&
             mu_tp1,mu_tp2,PI1,Ptotal,pp,p1,p2,PI1_phi,PI2_phi,PI1_time,PI2_time,PI1_sigma,&
             PI2_sigma,pp_sigma,p1_sigma,p2_sigma,pp_time,p1_time,p2_time,pp_phi,p1_phi,&
             p2_phi,p_mu,Delta,c_time,c_phi,PI2,mu,thetamax,thetamax_1
      integer :: t1,t2,count_num=1,i,j
      save :: kp_1,kt_1,lambda_1,q_1,sinobs_1,muobs_1,a_spin_1,robs_1,PI1_phi,PI2_phi,&
              PI1_time,PI2_time,PI1_sigma,PI2_sigma,Delta,c_time,c_phi,PI1,PI2,mobseqmtp,&
              AA,BB,mu_tp1,mu_tp2,Ptotal,thetamax_1
      logical :: mobseqmtp 

30    continue
      IF(count_num .eq. 1)THEN
          kp_1 = kp
          kt_1 = kt
          lambda_1 = lambda
          q_1 =q
          sinobs_1 = sinobs
          muobs_1 = muobs
          a_spin_1 = a_spin
          robs_1 = robs
          thetamax_1 = thetamax
          t1=0
          t2=0
          mobseqmtp = .false.
 
          IF(q.gt.zero)THEN
              AA = dsqrt((q+lambda*lambda)/q)
              BB = dsqrt(q)
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              If(kt.gt.zero)then
                  mu=sin(asin(muobs*AA)-p*BB*AA)/AA
              else
                  If(kt.eq.zero)then
                      mu=cos(p*AA*BB)*muobs
                  else
                      mu=sin(asin(muobs*AA)+p*AA*BB)/AA
                  endif
              endif
              mucos = mu  
          !****************************************************
              If(kt.ne.zero)then
                  mu_tp1=sqrt(q/(lambda**two+q))
                  mu_tp2=-mu_tp1
              else
                  mu_tp1=abs(muobs)
                  mu_tp2=-mu_tp1
                  mobseqmtp=.true.
              endif
              If(abs(muobs).eq.one)mobseqmtp=.true.
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  timet = zero
                  sigmat = zero
                  phyt = zero
                  return
              endif  
          !***************************************************
              PI1=(PI/two-asin(muobs/mu_tp1))*mu_tp1/BB
              Ptotal=PI*mu_tp1/BB
              PI2=Ptotal-PI1
              pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB
              p1=PI1-pp
              p2=Ptotal-p1 
              PI1_phi=zero
              PI2_phi=zero
              PI1_time=zero
              PI2_time=zero
              PI1_sigma=zero
              PI2_sigma=zero
              Do j=0,10**6
                  Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          else
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                      else
                          If(kt.lt.zero)then
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                          If(kt.gt.zero)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          endif
                      endif   
                      If(abs(p-p_mu).lt.1.D-6)goto 300
                  enddo
              enddo
300           continue
              !p0 part!
              Delta = robs*robs+a_spin*a_spin-two*robs
              c_phi = a_spin*(robs*(two)-(lambda*a_spin))/Delta
              c_time = ( (two)*robs**three+(two*a_spin*(a_spin-lambda))*robs )/Delta

              pp_sigma = a_spin*a_spin*mveone_int(muobs,mu,one/AA)/BB+robs*robs*pp
              pp_time = pp_sigma+c_time*pp
              pp_phi = lambda*schwatz_int(muobs,mu,AA)/BB+c_phi*pp
              !******p1 part***********************************
              IF(t1 .eq. 0)THEN
                  p1_sigma=zero
                  p1_time=zero
                  p1_phi=zero
              ELSE
                  IF(PI1_time .eq. zero .or. PI1_sigma .eq. zero)THEN
                      PI1_sigma = a_spin*a_spin*mveone_int(muobs,mu_tp1,one/AA)/BB+robs*robs*PI1
                      PI1_time = PI1_sigma+c_time*PI1 
                  ENDIF
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB+c_phi*PI1
                  ENDIF
                  p1_sigma = PI1_sigma-pp_sigma
                  p1_time = PI1_time-pp_time
                  p1_phi = PI1_phi-pp_phi
              ENDIF  
              !******p2 part***********************************
              IF(t2 .eq. 0)THEN
                  p2_sigma=zero
                  p2_time=zero
                  p2_phi=zero
              ELSE
                  IF(PI2_time .eq. zero .or. PI2_sigma .eq. zero)THEN
                      PI2_sigma = a_spin*a_spin*mveone_int(mu_tp2,muobs,one/AA)/BB+robs*robs*PI2
                      PI2_time = PI2_sigma+c_time*PI2 
                  ENDIF
                  IF(PI2_phi .eq. zero)THEN
                      PI2_phi = lambda*schwatz_int(mu_tp2,muobs,AA)/BB+c_phi*PI2
                  ENDIF
                  p2_sigma = PI2_sigma+pp_sigma
                  p2_time = PI2_time+pp_time
                  p2_phi = PI2_phi+pp_phi
              ENDIF
              !**********************************************  
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then  
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  else
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  endif  
              else
                 If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif
                 If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif
              endif
              If(thetamax.eq.zero.or.thetamax.eq.180.D0)phyt = phyt+(t1+t2)*PI
          ELSE
              !write(unit=6,fmt=*)'phyt_schwatz(): q<0, which is a affending',&
              !		'value, the program should be',&  
              !		'stoped! and q = ',q
              !stop
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyt = zero 
              timet = zero
              sigmat = zero
          ENDIF
      ELSE 
          IF(kp_1 .eq. kp.and.kt_1 .eq. kt.and.lambda_1 .eq. lambda.and.q_1 .eq.q.and.&
          sinobs_1 .eq. sinobs.and.muobs_1 .eq. muobs.and.a_spin_1 &
          .eq. a_spin.and.robs_1 .eq. robs .and.thetamax_1.eq.thetamax)THEN
      !*******************************************************************************
            IF(q.gt.zero)THEN 
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              If(kt.gt.zero)then
                  mu=sin(asin(muobs*AA)-p*BB*AA)/AA
              else
                  If(kt.eq.zero)then
                      mu=cos(p*AA*BB)*muobs
                  else
                      mu=sin(asin(muobs*AA)+p*AA*BB)/AA
                  endif
              endif
              mucos = mu   
          !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
              If(mu_tp1.eq.zero)then
              !photons are confined in the equatorial plane, 
              !so the integrations about !\theta are valished.
                  timet = zero
                  sigmat = zero
                  phyt = zero
                  return
              endif  
          !*************************************************** 
              pp=(asin(mu/mu_tp1)-asin(muobs/mu_tp1))*mu_tp1/BB
              p1=PI1-pp
              p2=Ptotal-p1  
              Do j=0,10**6
                  Do i=j,j+1
                      If(mobseqmtp)then
                          If(muobs.eq.mu_tp1)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          else
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                      else
                          If(kt.lt.zero)then
                              t1=i
                              t2=j
                              p_mu=pp+two*(t1*p1+t2*p2)
                          endif
                          If(kt.gt.zero)then
                              t1=j
                              t2=i
                              p_mu=-pp+two*(t1*p1+t2*p2)
                          endif
                      endif   
                      If(abs(p-p_mu).lt.1.D-6)goto 310
                  enddo
              enddo
310           continue
              !p0 part!   
              pp_sigma = a_spin*a_spin*mveone_int(muobs,mu,one/AA)/BB+robs*robs*pp
              pp_time = pp_sigma+c_time*pp
              pp_phi = lambda*schwatz_int(muobs,mu,AA)/BB+c_phi*pp
              !******p1 part***********************************
              IF(t1 .eq. 0)THEN
                  p1_sigma=zero
                  p1_time=zero
                  p1_phi=zero
              ELSE
                  IF(PI1_time .eq. zero .or. PI1_sigma .eq. zero)THEN
                      PI1_sigma = a_spin*a_spin*mveone_int(muobs,mu_tp1,one/AA)/BB+robs*robs*PI1
                      PI1_time = PI1_sigma+c_time*PI1 
                  ENDIF
                  IF(PI1_phi .eq. zero)THEN
                      PI1_phi = lambda*schwatz_int(muobs,mu_tp1,AA)/BB+c_phi*PI1
                  ENDIF
                  p1_sigma = PI1_sigma-pp_sigma
                  p1_time = PI1_time-pp_time
                  p1_phi = PI1_phi-pp_phi
              ENDIF  
              !******p2 part***********************************
              IF(t2 .eq. 0)THEN
                  p2_sigma=zero
                  p2_time=zero
                  p2_phi=zero
              ELSE
                  IF(PI2_time .eq. zero .or. PI2_sigma .eq. zero)THEN
                      PI2_sigma = a_spin*a_spin*mveone_int(mu_tp2,muobs,one/AA)/BB+robs*robs*PI2
                      PI2_time = PI2_sigma+c_time*PI2 
                  ENDIF
                  IF(PI2_phi .eq. zero)THEN
                      PI2_phi = lambda*schwatz_int(mu_tp2,muobs,AA)/BB+c_phi*PI2
                  ENDIF
                  p2_sigma = PI2_sigma+pp_sigma
                  p2_time = PI2_time+pp_time
                  p2_phi = PI2_phi+pp_phi
              ENDIF
              !**********************************************  
              If(mobseqmtp)then
                  If(muobs.eq.mu_tp1)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  else
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)
                  endif
              else
                 If(kt.lt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif
                 If(kt.gt.zero)then
!      equations (52) of Yang and Wang (2013).	
                      sigmat = -pp_sigma+two*(t1*p1_sigma+t2*p2_sigma)
                      timet = -pp_time+two*(t1*p1_time+t2*p2_time)
                      phyt = -pp_phi+two*(t1*p1_phi+t2*p2_phi)
                 endif
              endif
              If(thetamax.eq.zero.or.thetamax.eq.180.D0)phyt = phyt+(t1+t2)*PI
            ELSE 
              mucos=muobs 
              t1 = 0
              t2 = 0
              phyt = zero 
              timet = zero
              sigmat = zero
            ENDIF
          ELSE
              count_num = 1
              goto 30
          ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SPINZERO

!********************************************************************************************
      FUNCTION mveone_int(x,y,z0)
!******************************************************************************************** 
!*     PURPOSE:   To compute the integral \int^y_x z^2/sqrt(1-(z/z0)^2)dz, and z0<1. 
!*     INPUTS:    x,y-----------the integral limits.    
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************  
      USE constants
      implicit none
      double precision x,y,z0,mveone_int,xt,yt,Iy,Ix

      xt = x
      yt = y
      If(xt .eq. yt)THEN
          mveone_int = zero
          return
      ENDIF 
      Iy = z0**three*half*(dasin(yt/z0)-y/z0/z0*dsqrt(z0*z0-y*y))
      Ix = z0**three*half*(dasin(xt/z0)-x/z0/z0*dsqrt(z0*z0-x*x))
      mveone_int = Iy-Ix  
      return
      END FUNCTION mveone_int

!*****************************************************************************************************
      Subroutine lambdaq_sphericalm(r_sm,a_spin,lambda,q,theta_min)
!******************************************************************************************** 
!*     PURPOSE:   To compute the constants of motion for the spherical motion of a photon. 
!*     INPUTS:    r_sm--------------the radius of the spherical motion.
!*                theta_max---------the maximum or minimum value the \theta coordinate of the 
!*                                  spherical motion, which also the turning point of the motion
!*                                  in theta coordinate.  
!*                signs-------------which determine the direction of the motion of a photon with respect 
!*                                  to the black hole. signs>0 is for co-rotating/prograde orbiting,
!*                                  signs< is for counter-rotating/retrograde orbiting.
!*                a_spin------------the black hole spin.
!*     OUTPUTS:   lamda,q-----------constants of motion.    
!*                theta_min---------the minimum value of the theta coordinate the orbiting.                         
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  5 Jan 2012
!*     REVISIONS: ******************************************  
      implicit none
      Double precision :: r_sm,a_spin,zero,one,two,three,four,&
             sinthemax,costhemax,cotmax,dtor,AL,AE,Sigma,Delta,lambda,q,&
             a1,b1,c1,mu_starm,r_max,r_min,c_temp,mutp2,sintp2,theta_min
      parameter(dtor=asin(1.D0)/90.D0,zero=0.D0,one=1.D0,two=2.D0,three=3.D0,four=4.D0)
	
      r_max = radiusofsphericalmotion(a_spin,zero)
      c_temp = 90.D0
      r_min = radiusofsphericalmotion(a_spin,c_temp)
      If(a_spin.lt.zero)then
          c_temp = r_max
          r_max = r_min
          r_min = c_temp
      endif
      write(*,*)'r_min=',r_min,'   r_max=',r_max
      
      If (a_spin.ne.zero) then   
          If (r_sm.lt.r_min .or. r_sm.gt.r_max .or. r_sm-r_min.le.-1.D-6 .or. r_sm-r_max.ge.1.D-6) then
              write(*,*)'lambdaq_sphericalm(): For spin=',a_spin,' the allowed range for radial coordinate '&
              ,'of the spherical motion of a photon is between ',r_min,' and ',r_max,'. The radius you input is '&
              ,'out of the range and the code shall stop.'
              stop
          endif
      else
          If(r_sm.ne.three)then
              write(*,*)'For a=0, the radius of the spherical motion of a photon is 3 for all inclination anges of '&
                       ,'the orbit with respect to the equatorial plane.'
              r_sm = three
          endif
      endif
!*******************************************************************
      If (a_spin.ne.zero) then
          a1 = a_spin**four*(r_sm-one)**two
          b1 = two*a_spin**two*r_sm*(two*a_spin*a_spin-three*r_sm+r_sm**three)
          c1 = r_sm**three*(-four*a_spin*a_spin+r_sm*(r_sm-three)**two)
    
          mutp2 = (-b1+dsqrt(b1*b1-four*a1*c1))/two/a1 
          sintp2 = one-mutp2
          mu_starm = (-b1-dsqrt(b1*b1-four*a1*c1))/two/a1  
          theta_min = dacos(dsqrt(mutp2))/dtor
          write(*,*)'Theta_min is:',dacos(dsqrt(mutp2))/dtor 
      else
          write(*,*)'For a=0, we need you to input the minimum of the theta coordinate:'
          read(unit=5,fmt=*)theta_min      
          If(theta_min.ne.90.D0)then
              sintp2=sin(theta_min*dtor)**two
              mutp2=cos(theta_min*dtor)**two
          else
              sintp2=one
              mutp2=zero
          endif
      endif
!*******************************************************************
      cotmax=costhemax/sinthemax
      Sigma=r_sm**two+(a_spin)**two*mutp2
      Delta=r_sm**two-two*r_sm+a_spin**two

      AL=Delta*(three*r_sm**four+(a_spin*r_sm)**two*&
                    (mutp2+one)-a_spin**four*mutp2)
      AE=Delta*(r_sm**two-(a_spin)**two*mutp2)

      lambda = dsign(one,a_spin)*dsqrt(sintp2*AL/AE)
      q = mutp2*(AL/AE-a_spin*a_spin) 
      return
      End Subroutine lambdaq_sphericalm

!********************************************************************************
      function Denomenator_r(r,theta,a_spin)
!********************************************************************************
      use constants
      implicit none
      double precision r,theta,a_spin,Denomenator_r,sinthe,costhe

      If(theta.ne.90.D0)then
          sinthe=dsin(theta*dtor)
          costhe=dcos(theta*dtor)
      else
          sinthe=one
          costhe=zero
      endif
      Denomenator_r = ( -(  -three*r*r+(r+one)*a_spin*a_spin*costhe**two+two*a_spin*sinthe*& 
               dsqrt( r*(r*r-a_spin*a_spin*costhe**two) )   ) )**(one/three) 
      return
!********************************************************************************
      end function

!********************************************************************************
      function radiusofsphericalmotion(a,c)
!********************************************************************************
      implicit none
      double precision a,c,radiusofsphericalmotion,r1,r2,Dr

      r1 = 30.D0
      r2 = Denomenator_r(r1,c,a)
      Dr = r2-r1
      r1 = r2
      do while(dabs(Dr).ge.1.D-10)
          r2 = Denomenator_r(r1,c,a)
          Dr = r2-r1
          r1 = r2
      enddo
      radiusofsphericalmotion = r1
      return
      end function
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      end module sphericalmotion

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



!********************************************************************************************
      module pemfinding
!*******************************************************************************
!*     PURPOSE: This module aims on solving more general equation f(p)=0. For 
!*              detail definition of f(p), cf. Yang & Wang (2012).     
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  15 Jan 2012. 
!*******************************************************************************        
      use constants
      use rootsfinding
      use ellfunction
      USE BLcoordinate                
      implicit none

      contains
!******************************************************************************
      function pemfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&                              
                rin,rout,muup,mudown,phy1,phy2,caserange,Fp,paras,bisection)   
!******************************************************************************
!*     PURPOSE:  Searches for minimum root pem of equations f(p)=0.     
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_\phi, p_0. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               rin,rout-------Inner and outer radius of emission region or emission objects.
!*               muup,mudown----Boundary \mu coordinates of emission region or objects.
!*                              Or the maximum and minimum of \mu coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               phy1,phy2------Boundary \phi coordinates of emission region or objects.
!*                              Or the maximum and minimum of \phi coordinates of emission region or objects.
!*                              These parameter can be specified void value by setting
!*                              values of parameter caserange. 
!*               caserange------Tell routine whether muup, mudown and phy1, phy2 are provided.
!*                              caserange = 1, muup, mudown and phy1, phy2 are provided.
!*                              caserange = 2, muup, mudown are provided, but phy1, phy2 not.
!*                              caserange = 3, muup, mudown are not provided, phy1, phy2 are provided.
!*                              caserange = 1, muup, mudown and phy1, phy2 not are provided. 
!*                              Provided means corresponding parameters have been set specific value.
!*                              Not provided means corresponding parameters have not been set specific value,
!*                              but these parameter should also be given as dummy parameter.
!*               Fp-------------name of function f(p). This routine to compute f(p) should be prvided by
!*                              user, and the dummy variable of Fp should have following form:
!*                              Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras).
!*               paras(1:10)----Array of parameters to descirbe function f(p).  
!*               bisection------Logical variable, if TURE, then use bisection method to search the 
!*                              root of equation f(p)=0, else use Newton-Raphson method. 
!*               NN-------------In the subroutine pemfind there is an important parameter: NN, which is the number
!*                              of sections of the interval (p1 , p2 ) or (p3 , p4 ) has been divided when 
!*                              searching the roots. One needs to set a proper value for NN, since if 
!*                              NN is too small, the roots exit on the interval (p1 , p2 ) or (p3 , p4 ) 
!*                              maybe omitted, if NN is too big, the time consumed by the code will be
!*                              large.
!*     OUTPUTS:  pemfind--------value of root of equation f(p)=0 for p.  
!*                              pemdisk=-1.D0, if the photon goto infinity.
!*                              pemdisk=-2.D0, if the photon fall into event horizon.         
!*     REMARKS:  This routine will search root between interval (p1, p2). We will chose NN points on 
!*               this interval, and check one by one to wheter f(p) changing its sign, if so, the root
!*               must be on interval (p_{i-1}, p_{i}). Then we use Bisection or Newton-Raphson method 
!*               to find the roots. One should set NN propriately to guarantee no root missing and
!*               also quickly find the root, thus saving the running time of the code.
!*     ROUTINES CALLED: radiustp, r2p, rootfind, Sectionp.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision pemfind,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,rhorizon,NN,r_tp1,r_tp2,deltax,&
                        p1,p2,paras(10),f1234(4)
        Double precision rin,rout,muup,mudown,phy1,phy2
        optional rin,rout,muup,mudown,phy1,phy2
        complex*16  bb(1:4)
        integer  reals,cases,cases_of_tp,caserange,tr1,tr2
        Double precision,external :: Fp!,Bisection,rootfind,Sectionp
        parameter(deltax=5.D-5)
        logical :: robs_eq_rtp,indrhorizon,bisection

        rhorizon=one+sqrt(one-a_spin**two)        
        call radiustp(f1234(1),a_spin,robs,lambda,q,r_tp1,r_tp2,&
                reals,robs_eq_rtp,indrhorizon,cases_of_tp,bb)
!NN is important paramerter in the searching of roods by Bisection or Newton-Raphson algorithm.
!Which divides the internal [p1,p2] or [p3,p4] into NN parts. If NN is too small then one may
!miss the roots, while if NN is too big, the speed of code will very slow. Thus one should to
!chose a apropriate one. 
        IF(r_tp1.ge. 10.D0)THEN
            NN=30.D0
        ELSE
            NN=30.D0
        ENDIF 
! To classify the geodesic according its relationship with a shell.
! case 1, 2, 3, 4 correspond to A, B, C, D discussed in Yang & Wang (2012).
        If(present(rout).and.present(rin))then
            If(rin.gt.rhorizon)then
                If(r_tp1.ge.rout)then 
                    cases=1 
                    pemfind=-one 
                else
                    If(r_tp1.gt.rin)then
                        cases=2
                    else
                        If(r_tp1.gt.rhorizon)then
                            cases=3
                        else
                            cases=4        
                        endif
                    endif        
                endif
            else
                cases=5 
            endif
        else
            cases=5
        endif
        !write(*,*)'cases=',cases,r_tp1,rout,rin,f1234(2)
        select case(cases)
        case(1)               
        case(2)
            tr1=0
            tr2=0
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)
            tr1=1        
            p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)  
!write(*,*)'p1,p2=',p1,p2 
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,bisection)
            endif 
        case(3)
            tr1=0
            tr2=0        
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            p2=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2)        
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,bisection)
                !write(*,*)'ff=',pemfind
            endif
            If(pemfind.eq.-one)then
                tr1=1
                p1=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2)        
                p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)  
                p2 = p2!+1.D-2	     
                If(caserange.eq.4)then        
                    pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,bisection)        
                else        
                    pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                                tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,bisection)
                    !write(*,*)'ff=',p1,p2,pemfind
                endif                
            endif                        
        case(4)        
            tr1=0
            tr2=0
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            p2=r2p(f1234(1),rin,lambda,q,a_spin,robs,scal,tr1,tr2)     
            If(caserange.eq.4)then        
                pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,bisection)        
            else        
                pemfind=Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        tr1,tr2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,bisection)
            endif
!the photon will fall into the black hole.
            If(pemfind.eq.-one)then 
                pemfind=-two
            endif        
        case(5)
            tr1=0        
            tr2=0           
            p1=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)        
            If(r_tp1.gt.rhorizon)then
                tr1=1
                p2=r2p(f1234(1),rout,lambda,q,a_spin,robs,scal,tr1,tr2)
            else
                p2=r2p(f1234(1),rhorizon,lambda,q,a_spin,robs,scal,tr1,tr2)
            endif
            pemfind=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,tr1,tr2,p1,p2,NN,Fp,paras,bisection)
        end select
        return        
      end function pemfind
!*****************************************************************************************************
      Function Sectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                t1,t2,p1,p2,muup,mudown,phy1,phy2,caserange,NN,Fp,paras,bisection)  
!*****************************************************************************************************
!*     PURPOSE:  To judge whether a geodesic intersects with the surface of object or emission region
!*               which are described by function f(p). The space range of surface of object or
!*               emission region is (rin, rout) in radius, (mudown,muup) in poloidal and 
!*               (phy1, phy2) in azimuthal. If geodesic has no intersections on those intervals
!*               then it will no intersects with the surface and emission region, a special       
!*               value -1.D0 will be returned.
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2).
!*               muup,mudown,phy1,phy2,
!*               caserange,Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  Sectionp-------value of root of equation f(p)=0 for p.  
!*                              If Sectionp=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision Sectionp,lambda,q,sinobs,muobs,a_spin,robs,&
                        scal,NN,mu1,mu2,paras(10),deltax,&
                        phya1,phya2,p1,p2,f1234(4)
        Double precision ,optional :: muup,mudown,phy1,phy2
        integer  caserange,t1,t2
        Double precision,external :: Fp!,Bisection,rootfind
        parameter(deltax=5.D-5)
        logical :: bisection

            p1=p1-deltax
            p2=p2+deltax
            If(caserange.eq.1 .or.caserange.eq.3)then        
                phya1=phi(p1,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
                phya2=phi(p2,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)
            endif                
            If(caserange.eq.1 .or.caserange.eq.2)then        
                mu1=mucos(p1,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
                mu2=mucos(p2,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal) 
                If((mu1-muup)*(mu1-mudown).lt.zero.or.(mu2-muup)*(mu2-mudown).lt.zero.or.&
                (muup-mu1)*(muup-mu2).lt.zero.or.(mudown-mu1)*(mudown-mu2).lt.zero)then
                    If(caserange.eq.1 .or.caserange.eq.3)then 
                        If((phya1-phy1)*(phya1-phy2).lt.zero.or.(phya2-phy1)*(phya2-phy2).lt.zero.or.&
                        (phy1-phya1)*(phy1-phya2).lt.zero.or.(phy2-phya1)*(phy2-phya2).lt.zero)then
! the geodesic intersecte the zone defined by r1,r2,mu1,mu2,phy1,phy2,so it has the 
!possibility to hit the surface of the object.        
                            Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                        robs,scal,t1,t2,p1,p2,NN,Fp,paras,bisection) !(1,1)
!write(*,*)'ff=',p1,p2,Sectionp
                        else 
!the (phy1,phy2) and (phya1,phya2) has no public point,so we don't consider it.
                            Sectionp=-one  !(1,1)                
                        endif
                    else
                        Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,t1,t2,p1,p2,NN,Fp,paras,bisection) !(1,0)
                    endif
                else 
!the internal of (mu1,mu2) and (muup,mudown) does not overfold each other,so the geodesic will not hit the
!surface of the object at the internal (p_rout,p_rout2),which also means the geodesic will never hit the object.
                    Sectionp=-one          ! so nothing needed to be done.                        
                endif
            else
                If(caserange.eq.1 .or.caserange.eq.3)then
                    If((phya1-phy1)*(phya1-phy2).lt.zero.or.(phya2-phy1)*(phya2-phy2).lt.zero.or.&
                    (phy1-phya1)*(phy1-phya2).lt.zero.or.(phy2-phya1)*(phy2-phya2).lt.zero)then
                        Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,NN,Fp,paras,bisection) !(0,1)
                    else 
!the (phy1,phy2) and (phya1,phya2) has no public points,so we don't consider it.
                        Sectionp=-one  !(0,1)                
                    endif
                else              
                    Sectionp=rootfind(f1234,lambda,q,sinobs,muobs,a_spin,&
                                robs,scal,t1,t2,p1,p2,NN,Fp,paras,bisection)        ! (0,0)        
                endif        
            endif
        return
      End Function Sectionp


!********************************************************************************* 
      Function rootfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                           t1,t2,p1,p2,NN,Fp,paras,bisection)
!********************************************************************************* 
!*     PURPOSE:  To search roots on interval (p1, p2). If no roots were found 
!*               a special value -1.D0 will return. 
!* 
!*     INPUTS:   f1234(1:4)-----array of p_r, p_theta, p_phi, p_t which are the four momentum of the photon. 
!*               lambda,q-------motion constants, defined by lambda=L_z/E, q=Q/E^2. 
!*               sinobs,muobs---sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
!*                              \theta_{obs} is the inclination angle of the observer.
!*               a_spin---------spin of black hole, on interval (-1,1).  
!*               robs-----------radial coordinate of observer or initial position of photon. 
!*               scal-----------a dimentionless parameter to control the size of the images.
!*                              Which is usually be set to 1.D0. 
!*               t1,t2----------Number of photon meets turning points r_tp1 and r_tp2 respectively.
!*               p1,p2----------roots may hidden inside interval (p1,p2). 
!*               Fp,paras,bisection--------see instructions in pemfind.  
!*
!*     OUTPUTS:  rootfind-------value of root of equation f(p)=0 for p.  
!*                              If rootfind=-1.D0, No roots were found no interval (p1, p2).     
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision rootfind,lambda,q,sinobs,muobs,a_spin,&
                        robs,scal,p1,p2,NN,sp1,sp2,&
                        f_p,p,paras(10),f1234(4),deltap,dp
        Double precision ,external :: Fp
        parameter (dp=1.D-5)
        integer NNf,k,t1,t2
        logical :: bisection
 
        !p1 = p1 + dp
        deltap=(p2-p1)/NN
        NNf=floor(NN)
        p=p1
        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
        !write(*,*)'NNf=',f_p
        !write(unit=6,fmt=*)p1,p2,f_p
        If(f_p.eq.0.D0)then
                rootfind=p
                return
        endif 
        If(f_p.lt.zero)then
                k=0
                Do while(.true.)
                        If(f_p.eq.zero)then
                            rootfind=p
                            return
                        endif 
                        If(f_p.gt.zero .or.k.gt.NNf)exit 
                        k=k+1
                        p=p1+deltap*k
                        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
                Enddo                                        
        else
                k=0        
                Do while(.true.)
                        If(f_p.eq.zero)then
                            rootfind=p
                            return
                        endif
                        If(f_p.lt.zero.or.k.gt.NNf)exit 
                        k=k+1
                        p=deltap*k+p1
                        f_p=Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)        
                Enddo
        endif
 !write(unit=6,fmt=*)'f_p=',f_p,k,deltap,p-deltap,p
                If(k.le.NNf)then
                    sp1=p-deltap
                    sp2=p      
! Using bisection or Newton Raphson method to find roots on interval (sp1, sp2).   
                    IF(bisection)THEN
                        rootfind=Bisectionp(f1234,lambda,q,sinobs,muobs,&
                                a_spin,robs,scal,t1,t2,sp1,sp2,Fp,paras)
                    else
                        rootfind=NewRapson(f1234,lambda,q,sinobs,muobs,&
                                a_spin,robs,scal,t1,t2,sp1,sp2,Fp,paras) 
                    endif
                else
!In (p1,p2) no roots were found!
                        rootfind=-1.D0        
                endif
        return
      End Function rootfind
!*****************************************************************************************************
      Function Bisectionp(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,Fp,paras)
!*****************************************************************************************************
!*     PURPOSE:  To search roots on interval (p1, p2) by bisenction method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  Bisectionp-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*     REVISIONS: ****************************************** 
        use BLcoordinate
        implicit none
        Double precision Bisectionp,lambda,q,p1,p2,sinobs,muobs,&
                        a_spin,robs,scal,pc,f1,&
                        f2,fc,counter,paras(10),f1234(4)
        Double precision,external :: Fp
        integer t1,t2
        
        pc=(p1+p2)/two
        f1=Fp(p1,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
        f2=Fp(p2,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
        fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)        
!write(unit=6,fmt=*)'fc=',f1,f2,fc
        counter=0
        Do while(abs(p2-p1).gt.1.D-4)
             If(f1*fc.gt.zero)then
                p1=pc
                f1=fc
                pc=(p1+p2)/two 
                fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
             Endif        
             If(f2*fc.gt.zero)then
                p2=pc
                f2=fc
                pc=(p1+p2)/two
                fc=Fp(pc,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
             Endif
             IF(fc.eq.zero)Then 
                 exit
             Endif
             If(counter.gt.200)exit
!write(unit=6,fmt=*)'fc=',f1,fc,f2,p1,pc,p2
             counter=counter+1                
        Enddo
        Bisectionp=pc        
        return
      End Function Bisectionp 
!*****************************************************************************************************
      Function NewRapson(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,p1,p2,Fp,paras)
!*****************************************************************************************************
!*     PURPOSE:  To search roots on interval (p1, p2) by Newton Raphson method for 
!*               equation Fp(p)=0.
!*     INPUTS:   Parameters descriptions ---------see instructions in pemfind.  
!*
!*     OUTPUTS:  NewRapson-------value of root of equation f(p)=0 for p.        
!*     ROUTINES CALLED: phi, mucos, rootfind. 
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012
!*****************************************************************
        use BLcoordinate
        implicit none
        Double precision NewRapson,lambda,q,sinobs,muobs,a_spin,robs,&
               scal,p1,p2,deltap,f_p,f_pj,paras(10),f1234(4),dfp,&
               dp,pacc,EPS,temp,h,pj,ptn
        Double precision ,external :: Fp
        integer k,kMax,t1,t2
        parameter (kmax=30, deltap=1.D-6, pacc=1.D-2, EPS=1.D-5)

        ptn=(p2+p1)/two
        Do k=1,kmax
            temp=ptn
            h=EPS*abs(temp)
            if(h.eq.zero)h=EPS        
            pj=temp+h
            h=pj-temp
            f_p=Fp(temp,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
            f_pj=Fp(pj,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
            dfp=(f_pj-f_p)/h        
            dp=f_p/dfp
            ptn=ptn-dp
!If((ptn-p1)*(p2-ptn).lt.zero)then
!write(unit=6,fmt=*)'ptn jumps out of brakets [p1, p2]!'
            If(abs(dp).lt.pacc)then
                NewRapson=ptn
                return
            endif  
        Enddo
      End Function NewRapson
!******************************************************* 
      end module pemfinding 



!*******************************************************************************
      module obsemitter
!*******************************************************************************
!*     PURPOSE: This module aims on determining geodesic connecting observer
!*              and emitter, i.e. to solve equations (121)-(123) in  
!*              Yang & Wang (2012) by Newton Raphson method.     
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012).  
!*     DATE WRITTEN:  15 Jan 2012. 
!*******************************************************************************  
      use constants
      use rootsfinding
      use ellfunction
      use pemfinding
      use BLcoordinate         
      implicit none
!*****************************************************************************************************
      contains
!*****************************************************************************************************
      SUBROUTINE alphabetap(a0,B0,alen,Blen,sinobs,muobs,a_spin,robs,scal,&
                              obs_V,rmuphy_em,abp,func1,func2,func3)
      IMPLICIT NONE
!*******************************************************************************
!*    PURPOSE:   This subroutine aim to determine motion constants of a geodesic which connecting the 
!*               emmiter and the observer. The Boyer-Lindquist coordinates of them (the emmiter and 
!*               observer) have given.
!************************************************************************
!*     INPUTS:    (a0,B0)--------is the coordinates of one corner of the rectangle on the 
!*                               screen plane of the observer.
!*                alen,Blen------size of the rectangle. 
!*                sinobs,muobs---sinobs=sin(theta_obs) and muobs=cos(theta_obs), theta_obs is the 
!*                               theta coordinate of the observer.
!*                a_spin---------spin of black hole.                       
!*                robs-----------The radiual coordinate of the observer.
!*                scal-----------a parameter to control the size of the image.   
!*                               The value of which usually was set to be 1.D0.
!*                obs_V(1:3)-----array of the velocity of the observer respect to the ZAMO( or LNRF).    
!*                rmuphy_em------array of Boyer-Lindquist coordinates of the emmiter. 
!*                func1,func2,func3-----Name of functions: r(p), \mu(p), \phi(p).
!*     OUTPUTS:   abp(1:3)-------array of roots of equations (121)-(123) in Yang & Wang (2012),  
!*                               \alpha_em, \beta_em, p_em.
!*     ROUTINES CALLED: lambdaq, mnewt2.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012       
!*******************************************************************************  
      INTEGER n,ntrial,i,j
      INTEGER tr1,tr2
      PARAMETER(n=3,ntrial=1000)
      DOUBLE PRECISION a0,B0,a_spin,robs,scal,x(n),xem(n),tolx,tolf,sinobs,muobs,a1,B1,&
                   zero,one,two,three,four,rmss,&
                   p,lambda,q,abp(3),rmuphy_em(3),&
                   f1234(4),obs_V(3),alpha1,Beta1,&
                   alen,Blen,adel,Bdel        
      PARAMETER(two=2.D0,three=3.D0,four=4.D0,one=1.D0,zero=0.D0)
      DOUBLE PRECISION,EXTERNAL :: func1,func2,func3
      LOGICAL :: err

      p=one
      tr1=0
      tr2=0
      rmss=rms(a_spin) 
      tolx=1.D-4
      tolf=1.D-4
      xem(1)=rmuphy_em(1)
      xem(2)=rmuphy_em(2)
      rmuphy_em(3)=DMOD(rmuphy_em(3),twopi)
      If(rmuphy_em(3).lt.zero)rmuphy_em(3)=rmuphy_em(3)+twopi
      xem(3)=rmuphy_em(3)
      adel=0.1D0
      Bdel=0.1D0         
        B1=-sqrt(xem(1)**two+a_spin**two)*cos(xem(3))*zero
        a1=-sqrt(xem(1)**two+a_spin**two)*sin(xem(3))*zero
        Do i=0,floor(Blen/Bdel)
            Beta1=dabs(B0)-i*Bdel
            Do j=0,floor(alen/adel)
            alpha1=-dabs(a0)+j*adel
                x(1)=alpha1
                x(2)=Beta1
                tr1=0
                tr2=0
                call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
                x(3)=r2p(f1234(1),xem(1),lambda,q,a_spin,robs,scal,tr1,tr2)
                If(x(3).lt.zero)cycle        
                call mnewt2(ntrial,x,n,tolx,tolf,sinobs,muobs,a_spin,robs,scal,&
                                                 obs_V,xem,err,func1,func2,func3) 
                !write(*,*)'i,j',i,j,alpha1,Beta1,lambda,q,x,err
                If(.not.err)then
                !write(*,*)'ggg',x,err
                        abp(1)=x(1)
                        abp(2)=x(2)
                        abp(3)=x(3)
                        goto 100
                endif           
            Enddo   
        Enddo   
100     continue
        return
      END subroutine alphabetap
!********************************************************************************* 
      SUBROUTINE mnewt2(ntrial,x,n,tolx,tolf,sinobs,muobs,a_spin,robs,&
                                 scal,obs_V,xend,err,func1,func2,func3) 
!************************************************************************************** 
!*    PURPOSE:    Given an initial guess x for a root in n dimensions,take ntrial 
!*                Newton-Raphson steps to improve the root.Stop if the root converges 
!*                in either summed absolute variable imcrements tolx
!*                or summed absolute function values tolf.
!*         
!*     INPUTS:    
!*     ROUTINES CALLED:  lubksb,ludcmp,usrfun.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:   
!*******************************************************************************  
        IMPLICIT NONE
        INTEGER n,ntrial,NP
        
        DOUBLE PRECISION tolx,tolf,x(n),zero,xend(n)
        PARAMETER(NP=3,zero=0.D0)     !Up to NP variables.
        INTEGER i,k,indx(NP)
        DOUBLE PRECISION d,errf,errx,fjac(NP,NP),fvec(NP),p(NP),&
                         sinobs,muobs,a_spin,robs,scal,obs_V(3)
        DOUBLE PRECISION, EXTERNAL :: func1,func2,func3
        LOGICAL err

        do k=1,ntrial 
            call usrfun2(x,n,NP,fvec,fjac,sinobs,muobs,a_spin,robs,scal,&
                                          obs_V,xend,err,func1,func2,func3) 
            If(err)then
                return
            endif
!User subroutine supplies function valuse at x in fvec and
!Jacobian matrix in fjac. Check function convergence 
            errf=zero                           
            do i=1,n                            
                errf=errf+abs(fvec(i))
            enddo 
            if(errf.le.tolf)return
            do i=1,n
                p(i)=-fvec(i)
            enddo        
!Solve linear equations using LU decomposition.
            call ludcmp(fjac,n,NP,indx,d,err)  
                If(err)then
!write(unit=6,fmt=*)'ludcmp(): err=',err,x(3)
                    return        
                endif
            call lubksb(fjac,n,NP,indx,p) 
!Check root convergence.
            errx=zero                                 
            do i=1,n
                errx=errx+abs(p(i))
                x(i)=x(i)+p(i)
            enddo
            If(x(3).lt.zero)then
                err=.true.
                return
            endif        
!write(*,*)'ntrials=',x
            if(errx.le.tolx)return
        enddo
        return
      END SUBROUTINE mnewt2
!********************************************************************* 
      SUBROUTINE usrfun2(x,n,NP,fvec,fjac,sinobs,muobs,a_spin,robs,&
                                   scal,obs_V,xend,err,func1,func2,func3)
!*********************************************************************
!*    PURPOSE:   This subroutine aim to compute array: fvec and fjac.
!************************************************************************
!*     INPUTS:     
!*     OUTPUTS:  fvec and fjac 
!*     ROUTINES CALLED: lambdaq, mnewt2.
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Yang & Wang (2012)  
!*     DATE WRITTEN:  6 Jan 2012       
!*******************************************************************************  
      USE BLcoordinate 
      IMPLICIT NONE                
      INTEGER  n,NP,i,j
      DOUBLE PRECISION x(n),fjac(NP,NP),fvec(NP),sinobs,muobs,a_spin,robs,scal,&
                   deltax,f(NP),h,EPS,&
                   temp,xend(n),rhorizon,lambda,q,f1234(4),obs_V(3)
      PARAMETER(deltax=1.D-6,EPS=1.D-5)
      DOUBLE PRECISION,EXTERNAL :: func1,func2,func3
      LOGICAL err
        
      err=.false.
      rhorizon=one+sqrt(one-a_spin**two)
      call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
      fvec(1)=func1(x(3),f1234(1),lambda,q,a_spin,robs,scal)-xend(1)
      fvec(2)=func2(x(3),f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)-xend(2)
      fvec(3)=func3(x(3),f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)-xend(3)

      If(fvec(1)+xend(1).ge.robs.or.fvec(1)+xend(1).le.rhorizon)then
          err=.true.
          return
      Endif
                
      Do j=1,n
          temp=x(j)
          h=EPS*ABS(temp)
          IF(h.eq.zero)h=EPS        
          x(j)=temp+h
          call lambdaq(x(1),x(2),robs,sinobs,muobs,a_spin,scal,obs_V,f1234,lambda,q)
          h=x(j)-temp
          f(1)=func1(x(3),f1234(1),lambda,q,a_spin,robs,scal)-xend(1)
          f(2)=func2(x(3),f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)-xend(2)
          f(3)=func3(x(3),f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)-xend(3)
          x(j)=temp
          Do i=1,n
              fjac(i,j)=(f(i)-fvec(i))/h
          Enddo
      Enddo        

      END SUBROUTINE usrfun2
!********************************************************************* 
      SUBROUTINE ludcmp(a,n,np,indx,d,err) 
!***********************************************************************************************
!*    PURPOSE:     PGiven a matrix a(1:n,1:n), with physical dimension np by np, this routine 
!*                 replaces it by the LU decomposition of a rowwise permutation of itself. a 
!*                 and n are input. a is output, arranged as in equation (2.3.14) above; 
!*                 indx(1:n) is an output vector that records the row permutation effected 
!*                 by the partial pivoting; d is output as +/-1 depending on whether
!*                 the number of row interchanges was even or odd, respectively. This 
!*                 routine is used in combination with lubksb to solve linear equations 
!*                 or invert a matrix. 
!*
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:  
!***********************************************************************************************
      IMPLICIT NONE
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY 
      PARAMETER (NMAX=500,TINY=1.0D-20)     !Largest expected n, and a small number.
      INTEGER i,imax,j,k
      LOGICAL  err
      DOUBLE PRECISION aamax,dum,sums,vv(NMAX)   !vv stores the implicit scaling of each row.
      d=one                            !No row interchanges yet.
      err=.false. 
      do i=1,n                         !Loop over rows to get the implicit scaling information.
          aamax=zero
          do j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
              !write(*,*)'aa==',aamax
          enddo
          if (aamax.eq.zero)then
              !write(unit=6,fmt=*)'pause singular matrix in ludcmp'
              !write(*,*)'ss==',a(i,1),a(i,2),a(i,3)
              err=.true.
              return
              !stop
          endif         !pause 'singular matrix in ludcmp'   !No nonzero largest element.
          vv(i)=one/aamax                                   !Save the scaling.
      enddo
      do  j=1,n                       !This is the loop over columns of Crout's method.
          do  i=1,j-1                 !This is equation (2.3.12) except for i = j.
              sums=a(i,j)
              do  k=1,i-1
                  sums=sums-a(i,k)*a(k,j)
              enddo 
              a(i,j)=sums
          enddo 
          aamax=zero      !Initialize for the search for largest pivot element.
          do  i=j,n            !This is i = j of equation (2.3.12) and i = j+1: ::N
              sums=a(i,j)                 !of equation (2.3.13).
              do  k=1,j-1
                  sums=sums-a(i,k)*a(k,j)
              enddo 
              a(i,j)=sums
              dum=vv(i)*abs(sums)                 !Figure of merit for the pivot.
              if (dum.ge.aamax) then         !Is it better than the best so far?
                  imax=i
                  aamax=dum
              endif
          enddo 
          if (j.ne.imax)then                 !Do we need to interchange rows?
              do  k=1,n                         !Yes, do so...
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
              enddo 
              d=-d                         !...and change the parity of d.
              vv(imax)=vv(j)                 !Also interchange the scale factor.
          endif
          indx(j)=imax
          if(a(j,j).eq.zero)a(j,j)=TINY
!*  If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
!*  For some applications on singular matrices, it is desirable to substitute TINY
!*  for zero.
          if(j.ne.n)then
!Now, finally, divide by the pivot element.                                
              dum=one/a(j,j)
              do  i=j+1,n
                  a(i,j)=a(i,j)*dum
              enddo 
          endif
!Go back for the next column in the reduction.
      enddo   
      return
      END  SUBROUTINE ludcmp
!**************************************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b) 
!**************************************************************************************
!*    PURPOSE:   Solves the set of n linear equations A \dot X = B. Here a is input, 
!*               not as the matrix A but rather as its LU decomposition, determined 
!*               by the routine ludcmp. indx is input as the permutation vector returned 
!*               by ludcmp. b(1:n) is input as the right-hand side vector B, and returns 
!*               with the solution vector X. a, n, np, and indx are not modified by this 
!*               routine and can be left in place for successive calls with diffierent 
!*               right-hand sides b. This routine takes into account the possibility that 
!*               b will begin with many zero elements, so it is efficient for use in 
!*               matrix inversion.
!*
!*     ACCURACY:   Machine.    
!*     AUTHOR:     Press et al. (2007)  
!*     DATE WRITTEN:  
!**************************************************************************************
      IMPLICIT NONE
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n),zero,sums
      INTEGER  i,ii,j,ll
      ii=0      !When ii is set to a positive value, it will become the index
                !of the first nonvanishing element of b. We now do
                !the forward substitution, equation (2.3.6). The only new
                !wrinkle is to unscramble the permutation as we go.
      do  i=1,n
          ll=indx(i)
          sums=b(ll)
          b(ll)=b(i)
          if (ii.ne.0)then
              do  j=ii,i-1
                  sums=sums-a(i,j)*b(j)
              enddo 
          else 
!A nonzero element was encountered, so from now on we will
!have to do the sums in the loop above.
              if (sums.ne.zero) then
                  ii=i                
              endif                
          endif
          b(i)=sums
      enddo 
!Now we do the backsubstitution, equation (2.3.7).
      do  i=n,1,-1                        
          sums=b(i)
          do  j=i+1,n
          sums=sums-a(i,j)*b(j)
          enddo 
!Store a component of the solution vector X.
          b(i)=sums/a(i,i)         
      enddo 
!All done!
      return                                 
      END  SUBROUTINE lubksb
!************************************************************ 
      ENd module obsemitter



