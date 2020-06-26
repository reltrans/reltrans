! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for \ell * g^{2+Gamma}
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) Ionisation profile
! 5) The dependence on incident angle (mimicked by tweaking the ionization)

!CURRENT BRANCH
! This branch is adding the xillverDCp table to go even higher densities.
  

  
include 'subroutines/header.h'


! !-----------------------------------------------------------------------
!       subroutine tdreltransD(ear,ne,param,ifl,photar)
!       implicit none
!       integer ne,ifl,Cp
!       real ear(0:ne),param(19),photar(ne)
!       Cp = 0
!       call genreltransD(Cp,ear,ne,param,ifl,photar)
!       return
!     end subroutine tdreltransD
! !-----------------------------------------------------------------------
        
! !-----------------------------------------------------------------------
!       subroutine tdreltrans(ear,ne,param,ifl,photar)
!       implicit none
!       integer ne,ifl,Cp
!       real ear(0:ne),param(19),photar(ne)
!       Cp = 0
!       call genreltrans(Cp,ear,ne,param,ifl,photar)
!       return
!     end subroutine tdreltrans
! !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine tdreltransDCp(ear,ne,param,ifl,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(20),photar(ne)
      Cp = 1
      call genreltransDCp(Cp,ear,ne,param,ifl,photar)
      return
    end subroutine tdreltransDCp
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine genreltransDCp(Cp,ear,ne,param,ifl,photar)
!!! ReltransD with xillverDCp
  !!! This routine calculates reltrans model starting form a pivoting
  !!! powerlaw and the rest-frame spectrum xillverDCp

!!! Arg:
  !!! 
        
  use dyn_gr
  use conv_mod
  implicit none

  double precision, parameter :: pi = acos(-1.d0)
  integer, intent(inout) :: ifl  
  integer nro,nphi, i,nf,ne,ReIm,nfsave
  integer verbose,mubin,rbin
  integer me,xe,Cp,j
  
  double precision a,h,Gamma,inc,rout,rmin,disco,muobs,rin
  double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
  double precision fhisave,flosave,rh,frrel,lens,fc
  real afac,param(20),ear(0:ne),gso
  real Afe,Ecut_s,Ecut_obs,logxi,lognep, xillparDCp(8),E,dE,earx(0:nex),Emax,Emin,dloge, kTe, kTe_s
  real reline(nex),imline(nex),photarx(nex)
  real reconvmu(nex),imconvmu(nex),mue,gsd
  real phase,ReS(ne),ImS(ne),photar(ne)
  real paramsave(20),contx(nex),absorbx(nex),photerx(nex)
  real ReGx(nex),ImGx(nex),Nh
  complex,dimension(:,:,:,:),allocatable :: transe,transea
  logical firstcall,needtrans,needconv, fftw
  integer myenv,Cpsave,gbin,check
  real, allocatable :: ReW0(:,:),ImW0(:,:),ReW1(:,:),ImW1(:,:),ReW2(:,:),ImW2(:,:)
  real, allocatable :: ReW3(:,:),ImW3(:,:)
  real, allocatable :: ReSraw(:,:),ImSraw(:,:),ReSrawa(:,:),ImSrawa(:,:)
  real, allocatable :: ReGrawa(:,:),ImGrawa(:,:),ReG(:,:),ImG(:,:)
  double precision, allocatable :: logxir(:),gsdr(:), logner(:)
  complex FTphotarx(4*nex),FTphotarx_delta(4*nex),FTreline(4*nex),FTimline(4*nex)
  complex FTreline_a(4*nex),FTimline_a(4*nex),FTreconv(4*nex),FTimconv(4*nex)
  complex FTphotarx_dlogxi(4*nex),sum
  real dyn,f,integral,phiA,logxi0,ImGbar(nex),ReGbar(nex),DelA,fhiHz,floHz,fac

! !variable for the grid reading
  integer :: irec,spin_dim,mu_dim
  double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      
!variable for non linear effects
  real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,DelAB,g
  real :: reline_a(nex),imline_a(nex),photarx_dlogxi(nex),dlogxi1,dlogxi2
  integer ionvar, DC, ionvariation
      
  data firstcall /.true./
  data Cpsave/2/
  data nfsave /-1/
  save firstcall,Emax,Emin,dloge,earx
  save lens,contx,me,xe
  save paramsave,fhisave,flosave,nfsave,nro,nphi
  save frobs,frrel,Cpsave
  save transe, transea, logxir, gsdr, logner
  save d,rnmax
  save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3
  save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG
  
  ifl = 1

! Parameter (DO SOMETHING WITH THIS)
  ionvar = 1
      
! Settings
  dlogf = 0.09 !0.0073  !This is a resolution parameter (base 10)
  dyn   = 1e-7
      
! Call environment variables
  verbose = myenv("REV_VERB",0)     !Set verbose level
      
! Initialise
  call initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans&
     ,nphi,nro,me,xe)

!Allocate dynamically the array to calculate the trasfer function          
  if (.not. allocated(re1)) allocate(re1(nphi,nro))
  if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
  if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

! Parameters
  h        = dble( param(1) )
  a        = dble( param(2) )
  inc      = dble( param(3) )
  rin      = dble( param(4) )
  rout     = dble( param(5) )
  zcos     = dble( param(6) )
  Gamma    = dble( param(7) )
  logxi    = param(8)
  Afe      = param(9)
  lognep   = param(10)
  kTe      = param(11)
  Nh       = param(12)
  afac     = param(13)
  Mass     = dble( param(14) )
  floHz    = param(15)
  fhiHz    = param(16)
  ReIm     = int( param(17) )
  DelA     = param(18)
  DelAB    = param(19)
  g        = param(20)

!work out the energy cut off for the powerlaw starting form kTe because we are using xillverDCp but also cut-off powerlaw for the pivoting effects
  Ecut_obs = 2.5 * kTe 

  honr = 0.d0
  muobs = cos( inc * pi / 180.d0 )
      
!Work out how many frequencies to average over
  fc = 0.5d0 * ( floHz + fhiHz )
  nf = ceiling( log10(fhiHz/floHz) / dlogf )
  if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
    fhiHz = 0.d0
    floHz = 0.d0
    nf    = 1
  end if
      
!Convert frequency bounds from Hz to c/Rg
  fhi   = dble(fhiHz) * 4.916d-6 * Mass
  flo   = dble(floHz) * 4.916d-6 * Mass

!Decide if this is the DC component or not
  if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
     DC     = 1
     ionvar = 0
     g      = 0.0
     DelAB  = 0.0
     DelA   = 0.0
     ReIm   = 1
  else
     DC     = 0
  end if
  
!Set minimum r (ISCO) and convert rin and h to rg
  if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
  rmin   = disco( a )
  if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
  rh     = 1.d0+sqrt(1.d0-a**2)
  if( h .lt. 0.d0 ) h = abs(h) * rh
  if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
  if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h
  if( rin .lt. rmin )then
     write(*,*)"Warning! rin<ISCO! Set to ISCO"
     rin = rmin
  end if
  if( h .lt. 1.5d0*rh )then
     write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
     h = 1.5d0 * rh
  end if

!Calculate source to observer g-factor and source frame Ecut
  gso    = real( dgsofac(a,h) )
  Ecut_s = real(1.d0+zcos) * Ecut_obs / gso
  kTe_s = real(1.d0+zcos) * kTe / gso
  ! if( verbose .gt. 0 )then
  !    if( Cp .eq. 0 )then
  !       write(*,*)"Ecut in source restframe (keV)=",Ecut_s
  !    else
  !       write(*,*)"kTe in source restframe (keV)=",Ecut_s
  !    end if
  ! end if
  
!Determine if I need to calculate the kernel
  if( .not. needtrans )then
     do i = 1,8
        if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
     end do
     if( abs( param(10) - paramsave(10) ) .gt. 1e-7 )needtrans=.true.
     ! if( abs( param(11) - paramsave(11) ) .gt. 1e-7 )needtrans=.true. 
     if( nf .ne. nfsave ) needtrans = .true.
     if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
     if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
   end if
   
! Allocate arrays that depend on frequency
  if( nf .ne. nfsave )then
     if( allocated(transe ) ) deallocate(transe )
     if( allocated(transea) ) deallocate(transea)
     allocate(  transe(nex,nf,me,xe) )
     allocate( transea(nex,nf,me,xe) )
     if( allocated(ReW0) ) deallocate(ReW0)
     if( allocated(ImW0) ) deallocate(ImW0)
     if( allocated(ReW1) ) deallocate(ReW1)
     if( allocated(ImW1) ) deallocate(ImW1)
     if( allocated(ReW2) ) deallocate(ReW2)
     if( allocated(ImW2) ) deallocate(ImW2)
     if( allocated(ReW3) ) deallocate(ReW3)
     if( allocated(ImW3) ) deallocate(ImW3)
     allocate( ReW0(nex,nf) )
     allocate( ImW0(nex,nf) )
     allocate( ReW1(nex,nf) )
     allocate( ImW1(nex,nf) )
     allocate( ReW2(nex,nf) )
     allocate( ImW2(nex,nf) )
     allocate( ReW3(nex,nf) )
     allocate( ImW3(nex,nf) )
     if( allocated(ReSraw) ) deallocate(ReSraw)
     if( allocated(ImSraw) ) deallocate(ImSraw)
     allocate( ReSraw(nex,nf) )
     allocate( ImSraw(nex,nf) )
     if( allocated(ReSrawa) ) deallocate(ReSrawa)
     if( allocated(ImSrawa) ) deallocate(ImSrawa)
     allocate( ReSrawa(nex,nf) )
     allocate( ImSrawa(nex,nf) )
     if( allocated(ReGrawa) ) deallocate(ReGrawa)
     if( allocated(ImGrawa) ) deallocate(ImGrawa)
     allocate( ReGrawa(nex,nf) )
     allocate( ImGrawa(nex,nf) )
     if( allocated(ReG) ) deallocate(ReG)
     if( allocated(ImG) ) deallocate(ImG)
     allocate( ReG(nex,nf) )
     allocate( ImG(nex,nf) )
  end if
     
  if( needtrans )then
     !Allocate arrays for kernels     
     if( .not. allocated(logxir) ) allocate( logxir(xe) )
     if( .not. allocated(gsdr)   ) allocate( gsdr  (xe) )
     if( .not. allocated(logner) ) allocate( logner(xe) )
     !Calculate the Kernel for the given parameters
     status_re_tau = .true.
     call rtrans(a,h,muobs,Gamma,rin,rout,honr,d,rnmax,zcos,nro,nphi,nex,dloge,&
          nf,fhi,flo,me,xe,logxi, lognep, transe,transea,frobs,frrel,lens,logxir,gsdr, logner)
  end if
  
  if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
  if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

!Determine if I need to convolve with the restframe reflection spectrum
  needconv = .false.
  if( needtrans ) needconv = .true.
  if( abs( param(9) - paramsave(9) ) .gt. 1e-7 ) needconv = .true.
  if( abs( param(11) - paramsave(11) ) .gt. 1e-7 ) needconv = .true.
  ! if( Cp .ne. Cpsave ) needconv = .true.
  
  if( needconv )then     
     !Initialize arrays for transfer functions
     ReW0 = 0.0
     ImW0 = 0.0
     ReW1 = 0.0
     ImW1 = 0.0
     ReW2 = 0.0
     ImW2 = 0.0
     ReW3 = 0.0
     ImW3 = 0.0
     DeltaGamma = 0.01
     Gamma1 = real(Gamma) - 0.5*DeltaGamma
     Gamma2 = real(Gamma) + 0.5*DeltaGamma

     !Get continuum spectrum 
     call getcontDCp(nex, earx, Gamma, Afe, kTe, lognep, logxi, contx, xillparDCp)
     ! if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
     
     !Get logxi values corresponding to Gamma1 and Gamma2
     call xilimits(nex,earx,contx,DeltaGamma,gso,real(zcos),dlogxi1,dlogxi2)

     !Now reflection
     xillparDCp(8) = -1.0       !reflection fraction of 1             

!Set the ion-variation to 1, there is an if inside the radial loop to check if either the ionvar is 0 or the logxi is 0 to set ionvariation to 0
! it is important that ionvariation is different than ionvar because ionvar is used also later in rawS routine to calculate the cross-spectrum
     ionvariation = 1

     write(*,*) 'RADIAL ZONES', xe
     write(*,*) 'MU ZONES', me
     !Loop over radius, emission angle and frequency
     do rbin = 1, xe  !Loop over radial zones

!Remember: xillparDCp(3) is kTe so we need to convert the Ecut_s into kTe
        
        xillparDCp(3) = real( gsdr(rbin) ) * kTe_s
        xillparDCp(4) = real( logner(rbin) )
        logxi0     = real( logxir(rbin) )
        if( xe .eq. 1 )then
           xillparDCp(3) = kTe_s
           xillparDCp(4) = lognep
           logxi0     = logxi
        end if
        
!Avoid negative values of the ionisation parameter 
        if (logxi0 .eq. 0.0 .or. ionvar .eq. 0) then
           ionvariation = 0.0
        endif

        do mubin = 1, me      !loop over emission angle zones
           !Calculate input inclination angle
           mue = ( real(mubin) - 0.5 ) / real(me)
           xillparDCp(7) = acos( mue ) * 180.0 / real(pi)
           if( me .eq. 1 ) xillparDCp(7) = real( inc )
           !Call xillver
           xillparDCp(1) = real(Gamma)
           xillparDCp(5) = logxi0
           call myxillDCp(earx, nex, xillparDCp, ifl, photarx)
           
           if (DC .eq. 0) then 
!NON LINEAR EFFECTS
              !Gamma variations
              xillparDCp(1) = Gamma1
              xillparDCp(5) = logxi0 + ionvariation * dlogxi1
              call myxillDCp(earx,nex,xillparDCp,ifl,photarx_1)
              xillparDCp(1) = Gamma2
              xillparDCp(5) = logxi0 + ionvariation * dlogxi2
              call myxillDCp(earx,nex,xillparDCp,ifl,photarx_2)
              photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
              !xi variations
              xillparDCp(1) = real(Gamma)
              xillparDCp(5) = logxi0 + ionvariation * dlogxi1
              call myxillDCp(earx,nex,xillparDCp,ifl,photarx_1)
              xillparDCp(1) = real(Gamma)
              xillparDCp(5) = logxi0 + ionvariation * dlogxi2
              call myxillDCp(earx,nex,xillparDCp,ifl,photarx_2)
              photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10           

        endif
        
           !Loop through frequencies
           do j = 1,nf
                 do i = 1,nex
                    reline(i)   = real(  transe(i,j,mubin,rbin) )
                    imline(i)   = aimag( transe(i,j,mubin,rbin) )
                    reline_a(i) = real(  transea(i,j,mubin,rbin) )
                    imline_a(i) = aimag( transea(i,j,mubin,rbin) )
                 end do
              
                 call conv_all_FFTw(dyn, photarx, photarx_delta, reline, imline, reline_a , imline_a,&
                         photarx_dlogxi, ReW0(:,j), ImW0(:,j), ReW1(:,j), ImW1(:,j), &
                         ReW2(:,j), ImW2(:,j), ReW3(:,j), ImW3(:,j), DC)

                 

                 end do !end of the frequency loop 

              end do
           end do

        end if

! Calculate raw FT of the full spectrum without absorption
  call rawS(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,afac,real(zcos),&
                gso,real(lens),real(Gamma),ionvar,DC,ReSraw,ImSraw)


! Calculate absorption and multiply by the raw FT
 ! call FNINIT

  call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
  
  do j = 1, nf
     do i = 1, nex
        ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
        ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
     end do
  end do

! Average over the frequency range
  if( DC .eq. 1 )then
     do i = 1, nex
        ReGbar(i) = ReSrawa(i,1)
!        ImGbar(i) = ImSrawa(i,1)  !No need for the immaginary part in DC
     end do
  else

     ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
     if (ReIm .gt. 0.0) then
        call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
     else
        call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
     endif
     
! Apply phase correction parameter to the cross-spectral model (for bad calibration)
     do j = 1,nf
        do i = 1,nex
           ReG(i,j) = cos(DelA) * ReGrawa(i,j) - sin(DelA) * ImGrawa(i,j)
           ImG(i,j) = cos(DelA) * ImGrawa(i,j) + sin(DelA) * ReGrawa(i,j)
        end do
     end do

     ReGbar = 0.0
     ImGbar = 0.0
     fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
     do j = 1,nf
        f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
        do i = 1,nex
           ReGbar(i) = ReGbar(i) + ReG(i,j) / f
           ImGbar(i) = ImGbar(i) + ImG(i,j) / f
        end do
     end do
     ReGbar = ReGbar * fac
     ImGbar = ImGbar * fac
  end if
     
! Write output depending on ReIm parameter
!  if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) ) ReIm = 1
  if( abs(ReIm) .le. 4 )then
     call crebin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is in photar form
     if( abs(ReIm) .eq. 1 )then        !Real part
        photar = ReS
     else if( abs(ReIm) .eq. 2 )then   !Imaginary part
        photar = ImS
     else if( abs(ReIm) .eq. 3 )then   !Modulus
        photar = sqrt( ReS**2 + ImS**2 )
        write(*,*) "Warning ReIm=3 should not be used for fitting!"
     else if( abs(ReIm) .eq. 4 )then   !Time lag (s)
        do i = 1,ne
           dE = ear(i) - ear(i-1)
           photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
        end do
        write(*,*)"Warning ReIm=4 should not be used for fitting!"
     end if
  else
     call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is count rate
     if( abs(ReIm) .eq. 5 )then        !Modulus
        do i = 1, ne
           dE = ear(i) - ear(i-1)
           photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
        end do
     else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
        do i = 1, ne
           dE = ear(i) - ear(i-1)
           photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
        end do
     end if
  end if
 
  fhisave   = fhi
  flosave   = flo
  nfsave    = nf
  paramsave = param
   ! Cpsave    = Cp

end subroutine genreltransDCp
!-----------------------------------------------------------------------



! !-----------------------------------------------------------------------
! subroutine genreltransD(Cp,ear,ne,param,ifl,photar)
!   use dyn_gr
!   use conv_mod
!   implicit none
!   integer nro,nphi, i,nf,ifl,ne,ReIm,nfsave
!   integer verbose,mubin,rbin
!   integer me,ge,xe,Cp,j
!   ! integer nex
!   ! parameter (nex=2**12)
!   double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
!   double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
!   double precision fhisave,flosave,rh,frrel,lens,fc
!   real afac,param(19),ear(0:ne),gso
!   real Afe,Ecut_s,Ecut_obs,logxi,lognep, xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
!   real reline(nex),imline(nex),photarx(nex)
!   real reconvmu(nex),imconvmu(nex),mue,gsd
!   real phase,ReS(ne),ImS(ne),photar(ne)
!   real paramsave(19),contx(nex),absorbx(nex),photerx(nex)
!   real ReGx(nex),ImGx(nex),Nh
!   complex,dimension(:,:,:,:),allocatable :: transe,transea
!   logical firstcall,needtrans,needconv, fftw
!   integer myenv,Cpsave,gbin,check
!   real, allocatable :: ReW0(:,:),ImW0(:,:),ReW1(:,:),ImW1(:,:),ReW2(:,:),ImW2(:,:)
!   real, allocatable :: ReW3(:,:),ImW3(:,:)
!   real, allocatable :: ReSraw(:,:),ImSraw(:,:),ReSrawa(:,:),ImSrawa(:,:)
!   real, allocatable :: ReGrawa(:,:),ImGrawa(:,:),ReG(:,:),ImG(:,:)
!   double precision, allocatable :: logxir(:),gsdr(:), logner(:)
!   complex FTphotarx(4*nex),FTphotarx_delta(4*nex),FTreline(4*nex),FTimline(4*nex)
!   complex FTreline_a(4*nex),FTimline_a(4*nex),FTreconv(4*nex),FTimconv(4*nex)
!   complex FTphotarx_dlogxi(4*nex),sum
!   real dyn,f,integral,phiA,logxi0,ImGbar(nex),ReGbar(nex),DelA,fhiHz,floHz,fac

! ! !variable for the grid reading
!   integer :: irec,spin_dim,mu_dim
!   double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      
! !variable for non linear effects
!   real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,DelAB,g
!   real :: reline_a(nex),imline_a(nex),photarx_dlogxi(nex),dlogxi1,dlogxi2
!   integer ionvar, DC, ionvariation
      
!   data firstcall /.true./
!   data Cpsave/2/
!   data nfsave /-1/
!   save firstcall,Emax,Emin,dloge,earx
!   save lens,contx,me,ge,xe
!   save paramsave,fhisave,flosave,nfsave,nro,nphi
!   save frobs,frrel,Cpsave
!   save transe, transea, logxir, gsdr, logner
!   save check,d,rnmax
!   save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,logxir,gsdr
!   save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG
  
!   pi = acos(-1.d0)
!   ifl = 1

! ! Parameter (DO SOMETHING WITH THIS)
!   ionvar = 1
      
! ! Settings
!   dlogf = 0.09 !0.0073  !This is a resolution parameter (base 10)
!   dyn   = 1e-7
      
! ! Call environment variables
!   verbose = myenv("REV_VERB",0)     !Set verbose level
      
! ! Initialise
!   call initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans,check&
!      ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim,me,ge,xe)

! !Allocate dynamically the array to calculate the trasfer function          
!   if (.not. allocated(re1)) allocate(re1(nphi,nro))
!   if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
!   if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

! ! Parameters
!   h        = dble( param(1) )
!   a        = dble( param(2) )
!   inc      = dble( param(3) )
!   rin      = dble( param(4) )
!   rout     = dble( param(5) )
!   zcos     = dble( param(6) )
!   Gamma    = dble( param(7) )
!   logxi    = param(8)
!   Afe      = param(9)
! !   Ecut_obs = param(10)
!   lognep   = param(10)
!   Nh       = param(11)
!   afac     = param(12)
!   Mass     = dble( param(13) )
!   floHz    = param(14)
!   fhiHz    = param(15)
!   ReIm     = int( param(16) )
!   DelA     = param(17)
!   DelAB    = param(18)
!   g        = param(19)

!   Ecut_obs = 300.0 ! It is frozen in xillverD 

!   honr = 0.d0
!   muobs = cos( inc * pi / 180.d0 )

! !check if the grid values are the same one of the model
!   if( check .ne. 0 .and. honr_grid .ne. honr ) then
!      write(*,*) 'grid has a different honr!'
!      write(*,*) 'honr of the grid is ', honr_grid
!      stop
!   endif
      
! !Work out how many frequencies to average over
!   fc = 0.5d0 * ( floHz + fhiHz )
!   nf = ceiling( log10(fhiHz/floHz) / dlogf )
!   if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
!     fhiHz = 0.d0
!     floHz = 0.d0
!     nf    = 1
!   end if
      
! !Convert frequency bounds from Hz to c/Rg
!   fhi   = dble(fhiHz) * 4.916d-6 * Mass
!   flo   = dble(floHz) * 4.916d-6 * Mass

! !Decide if this is the DC component or not
!   if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
!      DC     = 1
!      ionvar = 0
!      g      = 0.0
!      DelAB  = 0.0
!      DelA   = 0.0
!      ReIm   = 1
!   else
!      DC     = 0
!   end if
  
! !Set minimum r (ISCO) and convert rin and h to rg
!   if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
!   rmin   = disco( a )
!   if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
!   rh     = 1.d0+sqrt(1.d0-a**2)
!   if( h .lt. 0.d0 ) h = abs(h) * rh
!   if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
!   if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h
!   if( rin .lt. rmin )then
!      write(*,*)"Warning! rin<ISCO! Set to ISCO"
!      rin = rmin
!   end if
!   if( h .lt. 1.5d0*rh )then
!      write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
!      h = 1.5d0 * rh
!   end if

! !Calculate source to observer g-factor and source frame Ecut
!   gso    = real( dgsofac(a,h) )
!   Ecut_s = real(1.d0+zcos) * Ecut_obs / gso
!   if( verbose .gt. 0 )then
!      if( Cp .eq. 0 )then
!         write(*,*)"Ecut in source restframe (keV)=",Ecut_s
!      else
!         write(*,*)"kTe in source restframe (keV)=",Ecut_s
!      end if
!   end if
  
! !Determine if I need to calculate the kernel
!   if( .not. needtrans )then
!      do i = 1,8
!         if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
!      end do
!      if( abs( param(10) - paramsave(10) ) .gt. 1e-7 )needtrans=.true.
!      if( nf .ne. nfsave ) needtrans = .true.
!      if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
!      if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
!    end if
   
! ! Allocate arrays that depend on frequency
!   if( nf .ne. nfsave )then
!      if( allocated(transe ) ) deallocate(transe )
!      if( allocated(transea) ) deallocate(transea)
!      allocate(  transe(nex,nf,me,xe) )
!      allocate( transea(nex,nf,me,xe) )
!      if( allocated(ReW0) ) deallocate(ReW0)
!      if( allocated(ImW0) ) deallocate(ImW0)
!      if( allocated(ReW1) ) deallocate(ReW1)
!      if( allocated(ImW1) ) deallocate(ImW1)
!      if( allocated(ReW2) ) deallocate(ReW2)
!      if( allocated(ImW2) ) deallocate(ImW2)
!      if( allocated(ReW3) ) deallocate(ReW3)
!      if( allocated(ImW3) ) deallocate(ImW3)
!      allocate( ReW0(nex,nf) )
!      allocate( ImW0(nex,nf) )
!      allocate( ReW1(nex,nf) )
!      allocate( ImW1(nex,nf) )
!      allocate( ReW2(nex,nf) )
!      allocate( ImW2(nex,nf) )
!      allocate( ReW3(nex,nf) )
!      allocate( ImW3(nex,nf) )
!      if( allocated(ReSraw) ) deallocate(ReSraw)
!      if( allocated(ImSraw) ) deallocate(ImSraw)
!      allocate( ReSraw(nex,nf) )
!      allocate( ImSraw(nex,nf) )
!      if( allocated(ReSrawa) ) deallocate(ReSrawa)
!      if( allocated(ImSrawa) ) deallocate(ImSrawa)
!      allocate( ReSrawa(nex,nf) )
!      allocate( ImSrawa(nex,nf) )
!      if( allocated(ReGrawa) ) deallocate(ReGrawa)
!      if( allocated(ImGrawa) ) deallocate(ImGrawa)
!      allocate( ReGrawa(nex,nf) )
!      allocate( ImGrawa(nex,nf) )
!      if( allocated(ReG) ) deallocate(ReG)
!      if( allocated(ImG) ) deallocate(ImG)
!      allocate( ReG(nex,nf) )
!      allocate( ImG(nex,nf) )
!   end if
     
!   if( needtrans )then
!      !Allocate arrays for kernels     
!      if( .not. allocated(logxir) ) allocate( logxir(xe) )
!      if( .not. allocated(gsdr)   ) allocate( gsdr  (xe) )
!      if( .not. allocated(logner) ) allocate( logner(xe) )
!      !Calculate the Kernel for the given parameters
!      status_re_tau = .true.
!      call rtrans(a,h,muobs,Gamma,rin,rout,honr,d,rnmax,zcos,nro,nphi,nex,dloge,&
!           nf,fhi,flo,me,xe,logxi, lognep, transe,transea,frobs,frrel,lens,logxir,gsdr, logner)
!   end if
  
!   if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
!   if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

! !Determine if I need to convolve with the restframe reflection spectrum
!   needconv = .false.
!   if( needtrans ) needconv = .true.
!   do i = 8,10
!     if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needconv = .true.
!   end do
!   if( Cp .ne. Cpsave ) needconv = .true.
  
!   if( needconv )then     
!      !Initialize arrays for transfer functions
!      ReW0 = 0.0
!      ImW0 = 0.0
!      ReW1 = 0.0
!      ImW1 = 0.0
!      ReW2 = 0.0
!      ImW2 = 0.0
!      ReW3 = 0.0
!      ImW3 = 0.0
!      DeltaGamma = 0.01
!      Gamma1 = real(Gamma) - 0.5*DeltaGamma
!      Gamma2 = real(Gamma) + 0.5*DeltaGamma
!      !Get continuum spectrum (D is the high density xillver)
! !      call getcont(nex,earx,Gamma,Afe,Ecut_obs,logxi,Cp,contx,xillpar)
!      call getcontD(nex, earx, Gamma, Afe, lognep, logxi, Cp, contx, xillpar)
!      if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
!      !Get logxi values corresponding to Gamma1 and Gamma2
!      call xilimits(nex,earx,contx,DeltaGamma,gso,real(zcos),dlogxi1,dlogxi2)
!      !Now reflection
!      xillpar(7) = -1.0       !reflection fraction of 1             

! !Set the ion-variation to 1, there is an if inside the radial loop to check if either the ionvar is 0 or the logxi is 0 to set ionvariation to 0
! ! it is important that ionvariation is different than ionvar because ionvar is used also later in rawS routine to calculate the cross-spectrum
!      ionvariation = 1

!      !Loop over radius, emission angle and frequency
!      do rbin = 1, xe  !Loop over radial zones

! !Remember: xillpar(3) is Ecut when xillver is with density fixed to 10^15 
! !          whereas in xillverD the parameter 3 is logN 

! !       xillpar(3) = real( gsdr(rbin) ) * Ecut_s
!         xillpar(3) = real( logner(rbin) )
!         logxi0     = real( logxir(rbin) )
!         if( xe .eq. 1 )then
! !            xillpar(3) = Ecut_s
!            xillpar(3) = lognep
!            logxi0     = logxi
!         end if
        
! !Avoid negative values of the ionisation parameter 
!         if (logxi0 .eq. 0.0 .or. ionvar .eq. 0) then
!            ionvariation = 0.0
!         endif

!         do mubin = 1,me      !loop over emission angle zones
!            !Calculate input inclination angle
!            mue = ( real(mubin) - 0.5 ) / real(me)
!            xillpar(6) = acos( mue ) * 180.0 / real(pi)
!            if( me .eq. 1 ) xillpar(6) = real( inc )
!            !Call xillver
!            xillpar(1) = real(Gamma)
!            xillpar(4) = logxi0
!           ! write(*,*) xillpar(4), xillpar(3)            
! !           call myxill   (earx,nex,xillpar,ifl,Cp,photarx)
!            call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx)
           
!            if (DC .eq. 0) then 
! !NON LINEAR EFFECTS
!               !Gamma variations
!               xillpar(1) = Gamma1
!               xillpar(4) = logxi0 + ionvariation * dlogxi1
!               ! call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
!               call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_1)
!               xillpar(1) = Gamma2
!               xillpar(4) = logxi0 + ionvariation * dlogxi2
!               ! call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
!               call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_2)
!               photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
!               !xi variations
!               xillpar(1) = real(Gamma)
!               xillpar(4) = logxi0 + ionvariation * dlogxi1
!               ! call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
!               call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_1)
!               xillpar(1) = real(Gamma)
!               xillpar(4) = logxi0 + ionvariation * dlogxi2
!               ! call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
!               call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_2)
!               photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10           

!         endif
        
!            !Loop through frequencies
!            do j = 1,nf
!                  do i = 1,nex
!                     reline(i)   = real(  transe(i,j,mubin,rbin) )
!                     imline(i)   = aimag( transe(i,j,mubin,rbin) )
!                     reline_a(i) = real(  transea(i,j,mubin,rbin) )
!                     imline_a(i) = aimag( transea(i,j,mubin,rbin) )
!                  end do
              
!                  call conv_all_FFTw(dyn, photarx, photarx_delta, reline, imline, reline_a , imline_a,&
!                          photarx_dlogxi, ReW0(:,j), ImW0(:,j), ReW1(:,j), ImW1(:,j), &
!                          ReW2(:,j), ImW2(:,j), ReW3(:,j), ImW3(:,j), DC)

                 
! !*************************************************************************!
! !Convolve with line profile. This convolution is done with four1
!                     ! !First FFTs
!                     ! call pad4FFT(nex,photarx,FTphotarx)
!                     ! call pad4FFT(nex,reline,FTreline)
!                     ! call pad4FFT(nex,imline,FTimline)
!                     ! call pad4FFT(nex,photarx_delta,FTphotarx_delta)
!                     ! call pad4FFT(nex,reline_a,FTreline_a)
!                     ! call pad4FFT(nex,imline_a,FTimline_a)
!                     ! call pad4FFT(nex,photarx_dlogxi,FTphotarx_dlogxi)

!                     ! !Then the multiplications and inverse FFTs
!                     ! FTreconv = FTreline * FTphotarx
!                     ! FTimconv = FTimline * FTphotarx
!                     ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
!                     ! call pad4invFFT(dyn,nex,FTimconv,imconvmu) 
!                     ! do i = 1,nex
!                     !    ReW0(i,j) = ReW0(i,j) + reconvmu(i)
!                     !    ImW0(i,j) = ImW0(i,j) + imconvmu(i)
!                     ! end do

!                     ! FTreconv = FTreline_a * FTphotarx
!                     ! FTimconv = FTimline_a * FTphotarx
!                     ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
!                     ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
!                     ! do i = 1,nex
!                     !    ReW1(i,j) = ReW1(i,j) + reconvmu(i)
!                     !    ImW1(i,j) = ImW1(i,j) + imconvmu(i)
!                     ! end do
!                     ! FTreconv = FTreline * FTphotarx_delta
!                     ! FTimconv = FTimline * FTphotarx_delta
!                     ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
!                     ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
!                     ! do i = 1,nex
!                     !    ReW2(i,j) = ReW2(i,j) + reconvmu(i)
!                     !    ImW2(i,j) = ImW2(i,j) + imconvmu(i)
!                     ! end do
!                     ! FTreconv = FTreline * FTphotarx_dlogxi
!                     ! FTimconv = FTimline * FTphotarx_dlogxi
!                     ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
!                     ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
!                     ! do i = 1,nex
!                     !    ReW3(i,j) = ReW3(i,j) + reconvmu(i)
!                     !    ImW3(i,j) = ImW3(i,j) + imconvmu(i)
!                     ! end do
! !*************************************************************************!

!                  end do !end of the frequency loop 

!               end do
!            end do

!         end if

! ! Calculate raw FT of the full spectrum without absorption
!   call rawS(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,afac,real(zcos),&
!                 gso,real(lens),real(Gamma),ionvar,DC,ReSraw,ImSraw)


! ! Calculate absorption and multiply by the raw FT
!  ! call FNINIT

!   call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
  
!   do j = 1, nf
!      do i = 1, nex
!         ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
!         ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
!      end do
!   end do

! ! Average over the frequency range
!   if( DC .eq. 1 )then
!      do i = 1, nex
!         ReGbar(i) = ReSrawa(i,1)
! !        ImGbar(i) = ImSrawa(i,1)  !No need for the immaginary part in DC
!      end do
!   else

!      ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
!      if (ReIm .gt. 0.0) then
!         call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
!      else
!         call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
!      endif
     
! ! Apply phase correction parameter to the cross-spectral model (for bad calibration)
!      do j = 1,nf
!         do i = 1,nex
!            ReG(i,j) = cos(DelA) * ReGrawa(i,j) - sin(DelA) * ImGrawa(i,j)
!            ImG(i,j) = cos(DelA) * ImGrawa(i,j) + sin(DelA) * ReGrawa(i,j)
!         end do
!      end do

!      ReGbar = 0.0
!      ImGbar = 0.0
!      fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
!      do j = 1,nf
!         f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
!         do i = 1,nex
!            ReGbar(i) = ReGbar(i) + ReG(i,j) / f
!            ImGbar(i) = ImGbar(i) + ImG(i,j) / f
!         end do
!      end do
!      ReGbar = ReGbar * fac
!      ImGbar = ImGbar * fac
!   end if
     
! ! Write output depending on ReIm parameter
! !  if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) ) ReIm = 1
!   if( abs(ReIm) .le. 4 )then
!      call crebin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is in photar form
!      if( abs(ReIm) .eq. 1 )then        !Real part
!         photar = ReS
!      else if( abs(ReIm) .eq. 2 )then   !Imaginary part
!         photar = ImS
!      else if( abs(ReIm) .eq. 3 )then   !Modulus
!         photar = sqrt( ReS**2 + ImS**2 )
!         write(*,*) "Warning ReIm=3 should not be used for fitting!"
!      else if( abs(ReIm) .eq. 4 )then   !Time lag (s)
!         do i = 1,ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
!         end do
!         write(*,*)"Warning ReIm=4 should not be used for fitting!"
!      end if
!   else
!      call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is count rate
!      if( abs(ReIm) .eq. 5 )then        !Modulus
!         do i = 1, ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
!         end do
!      else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
!         do i = 1, ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
!         end do
!      end if
!   end if
 
!   fhisave   = fhi
!   flosave   = flo
!   nfsave    = nf
!   paramsave = param
!   Cpsave    = Cp

! end subroutine genreltransD
! !-----------------------------------------------------------------------



! !-----------------------------------------------------------------------
! subroutine genreltrans(Cp,ear,ne,param,ifl,photar)
!   use dyn_gr
!   use conv_mod
!   implicit none
!   integer nro,nphi, i,nf,ifl,ne,ReIm,nfsave
!   integer verbose,mubin,rbin
!   integer me,ge,xe,Cp,j
!   ! integer nex
!   ! parameter (nex=2**12)
!   double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
!   double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
!   double precision fhisave,flosave,rh,frrel,lens,fc
!   real afac,param(19),ear(0:ne),gso
!   real Afe,Ecut_s,Ecut_obs,logxi,lognep, xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
!   real reline(nex),imline(nex),photarx(nex)
!   real reconvmu(nex),imconvmu(nex),mue,gsd
!   real phase,ReS(ne),ImS(ne),photar(ne)
!   real paramsave(19),contx(nex),absorbx(nex),photerx(nex)
!   real ReGx(nex),ImGx(nex),Nh
!   complex,dimension(:,:,:,:),allocatable :: transe,transea
!   logical firstcall,needtrans,needconv, fftw
!   integer myenv,Cpsave,gbin,check
!   real, allocatable :: ReW0(:,:),ImW0(:,:),ReW1(:,:),ImW1(:,:),ReW2(:,:),ImW2(:,:)
!   real, allocatable :: ReW3(:,:),ImW3(:,:)
!   real, allocatable :: ReSraw(:,:),ImSraw(:,:),ReSrawa(:,:),ImSrawa(:,:)
!   real, allocatable :: ReGrawa(:,:),ImGrawa(:,:),ReG(:,:),ImG(:,:)
!   double precision, allocatable :: logxir(:),gsdr(:), logner(:)
!   complex FTphotarx(4*nex),FTphotarx_delta(4*nex),FTreline(4*nex),FTimline(4*nex)
!   complex FTreline_a(4*nex),FTimline_a(4*nex),FTreconv(4*nex),FTimconv(4*nex)
!   complex FTphotarx_dlogxi(4*nex),sum
!   real dyn,f,integral,phiA,logxi0,ImGbar(nex),ReGbar(nex),DelA,fhiHz,floHz,fac

! ! !variable for the grid reading
!   integer :: irec,spin_dim,mu_dim
!   double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      
! !variable for non linear effects
!   real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,DelAB,g
!   real :: reline_a(nex),imline_a(nex),photarx_dlogxi(nex),dlogxi1,dlogxi2
!   integer ionvar,DC
      
!   data firstcall /.true./
!   data Cpsave/2/
!   data nfsave /-1/
!   save firstcall,Emax,Emin,dloge,earx
!   save lens,contx,me,ge,xe
!   save paramsave,fhisave,flosave,nfsave,nro,nphi
!   save frobs,frrel,Cpsave
!   save transe,transea
!   save check,d,rnmax
!   save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,logxir,gsdr
!   save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG
  
!   pi = acos(-1.d0)
!   ifl = 1

! ! Parameter (DO SOMETHING WITH THIS)
!   ionvar = 1
      
! ! Settings
!   dlogf = 0.09 !0.0073  !This is a resolution parameter (base 10)
!   dyn   = 1e-7
      
! ! Call environment variables
!   verbose = myenv("REV_VERB",0)     !Set verbose level
      
! ! Initialise
!   call initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans,check&
!      ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim,me,ge,xe)

! !Allocate dynamically the array to calculate the trasfer function          
!   if (.not. allocated(re1)) allocate(re1(nphi,nro))
!   if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
!   if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

! ! Parameters
!   h        = dble( param(1) )
!   a        = dble( param(2) )
!   inc      = dble( param(3) )
!   rin      = dble( param(4) )
!   rout     = dble( param(5) )
!   zcos     = dble( param(6) )
!   Gamma    = dble( param(7) )
!   logxi    = param(8)
!   Afe      = param(9)
!   Ecut_obs = param(10)
!   Nh       = param(11)
!   afac     = param(12)
!   Mass     = dble( param(13) )
!   floHz    = param(14)
!   fhiHz    = param(15)
!   ReIm     = int( param(16) )
!   DelA     = param(17)
!   DelAB    = param(18)
!   g        = param(19)

 
!   honr = 0.d0
!   muobs = cos( inc * pi / 180.d0 )

! !check if the grid values are the same one of the model
!   if( check .ne. 0 .and. honr_grid .ne. honr ) then
!      write(*,*) 'grid has a different honr!'
!      write(*,*) 'honr of the grid is ', honr_grid
!      stop
!   endif
      
! !Work out how many frequencies to average over
!   fc = 0.5d0 * ( floHz + fhiHz )
!   nf = ceiling( log10(fhiHz/floHz) / dlogf )
!   if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
!     fhiHz = 0.d0
!     floHz = 0.d0
!     nf    = 1
!   end if
      
! !Convert frequency bounds from Hz to c/Rg
!   fhi   = dble(fhiHz) * 4.916d-6 * Mass
!   flo   = dble(floHz) * 4.916d-6 * Mass

! !Decide if this is the DC component or not
!   if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
!      DC     = 1
!      ionvar = 0
!      g      = 0.0
!      DelAB  = 0.0
!      DelA   = 0.0
!      ReIm   = 1
!   else
!      DC     = 0
!   end if
  
! !Set minimum r (ISCO) and convert rin and h to rg
!   if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
!   rmin   = disco( a )
!   if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
!   rh     = 1.d0+sqrt(1.d0-a**2)
!   if( h .lt. 0.d0 ) h = abs(h) * rh
!   if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
!   if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h
!   if( rin .lt. rmin )then
!      write(*,*)"Warning! rin<ISCO! Set to ISCO"
!      rin = rmin
!   end if
!   if( h .lt. 1.5d0*rh )then
!      write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
!      h = 1.5d0 * rh
!   end if

! !Calculate source to observer g-factor and source frame Ecut
!   gso    = real( dgsofac(a,h) )
!   Ecut_s = real(1.d0+zcos) * Ecut_obs / gso
!   if( verbose .gt. 0 )then
!      if( Cp .eq. 0 )then
!         write(*,*)"Ecut in source restframe (keV)=",Ecut_s
!      else
!         write(*,*)"kTe in source restframe (keV)=",Ecut_s
!      end if
!   end if
  
! !Determine if I need to calculate the kernel
!   if( .not. needtrans )then
!      do i = 1,8
!         if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
!      end do
!      if( nf .ne. nfsave ) needtrans = .true.
!      if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
!      if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
!    end if
   
! ! Allocate arrays that depend on frequency
!   if( nf .ne. nfsave )then
!      if( allocated(transe ) ) deallocate(transe )
!      if( allocated(transea) ) deallocate(transea)
!      allocate(  transe(nex,nf,me,xe) )
!      allocate( transea(nex,nf,me,xe) )
!      if( allocated(ReW0) ) deallocate(ReW0)
!      if( allocated(ImW0) ) deallocate(ImW0)
!      if( allocated(ReW1) ) deallocate(ReW1)
!      if( allocated(ImW1) ) deallocate(ImW1)
!      if( allocated(ReW2) ) deallocate(ReW2)
!      if( allocated(ImW2) ) deallocate(ImW2)
!      if( allocated(ReW3) ) deallocate(ReW3)
!      if( allocated(ImW3) ) deallocate(ImW3)
!      allocate( ReW0(nex,nf) )
!      allocate( ImW0(nex,nf) )
!      allocate( ReW1(nex,nf) )
!      allocate( ImW1(nex,nf) )
!      allocate( ReW2(nex,nf) )
!      allocate( ImW2(nex,nf) )
!      allocate( ReW3(nex,nf) )
!      allocate( ImW3(nex,nf) )
!      if( allocated(ReSraw) ) deallocate(ReSraw)
!      if( allocated(ImSraw) ) deallocate(ImSraw)
!      allocate( ReSraw(nex,nf) )
!      allocate( ImSraw(nex,nf) )
!      if( allocated(ReSrawa) ) deallocate(ReSrawa)
!      if( allocated(ImSrawa) ) deallocate(ImSrawa)
!      allocate( ReSrawa(nex,nf) )
!      allocate( ImSrawa(nex,nf) )
!      if( allocated(ReGrawa) ) deallocate(ReGrawa)
!      if( allocated(ImGrawa) ) deallocate(ImGrawa)
!      allocate( ReGrawa(nex,nf) )
!      allocate( ImGrawa(nex,nf) )
!      if( allocated(ReG) ) deallocate(ReG)
!      if( allocated(ImG) ) deallocate(ImG)
!      allocate( ReG(nex,nf) )
!      allocate( ImG(nex,nf) )
!   end if
     
!   if( needtrans )then
!      !Allocate arrays for kernels     
!      if( .not. allocated(logxir) ) allocate( logxir(xe) )
!      if( .not. allocated(gsdr)   ) allocate( gsdr  (xe) )
!      if( .not. allocated(logner) ) allocate( logner(xe) )
!      !Calculate the Kernel for the given parameters
!      status_re_tau = .true.
!      call rtrans(a, h, muobs, Gamma, rin, rout, honr, d, rnmax,&
!           zcos, nro, nphi, nex, dloge, nf, fhi, flo, me, xe, &
!           logxi, lognep, transe, transea, frobs, frrel, lens, &
!           logxir, gsdr, logner)
!   end if
  
!   if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
!   if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

! !Determine if I need to convolve with the restframe reflection spectrum
!   needconv = .false.
!   if( needtrans ) needconv = .true.
!   do i = 8,10
!     if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needconv = .true.
!   end do
!   if( Cp .ne. Cpsave ) needconv = .true.
  
!   if( needconv )then     
!      !Initialize arrays for transfer functions
!      ReW0 = 0.0
!      ImW0 = 0.0
!      ReW1 = 0.0
!      ImW1 = 0.0
!      ReW2 = 0.0
!      ImW2 = 0.0
!      ReW3 = 0.0
!      ImW3 = 0.0
!      DeltaGamma = 0.01
!      Gamma1 = real(Gamma) - 0.5*DeltaGamma
!      Gamma2 = real(Gamma) + 0.5*DeltaGamma

!      call getcont(nex, earx, Gamma, Afe, Ecut_obs, logxi, Cp, contx, xillpar)
!      if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
!      !Get logxi values corresponding to Gamma1 and Gamma2
!      call xilimits(nex,earx,contx,DeltaGamma,gso,real(zcos),dlogxi1,dlogxi2)
!      !Now reflection
!      xillpar(7) = -1.0       !reflection fraction of 1             

!      !Loop over radius, emission angle and frequency
!      do rbin = 1, xe  !Loop over radial zones

!         xillpar(3) = real( gsdr(rbin) ) * Ecut_s
!         logxi0     = real( logxir(rbin) )
!         if( xe .eq. 1 )then
!            xillpar(3) = Ecut_s
!            logxi0     = logxi
!         end if
        
! !Avoid negative values of the ionisation parameter 
!         do mubin = 1,me      !loop over emission angle zones
!            !Calculate input inclination angle
!            mue = ( real(mubin) - 0.5 ) / real(me)
!            xillpar(6) = acos( mue ) * 180.0 / real(pi)
!            if( me .eq. 1 ) xillpar(6) = real( inc )
!            !Call xillver
!            xillpar(1) = real(Gamma)
!            xillpar(4) = logxi0
! !           write(*,*) xillpar            
!            call myxill   (earx,nex,xillpar,ifl,Cp,photarx)
! !           call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx)
           
!            if (DC .eq. 0) then 
! !NON LINEAR EFFECTS
!               !Gamma variations
!               xillpar(1) = Gamma1
!               xillpar(4) = logxi0 + ionvar * dlogxi1
!               call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
! !              call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_1)
!               xillpar(1) = Gamma2
!               xillpar(4) = logxi0 + ionvar * dlogxi2
!               call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
! !              call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_2)
!               photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
!               !xi variations
!               xillpar(1) = real(Gamma)
!               xillpar(4) = logxi0 + ionvar * dlogxi1
!               call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
! !              call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_1)
!               xillpar(1) = real(Gamma)
!               xillpar(4) = logxi0 + ionvar * dlogxi2
!               call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
! !              call myxill_hD(earx,nex,xillpar,ifl,Cp,photarx_2)
!               photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10           

!         endif
        
!            !Loop through frequencies
!            do j = 1,nf
!                  do i = 1,nex
!                     reline(i)   = real(  transe(i,j,mubin,rbin) )
!                     imline(i)   = aimag( transe(i,j,mubin,rbin) )
!                     reline_a(i) = real(  transea(i,j,mubin,rbin) )
!                     imline_a(i) = aimag( transea(i,j,mubin,rbin) )
!                  end do
              
!                  call conv_all_FFTw(dyn, photarx, photarx_delta, reline, imline, reline_a , imline_a,&
!                          photarx_dlogxi, ReW0(:,j), ImW0(:,j), ReW1(:,j), ImW1(:,j), &
!                          ReW2(:,j), ImW2(:,j), ReW3(:,j), ImW3(:,j), DC)

                 
!                  end do !end of the frequency loop 

!               end do
!            end do

!         end if

! ! Calculate raw FT of the full spectrum without absorption
!   call rawS(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,afac,real(zcos),&
!                 gso,real(lens),real(Gamma),ionvar,DC,ReSraw,ImSraw)


! ! Calculate absorption and multiply by the raw FT
! !  call FNINIT

!   call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
  
!   do j = 1, nf
!      do i = 1, nex
!         ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
!         ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
!      end do
!   end do

! ! Average over the frequency range
!   if( DC .eq. 1 )then
!      do i = 1, nex
!         ReGbar(i) = ReSrawa(i,1)
! !        ImGbar(i) = ImSrawa(i,1)  !No need for the immaginary part in DC
!      end do
!   else

!      ! Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
!      if (ReIm .gt. 0.0) then
!         call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
!      else
!         call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
!      endif
     
! ! Apply phase correction parameter to the cross-spectral model (for bad calibration)
!      do j = 1,nf
!         do i = 1,nex
!            ReG(i,j) = cos(DelA) * ReGrawa(i,j) - sin(DelA) * ImGrawa(i,j)
!            ImG(i,j) = cos(DelA) * ImGrawa(i,j) + sin(DelA) * ReGrawa(i,j)
!         end do
!      end do

!      ReGbar = 0.0
!      ImGbar = 0.0
!      fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
!      do j = 1,nf
!         f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
!         do i = 1,nex
!            ReGbar(i) = ReGbar(i) + ReG(i,j) / f
!            ImGbar(i) = ImGbar(i) + ImG(i,j) / f
!         end do
!      end do
!      ReGbar = ReGbar * fac
!      ImGbar = ImGbar * fac
!   end if
     
! ! Write output depending on ReIm parameter
! !  if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) ) ReIm = 1
!   if( abs(ReIm) .le. 4 )then
!      call crebin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is in photar form
!      if( abs(ReIm) .eq. 1 )then        !Real part
!         photar = ReS
!      else if( abs(ReIm) .eq. 2 )then   !Imaginary part
!         photar = ImS
!      else if( abs(ReIm) .eq. 3 )then   !Modulus
!         photar = sqrt( ReS**2 + ImS**2 )
!         write(*,*) "Warning ReIm=3 should not be used for fitting!"
!      else if( abs(ReIm) .eq. 4 )then   !Time lag (s)
!         do i = 1,ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
!         end do
!         write(*,*)"Warning ReIm=4 should not be used for fitting!"
!      end if
!   else
!      call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is count rate
!      if( abs(ReIm) .eq. 5 )then        !Modulus
!         do i = 1, ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
!         end do
!      else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
!         do i = 1, ne
!            dE = ear(i) - ear(i-1)
!            photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
!         end do
!      end if
!   end if
 
!   fhisave   = fhi
!   flosave   = flo
!   nfsave    = nf
!   paramsave = param
!   Cpsave    = Cp

! end subroutine genreltrans
! !-----------------------------------------------------------------------




! *** Old convolution code
! This is slower because it does FFTs that don't need to be done.
! The new way takes ~68% of the time.
              !Convolve with line profile

              ! call padcnv(1e-7,nex,reline,photarx,reconvmu)
              ! call padcnv(1e-7,nex,imline,photarx,imconvmu)
              ! !Add to running sum
              ! do i = 1,nex
              !    ReW0(i,j) = ReW0(i,j) + reconvmu(i)
              !    ImW0(i,j) = ImW0(i,j) + imconvmu(i)
              ! end do

              ! call padcnv(1e-7,nex,reline_a,photarx,reconvmu)
              ! call padcnv(1e-7,nex,imline_a,photarx,imconvmu)
              ! !Add to running sum
              ! do i = 1,nex
              !    ReW1(i,j) = ReW1(i,j) + reconvmu(i)
              !    ImW1(i,j) = ImW1(i,j) + imconvmu(i)
              ! end do
              
              ! call padcnv(1e-7,nex,reline,photarx_delta,reconvmu)
              ! call padcnv(1e-7,nex,imline,photarx_delta,imconvmu)
              ! !Add to running sum
              ! do i = 1,nex
              !    ReW2(i,j) = ReW2(i,j) + reconvmu(i)
              !    ImW2(i,j) = ImW2(i,j) + imconvmu(i)
              ! end do

              ! call padcnv(1e-7,nex,reline,photarx_dlogxi,reconvmu)
              ! call padcnv(1e-7,nex,imline,photarx_dlogxi,imconvmu)
              ! !Add to running sum
              ! do i = 1,nex
              !    ReW3(i,j) = ReW3(i,j) + reconvmu(i)
              !    ImW3(i,j) = ImW3(i,j) + imconvmu(i)
              ! end do
    

!-----------------------------------------------------------------------
subroutine xilimits(nex,earx,contx,DeltaGamma,gso,z,dlogxi1,dlogxi2)
! Inputs: nex,earx,contx,DeltaGamma,gso
! Outputs: dlogxi1,dlogxi2
! logxi1 = logxi0 + dlogxi1
  implicit none
  integer nex,i
  real earx(0:nex),contx(nex),DeltaGamma,gso,z,dlogxi1,dlogxi2
  real num1,num2,den,logS1,logS2,gsoz,E
  gsoz = gso / (1.0+z)
  num1 = 0.0
  num2 = 0.0
  den = 0.0
  do i = 1,nex
     E   = 0.5 * ( earx(i) + earx(i-1) )
     num2 = num2 + E**(1.0-0.5*DeltaGamma) * contx(i)
     num1 = num1 + E**(1.0+0.5*DeltaGamma) * contx(i)
     den  = den  + E * contx(i)
  end do
  logS1 = log10(num1/den)
  logS2 = log10(num2/den)
  dlogxi1 = -0.5*DeltaGamma * log10(gsoz) + logS1
  dlogxi2 =  0.5*DeltaGamma * log10(gsoz) + logS2
  return
end subroutine xilimits
!-----------------------------------------------------------------------  
      

! !-----------------------------------------------------------------------
!       subroutine genreltrans(Cp,ear,ne,param,ifl,photar)
!       use dyn_gr
!       implicit none
!       integer nro,nphi,nex,i,nf,ifl,ne,ReIm,nfsave
!       integer verbose,mubin,sdbin
!       integer me,ge,xe,Cp
!       parameter (nex=2**12)
!       double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
!       double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
!       double precision fhisave,flosave,rh,frrel,lens,fc
!       real afac,param(19),ear(0:ne),gso,direct,ximin,ximax
!       real ldir,ReW,ImW,ReW1,ImW1
!       real Afe,Ecut_s,Ecut_obs,logxi,xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
!       real reline(nex),imline(nex),photarx(nex),reconv(nex),imconv(nex)
!       real reconvmu(nex),imconvmu(nex),mue,sdmin,sdmax,gsd
!       real phase,ReS(ne),ImS(ne),photar(ne)
!       real paramsave(19),contx(nex),phiA,absorbx(nex),photerx(nex)
!       real ReGx(nex),ImGx(nex),Nh
!       real contxabs(nex),reconvabs(nex),imconvabs(nex)
!       complex,dimension(:,:,:,:,:),allocatable :: transe,transea
!       logical firstcall,needtrans,needconv,needresp
!       integer xbin,xbinhi,myenv,Cpsave,gbin
!       real logxir,ghi,glo

! ! !variable for the grid reading
! !       integer :: irec,spin_dim,mu_dim,check
! !       integer :: s_lo,s_hi,m_lo,m_hi,xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh
! !       double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
! !       real :: sdmin_sl_ml,sdmax_sl_ml,sdmin_sl_mh,sdmax_sl_mh,sdmin_sh_ml,sdmax_sh_ml&
! !            ,sdmin_sh_mh,sdmax_sh_mh
! !       double precision :: frobs_sl_ml,frrel_sl_ml,frobs_sl_mh,frrel_sl_mh&
! !            ,frobs_sh_ml,frrel_sh_ml,frobs_sh_mh,frrel_sh_mh
! !       complex :: transe_1(nex),transe_2(nex),transe_(nex)
! !       complex,dimension(:,:,:,:),allocatable :: transe_11,transe_12,transe_21,transe_22
      
! !variable for non linear effects
!       complex :: transe_1a(nex),transe_2a(nex),transe_a(nex)
!       complex,dimension(:,:,:,:),allocatable :: transe_11a,transe_12a,transe_21a,transe_22a
!       real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,phiB,g
!       real :: reline_a(nex),imline_a(nex),reconvW1(nex),imconvW1(nex),reconvW1a(nex),imconvW1a(nex)
!       real :: reconvabsW1(nex),imconvabsW1(nex),reconvabsW1a(nex),imconvabsW1a(nex)

! !variables that can be usefull
! !      real :: sum,ReSx(nex),ImSx(nex),ReG(ne),ImG(ne),frac,
! !      integer :: nro_grid,nphi_grid,mesave,gesave,xesave,lrec,

! !time testing 
! !      real :: t0,t1,t2,t3
      
!       data firstcall /.true./
!       data needresp/.true./
!       data Cpsave/2/
!       save firstcall,Emax,Emin,dloge,earx
!       save lens,xbinhi,contx,needresp,me,ge,xe !,mesave,gesave,xesave
!       save paramsave,fhisave,flosave,nfsave,nro,nphi
!       save frobs,sdmin,sdmax,frrel,Cpsave !,hsave,rinsave
!       save transe,transea,transe_11,transe_12,transe_21,transe_22,transe_11a,transe_12a,transe_21a,transe_22a
!       save reconv,imconv,reconvW1,imconvW1,reconvW1a,imconvW1a,check,d,rnmax


!       pi = acos(-1.d0)
!       ifl = 1
      
! ! Settings
!       dlogf = 0.0073  !This is a resolution parameter (base 10)

! ! Call environment variables
!       verbose = myenv("REV_VERB",0)     !Set verbose level
      
! ! Initialise
! !      if( firstcall ) call FNINIT
!       call initialiser(firstcall,Emin,Emax,nex,dloge,earx,rnmax,d,needtrans,check&
!      ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim,me,ge,xe)


! !gridname = "/Users/Gullik/Dropbox/work/reflection_lag/Model_crossen/grid/grid_300_30x30_pem.dat"

! !Allocate dynamically the array to calculate the trasfer function          
!          if (.not. allocated(re1)) allocate(re1(nphi,nro))
!          if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
!          if (.not. allocated(pem1)) allocate(pem1(nphi,nro))

! ! Parameters
!       h        = dble( param(1) )
!       a        = dble( param(2) )
!       inc      = dble( param(3) )
!       rin      = dble( param(4) )
!       rout     = dble( param(5) )
!       zcos     = dble( param(6) )
!       Gamma    = dble( param(7) )
!       logxi    = param(8)
!       Afe      = param(9)
!       Ecut_obs = param(10)
!       Nh       = param(11)
!       afac     = param(12)
!       Mass     = dble( param(13) )
!       flo      = dble( param(14) )
!       fhi      = dble( param(15) )
!       ReIm     = int( param(16) )
!       phiA     = param(17)
!       phiB     = param(18)
!       g        = param(19)

!       honr = 0.d0
!       muobs = cos( inc * pi / 180.d0 )


! !check if the grid values are the same one of the model
!             if ( check .ne. 0 .and. honr_grid .ne. honr ) then
!             write(*,*) 'grid has a different honr!'
!             write(*,*) 'honr of the grid is ', honr_grid
!             stop
!          endif
      
! !Work out how many frequencies to average over
!       fc = 0.5d0 * ( flo + fhi )
!       nf = ceiling( log10(fhi/flo) / dlogf )
!       if( fhi .lt. tiny(fhi) .or. flo .lt. tiny(flo) )then
!         fhi = 0.d0
!         flo = 0.d0
!         nf  = 1
!       end if

! !Convert frequency bounds from Hz to c/Rg
!       fhi = fhi * 4.916d-6 * Mass
!       flo = flo * 4.916d-6 * Mass
      
! !Set minimum r (ISCO) and convert rin and h to rg
!       rmin   = disco( a )
!       if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
!       rh     = 1.d0+sqrt(1.d0-a**2)
!       if( h .lt. 0.d0 ) h = abs(h) * rh
!       if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
!       if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h
!       if( rin .lt. rmin )then
!         write(*,*)"Warning! rin<ISCO! Set to ISCO"
!         rin = rmin
!       end if

!       if( h .lt. 1.5d0*rh )then
!         write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
!         h = 1.5d0 * rh
!       end if

! !Calculate source to observer g-factor and source frame Ecut
!       gso    = real( dgsofac(a,h) )
!       Ecut_s = real(1.d0+zcos) * Ecut_obs / gso
!       if( verbose .gt. 0 )then
!         if( Cp .eq. 0 )then
!           write(*,*)"Ecut in source restframe (keV)=",Ecut_s
!         else
!           write(*,*)"kTe in source restframe (keV)=",Ecut_s
!         end if
!       end if
      
! !Determine if I need to calculate the kernel
!       if( .not. needtrans )then
!         do i = 1,8
!           if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
!         end do
!         if( abs( param(11) - paramsave(11) ) .gt. 1e-7 ) needtrans = .true.
!         if( nf .ne. nfsave ) needtrans = .true.
!         if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
!         if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
!       end if

!       ximin   = 0.0
!       ximax   = 4.7

! !      call CPU_TIME(t0)
      
!       if( needtrans )then
! !          if (check .ne. 0) then

! ! !Allocate dinamically the transfer functions
! !             if(.not. allocated(transe_11) ) allocate(transe_11(nex,me,ge,xe))
! !             if(.not. allocated(transe_12) ) allocate(transe_12(nex,me,ge,xe))
! !             if(.not. allocated(transe_21) ) allocate(transe_21(nex,me,ge,xe))
! !             if(.not. allocated(transe_22) ) allocate(transe_22(nex,me,ge,xe))
! !             if(.not. allocated(transe_11a) ) allocate(transe_11a(nex,me,ge,xe))
! !             if(.not. allocated(transe_12a) ) allocate(transe_12a(nex,me,ge,xe))
! !             if(.not. allocated(transe_21a) ) allocate(transe_21a(nex,me,ge,xe))
! !             if(.not. allocated(transe_22a) ) allocate(transe_22a(nex,me,ge,xe))
            
! ! !Choose the index of spin and inclination to extract from the grid            
! !             call chclose(a,spin_start,spin_end,spin_dim,s_lo,s_hi)
! !             call chclose(muobs,mu_start,mu_end,mu_dim,m_lo,m_hi)
! !             call ch_ind_val(spin_start,spin_end,spin_dim,s_lo,spin_lo)
! !             call ch_ind_val(spin_start,spin_end,spin_dim,s_hi,spin_hi)
! !             call ch_ind_val(mu_start,mu_end,mu_dim,m_lo,mu_lo)
! !             call ch_ind_val(mu_start,mu_end,mu_dim,m_hi,mu_hi)

! ! !Extraction from the binary grid according to the index and how they have been saved in the grid         
! !             irec = 3*mu_dim*(s_lo-1) + 3*m_lo - 1 
! !             read(98,rec=irec) re1
! !             read(98,rec=irec+1) taudo1
! !             read(98,rec=irec+2) pem1
! !             status_re_tau = .false.
! ! !Calculate the Kernel for the given parameters
! !             call strans(spin_lo,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
! !                  nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_ml,sdmax_sl_ml&
! !                  ,transe_11,transe_11a,frobs_sl_ml,frrel_sl_ml,xbinhi_sl_ml,lens)

! ! !Repeat this for all nthe combinations of spin and inclination (4 times)
! !             irec = 3*mu_dim*(s_lo-1) + 3*m_hi - 1
! !             read(98,rec=irec) re1
! !             read(98,rec=irec+1) taudo1
! !             read(98,rec=irec+2) pem1
! !             call strans(spin_lo,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
! !                  nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_mh,sdmax_sl_mh&
! !                  ,transe_21,transe_21a,frobs_sl_mh,frrel_sl_mh,xbinhi_sl_mh,lens)

! !             irec = 3*mu_dim*(s_hi-1) + 3*m_lo - 1
! !             read(98,rec=irec) re1
! !             read(98,rec=irec+1) taudo1
! !             read(98,rec=irec+2) pem1
! !             call strans(spin_hi,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
! !                  nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_ml,sdmax_sh_ml&
! !                  ,transe_12,transe_12a,frobs_sh_ml,frrel_sh_ml,xbinhi_sh_ml,lens)
            
! !             irec = 3*mu_dim*(s_hi-1) + 3*m_hi - 1
! !             read(98,rec=irec) re1
! !             read(98,rec=irec+1) taudo1
! !             read(98,rec=irec+2) pem1
! !             call strans(spin_hi,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
! !                  nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_mh,sdmax_sh_mh&
! !                  ,transe_22,transe_22a,frobs_sh_mh,frrel_sh_mh,xbinhi_sh_mh,lens)
   

! ! !Weighted average of the reflection fraction 
! !             frobs = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frobs_sl_ml,frobs_sl_mh,frobs_sh_ml,frobs_sh_mh)
! !             frrel = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frrel_sl_ml,frrel_sl_mh,frrel_sh_ml,frrel_sh_mh)
! ! !Determin the maximum value of the inonisation binning 
! !             xbinhi = max(xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh)
! !             sdmin = min(sdmin_sl_ml,sdmin_sl_mh,sdmin_sh_ml,sdmin_sh_mh)
! !             sdmax = max(sdmax_sl_ml,sdmax_sl_mh,sdmax_sh_ml,sdmax_sh_mh)
            
! !          else 

!          if( allocated(transe ) ) deallocate(transe )
!          if( allocated(transea) ) deallocate(transea)
!          allocate(transe(nf,nex,me,ge,xe))
!          allocate(transea(nf,nex,me,ge,xe))
! !Calculate the Kernel for the given parameters
!             status_re_tau = .true.
!             call strans(a,h,muobs,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
!                  nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin,sdmax,transe,transea,frobs,frrel,xbinhi,lens)
!          ! endif        
!       end if
      
!       if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
!       if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

! !Determine if I need to convolve with the restframe reflection spectrum
!       needconv = .false.
!       if( needtrans ) needconv = .true.
!       do i = 8,10
!         if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needconv = .true.
!       end do
!       if( Cp .ne. Cpsave ) needconv = .true.
      
!       if( needconv )then

! !Get continuum spectrum
!         call getcont(nex,earx,Gamma,Afe,Ecut_obs,logxi,Cp,contx,xillpar)
!         if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
! !Now reflection
!         xillpar(7) = -1.0       !reflection fraction of 1        
!         reconv = 0.0
!         imconv = 0.0
!         reconvW1 = 0.0
!         imconvW1 = 0.0
!         reconvW1a = 0.0
!         imconvW1a = 0.0
        
! !Loop over ionisatrion, emission angle and Ecut zones
!         do xbin = 1,xbinhi   !loop over ionisation zones
!           logxir = ximin + (xbin-0.5)*(ximax-ximin)/float(xe)
!           xillpar(4) = logxir
!           if( xe .eq. 1 ) xillpar(4) = logxi
!           do sdbin = 1,ge      !loop over Ecut zones
! !Calculate blueshift factor
!             gsd = sdmin + (sdbin-0.5) * (sdmax-sdmin)/float(ge)
!             !Set local high energy cut-off
!             xillpar(3) = gsd * Ecut_s
!             if( ge .eq. 1 ) xillpar(3) = Ecut_s
! !Loop through emission angles
!             do mubin = 1,me      !loop over emission angle zones
!               !Calculate input inclination angle
!               mue = ( real(mubin) - 0.5 ) / real(me)
!               xillpar(6) = acos( mue ) * 180.0 / real(pi)
!               if( me .eq. 1 ) xillpar(6) = real( inc )

! !Call xillver
!               xillpar(1) = real(Gamma)
!               call myxill(earx,nex,xillpar,ifl,Cp,photarx)

! !non linear effects
!                  DeltaGamma = 0.01
!                  Gamma1 = real(Gamma) - 0.5*DeltaGamma
!                  Gamma2 = real(Gamma) + 0.5*DeltaGamma

!                  xillpar(1) = Gamma1
!                  call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
!                  xillpar(1) = Gamma2
!                  call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
!                  photarx_delta = (photarx_2 - photarx_1)/DeltaGamma

                 
! !               if (check .ne. 0) then
! ! !interpolation of the calculated transfer functions 
! !                  call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_11(:,mubin,sdbin,xbin),transe_12(:,mubin,sdbin,xbin)&
! !                       ,transe_11a(:,mubin,sdbin,xbin),transe_12a(:,mubin,sdbin,xbin),transe_1,transe_1a)
! !                  call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_21(:,mubin,sdbin,xbin),transe_22(:,mubin,sdbin,xbin)&
! !                       ,transe_21a(:,mubin,sdbin,xbin),transe_22a(:,mubin,sdbin,xbin),transe_2,transe_2a)
! !                  call myinterp_complex2(muobs,mu_lo,mu_hi,nex,transe_1,transe_2,transe_1a,transe_2a,transe_,transe_a)

! !                  do i=1,nex
! !                     reline(i) = real(  transe_(i) ) 
! !                     imline(i) = aimag( transe_(i) )
! !                     reline_a(i) = real(  transe_a(i) ) 
! !                     imline_a(i) = aimag( transe_a(i) )
! !                  enddo

! !               else
                 
!                  do i = 1,nex
!                     reline(i) = real(  transe(i,mubin,sdbin,xbin) )
!                     imline(i) = aimag( transe(i,mubin,sdbin,xbin) )
!                     reline_a(i) = real(  transea(i,mubin,sdbin,xbin) )
!                     imline_a(i) = aimag( transea(i,mubin,sdbin,xbin) )
!                  end do
!               endif

! !Convolve with line profile
!               call padcnv(1e-7,nex,reline,photarx,reconvmu)
!               call padcnv(1e-7,nex,imline,photarx,imconvmu)
!               !Add to running sum
!               reconv = reconv + reconvmu
!               imconv = imconv + imconvmu

!               call padcnv(1e-7,nex,reline,photarx_delta,reconvmu)
!               call padcnv(1e-7,nex,imline,photarx_delta,imconvmu)
!               reconvW1 = reconvW1 + reconvmu
!               imconvW1 = imconvW1 + imconvmu

!               call padcnv(1e-7,nex,reline_a,photarx,reconvmu)
!               call padcnv(1e-7,nex,imline_a,photarx,imconvmu)
!               reconvW1a = reconvW1a + reconvmu
!               imconvW1a = imconvW1a + imconvmu
              
!             end do
!           end do
!         end do

!       end if

      
! ! Calculate absorption
!       call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
      
! ! Include absorption in the model
!       contxabs  = contx  * absorbx
!       reconvabs = reconv * absorbx
!       imconvabs = imconv * absorbx

!       reconvabsW1 = reconvW1 * absorbx
!       imconvabsW1 = imconvW1 * absorbx
!       reconvabsW1a = reconvW1a * absorbx
!       imconvabsW1a = imconvW1a * absorbx

!       ! do gbin = 1,nex
!       !   ghi = (1+zcos) * 10.0**( (gbin-nex/2)*dloge )
!       !   glo = (1+zcos) * 10.0**( (gbin-1-nex/2)*dloge )
!       !   write(410,*)6.4*0.5*(glo+ghi),real(transe(gbin,1,1,1))
!       ! end do
!       ! write(410,*)"no no"
      
! ! Calculate phiA from instrument response - if this option is set to on      
!       call phaseA(nex,earx,contxabs,reconvabs,imconvabs,gso,zcos,Gamma,afac,lens,phiA)
      
! !Add on continuum (and include boosting fudge factor) + consider the non-linear effects

!       !Stop unphysical parameters for the DC component
!       if( fhi .lt. tiny(fhi) .or. flo .lt. tiny(flo) )then
!         phiA = 0.0
!         g    = 0.0
!       end if

!       !Calculate cross-spectrum
!       do i = 1,nex
!         E  = 0.5 * ( earx(i) + earx(i-1) )
!         dE = earx(i) - earx(i-1)
!         !Direct energy-dependent component
! !        direct  = contxabs(i) / dE * real(lens) * ( gso / real(1.d0+zcos) )**real(2.d0+Gamma)   
!         direct  = contxabs(i) / dE * real(lens) * ( gso / real(1.d0+zcos) )**real(Gamma)   
!         !Factors
!         ldir = direct * log(E*(1+real(zcos))/gso)
!         ReW  = afac * reconvabs(i)/dE
!         ImW  = afac * imconvabs(i)/dE
!         ReW1 = afac * ( reconvabsW1(i) + reconvabsW1a(i) ) / dE
!         ImW1 = afac * ( imconvabsW1(i) + imconvabsW1a(i) ) / dE
!         !Real part
!         ReGx(i) = cos(phiA)*(direct+ReW) - sin(phiA)*ImW
!         ReGx(i) = ReGx(i) + g * ( cos(phiB)*(ldir-ReW1) + sin(phiB)*ImW1 )
!         !Imaginary part
!         ImGx(i) = cos(phiA)*ImW + sin(phiA)*(direct+ReW)
!         ImGx(i) = ImGx(i) + g * ( sin(phiB)*(ldir-ReW1) - cos(phiB)*ImW1 )        
!      end do

!         ! ReGx(i) = cos(phiA) * (direct + afac*reconvabs(i)/dE) - sin(phiA)*afac*imconvabs(i)/dE + &
!         ! g * (cos(phiB) * ( log(E*(1+real(zcos))/gso)*direct - afac*reconvabsW1a(i)/dE - afac*reconvabsW1(i)/dE)+&
!         ! sin(phiB)* (imconvabsW1a(i)/dE + imconvabsW1(i)/dE) )
!         ! ImGx(i) = sin(phiA) * (direct + afac*reconvabs(i)/dE) + cos(phiA)*afac*imconvabs(i)/dE + &
!         ! g * (sin(phiB) * ( log(E*(1+real(zcos))/gso)*direct - afac*reconvabsW1a(i)/dE - afac*reconvabsW1(i)/dE)-&
!         ! cos(phiB) * (afac*imconvabsW1a(i)/dE + afac*imconvabsW1(i)/dE) )
     
!       !Re-bin onto input grid
!       call rebinE(earx,ReGx,nex,ear,ReS,ne)
!       call rebinE(earx,ImGx,nex,ear,ImS,ne)
      
!       !Write output (depends on ReIm parameter)
!       if( ReIm .eq. 1 .or. flo .lt. 1e-10 )then   !Real part / DC part
!         do i = 1,ne
!           E = 0.5 * ( ear(i) + ear(i-1) )
!           dE = ear(i) - ear(i-1)
!           photar(i) = ReS(i) * dE
!         end do
!      else if( ReIm .eq. 2 )then   !Imaginary part
!         do i = 1,ne
!           E = 0.5 * ( ear(i) + ear(i-1) )
!           dE = ear(i) - ear(i-1)
!           photar(i) = ImS(i) * dE
!         end do
!       else if( ReIm .eq. 3 )then   !Modulus
!         do i = 1,ne
!           E = 0.5 * ( ear(i) + ear(i-1) )
!           dE = ear(i) - ear(i-1)
!           photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
!         end do
!         write(*,*)"Warning ReIm=3 should not be used for fitting!"
!      else if( ReIm .eq. 4 )then   !Time lag (seconds)
!         do i = 1,ne
!           E = 0.5 * ( ear(i) + ear(i-1) )
!           dE = ear(i) - ear(i-1)
!           phase = atan2( ImS(i) , ReS(i) )
!           photar(i) = phase / real(2.d0*pi*fc) * dE
!        end do
!        write(*,*)"Warning ReIm=4 should not be used for fitting!"

!     end if

!       !Write out options to fit for lags and amplitude
!       if( ReIm .gt. 4 )then
!         !Fold real and imaginary parts around the telescope response
!         call folder(nex,earx,ReGx,ImGx,ne,ear,ReS,ImS)
!         if( ReIm .eq. 5 )then   !Modulus
!           do i = 1,ne
!             E = 0.5 * ( ear(i) + ear(i-1) )
!             dE = ear(i) - ear(i-1)
!             photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
!           end do
!         else if( ReIm .eq. 6 )then   !Time lag (seconds)
!           do i = 1,ne
!             E = 0.5 * ( ear(i) + ear(i-1) )
!             dE = ear(i) - ear(i-1)
!             phase = atan2( ImS(i) , ReS(i) )
!             photar(i) = phase / real(2.d0*pi*fc) * dE
!           end do
!         end if
!       end if
      
!       fhisave   = fhi
!       flosave   = flo
!       nfsave    = nf
!       paramsave = param
!       Cpsave    = Cp
!       ! mesave    = me
!       ! gesave    = ge
!       ! xesave    = xe
      

!     end subroutine genreltrans
! !-----------------------------------------------------------------------



