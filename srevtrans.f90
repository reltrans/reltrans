! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for \ell * g^{2+Gamma}
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) Ionisation profile
! 5) The dependence on incident angle (mimicked by tweaking the ionization)

include 'subroutines/header.h'


!-----------------------------------------------------------------------
      subroutine tdreltrans(ear,ne,param,ifl,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(19),photar(ne)
      Cp = 0
      call genreltrans(Cp,ear,ne,param,ifl,photar)
      return
      end subroutine tdreltrans
!-----------------------------------------------------------------------
        

!-----------------------------------------------------------------------
      subroutine tdreltransCp(ear,ne,param,ifl,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(19),photar(ne)
      Cp = 1
      call genreltrans(Cp,ear,ne,param,ifl,photar)
      return
      end subroutine tdreltransCp
!-----------------------------------------------------------------------

      
!-----------------------------------------------------------------------
      subroutine tdreltransCpT(ear,ne,param,ifl,photar)
      implicit none
      integer ne,ifl,Cp
      real ear(0:ne),param(20),photar(ne)
      Cp = 1
      call genreltransT(Cp,ear,ne,param,ifl,photar)
      return
      end subroutine tdreltransCpT
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      subroutine genreltransT(Cp,ear,ne,param,ifl,photar)
      use dyn_gr
      implicit none
      integer nro,nphi,nex,i,nf,ifl,ne,ReIm,nfsave
      integer verbose,mubin,sdbin
      integer me,ge,xe,Cp
      parameter (nex=2**12)
      double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
      double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
      double precision fhisave,flosave,rh,frrel,lens,fc
      real afac,param(20),ear(0:ne),gso,direct,ximin,ximax
      real ldir,ReW,ImW,ReW1,ImW1
      real Afe,Ecut_s,Ecut_obs,logxi,xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
      real reline(nex),imline(nex),photarx(nex),reconv(nex),imconv(nex)
      real reconvmu(nex),imconvmu(nex),mue,sdmin,sdmax,gsd
      real phase,ReS(ne),ImS(ne),photar(ne)
      real paramsave(20),contx(nex),phiA,absorbx(nex),photerx(nex)
      real ReGx(nex),ImGx(nex),Nh
      real contxabs(nex),reconvabs(nex),imconvabs(nex)
      complex,dimension(:,:,:,:),allocatable :: transe,transea
      logical firstcall,needtrans,needconv,needresp
      integer xbin,xbinhi,myenv,Cpsave,gbin
      real logxir,ghi,glo,kTbb,kTbb_s,kTbb_r

!variable for the grid reading
      integer :: irec,spin_dim,mu_dim,check
      integer :: s_lo,s_hi,m_lo,m_hi,xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh
      double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      real :: sdmin_sl_ml,sdmax_sl_ml,sdmin_sl_mh,sdmax_sl_mh,sdmin_sh_ml,sdmax_sh_ml&
           ,sdmin_sh_mh,sdmax_sh_mh
      double precision :: frobs_sl_ml,frrel_sl_ml,frobs_sl_mh,frrel_sl_mh&
           ,frobs_sh_ml,frrel_sh_ml,frobs_sh_mh,frrel_sh_mh
      complex :: transe_1(nex),transe_2(nex),transe_(nex)
      complex,dimension(:,:,:,:),allocatable :: transe_11,transe_12,transe_21,transe_22
      
!variable for non linear effects
      complex :: transe_1a(nex),transe_2a(nex),transe_a(nex)
      complex,dimension(:,:,:,:),allocatable :: transe_11a,transe_12a,transe_21a,transe_22a
      real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,phiB,g
      real :: reline_a(nex),imline_a(nex),reconvW1(nex),imconvW1(nex),reconvW1a(nex),imconvW1a(nex)
      real :: reconvabsW1(nex),imconvabsW1(nex),reconvabsW1a(nex),imconvabsW1a(nex)

!variables that can be usefull
!      real :: sum,ReSx(nex),ImSx(nex),ReG(ne),ImG(ne),frac,
!      integer :: nro_grid,nphi_grid,mesave,gesave,xesave,lrec,

!time testing 
!      real :: t0,t1,t2,t3
      
      data firstcall /.true./
      data needresp/.true./
      data Cpsave/2/
      save firstcall,Emax,Emin,dloge,earx
      save lens,xbinhi,contx,needresp,me,ge,xe !,mesave,gesave,xesave
      save paramsave,fhisave,flosave,nfsave,nro,nphi
      save frobs,sdmin,sdmax,frrel,Cpsave !,hsave,rinsave
      save transe,transea,transe_11,transe_12,transe_21,transe_22,transe_11a,transe_12a,transe_21a,transe_22a
      save reconv,imconv,reconvW1,imconvW1,reconvW1a,imconvW1a,check,d,rnmax


      pi = acos(-1.d0)
      ifl = 1
      
! Settings
      dlogf = 0.0073  !This is a resolution parameter (base 10)

! Call environment variables
      verbose = myenv("REV_VERB",0)     !Set verbose level
      
! Initialise
      if( firstcall ) call FNINIT
      call initialiser(firstcall,Emin,Emax,nex,dloge,earx,rnmax,d,needtrans,check&
     ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim,me,ge,xe)


!gridname = "/Users/Gullik/Dropbox/work/reflection_lag/Model_crossen/grid/grid_300_30x30_pem.dat"

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
      kTbb     = param(10)
      Ecut_obs = param(11)
      Nh       = param(12)
      afac     = param(13)
      Mass     = dble( param(14) )
      flo      = dble( param(15) )
      fhi      = dble( param(16) )
      ReIm     = int( param(17) )
      phiA     = param(18)
      phiB     = param(19)
      g        = param(20)

      honr = 0.d0
      muobs = cos( inc * pi / 180.d0 )


!check if the grid values are the same one of the model
         if ( check .ne. 0 .and. honr_grid .ne. honr ) then
            write(*,*) 'grid has a different honr!'
            write(*,*) 'honr of the grid is ', honr_grid
            stop
         endif
      
!Work out how many frequencies to average over
      fc = 0.5d0 * ( flo + fhi )
      nf = ceiling( log10(fhi/flo) / dlogf )
      if( fhi .lt. 1d-10 .or. flo .lt. 1d-10 )then
        fhi = 0.d0
        flo = 0.d0
        nf  = 1
      end if

!Convert frequency bounds from Hz to c/Rg
      fhi = fhi * 4.916d-6 * Mass
      flo = flo * 4.916d-6 * Mass
      
!Set minimum r (ISCO) and convert rin and h to rg
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
      kTbb_s = real(1.d0+zcos) * kTbb / gso
      if( verbose .gt. 0 )then
         write(*,*)"kTe in source restframe (keV)=",Ecut_s
         write(*,*)"kTbb in source restframe (keV)=",kTbb_s
      end if
      
!Determine if I need to calculate the kernel
      if( .not. needtrans )then
        do i = 1,8
          if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
        end do
        if( nf .ne. nfsave ) needtrans = .true.
        if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
        if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
      end if
      ximin   = 0.0
      ximax   = 4.7
      
      if( needtrans )then
         if (check .ne. 0) then

!Allocate dynamically the transfer functions
            if(.not. allocated(transe_11) ) allocate(transe_11(nex,me,ge,xe))
            if(.not. allocated(transe_12) ) allocate(transe_12(nex,me,ge,xe))
            if(.not. allocated(transe_21) ) allocate(transe_21(nex,me,ge,xe))
            if(.not. allocated(transe_22) ) allocate(transe_22(nex,me,ge,xe))
            if(.not. allocated(transe_11a) ) allocate(transe_11a(nex,me,ge,xe))
            if(.not. allocated(transe_12a) ) allocate(transe_12a(nex,me,ge,xe))
            if(.not. allocated(transe_21a) ) allocate(transe_21a(nex,me,ge,xe))
            if(.not. allocated(transe_22a) ) allocate(transe_22a(nex,me,ge,xe))
            
!Choose the index of spin and inclination to extract from the grid            
            call chclose(a,spin_start,spin_end,spin_dim,s_lo,s_hi)
            call chclose(muobs,mu_start,mu_end,mu_dim,m_lo,m_hi)
            call ch_ind_val(spin_start,spin_end,spin_dim,s_lo,spin_lo)
            call ch_ind_val(spin_start,spin_end,spin_dim,s_hi,spin_hi)
            call ch_ind_val(mu_start,mu_end,mu_dim,m_lo,mu_lo)
            call ch_ind_val(mu_start,mu_end,mu_dim,m_hi,mu_hi)

!Extraction from the binary grid according to the index and how they have been saved in the grid         
            irec = 3*mu_dim*(s_lo-1) + 3*m_lo - 1 
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            status_re_tau = .false.
!Calculate the Kernel for the given parameters
            call strans(spin_lo,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_ml,sdmax_sl_ml&
                 ,transe_11,transe_11a,frobs_sl_ml,frrel_sl_ml,xbinhi_sl_ml,lens)

!Repeat this for all nthe combinations of spin and inclination (4 times)
            irec = 3*mu_dim*(s_lo-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_lo,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_mh,sdmax_sl_mh&
                 ,transe_21,transe_21a,frobs_sl_mh,frrel_sl_mh,xbinhi_sl_mh,lens)

            irec = 3*mu_dim*(s_hi-1) + 3*m_lo - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_ml,sdmax_sh_ml&
                 ,transe_12,transe_12a,frobs_sh_ml,frrel_sh_ml,xbinhi_sh_ml,lens)
            
            irec = 3*mu_dim*(s_hi-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_mh,sdmax_sh_mh&
                 ,transe_22,transe_22a,frobs_sh_mh,frrel_sh_mh,xbinhi_sh_mh,lens)
   

!Weighted average of the reflection fraction 
            frobs = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frobs_sl_ml,frobs_sl_mh,frobs_sh_ml,frobs_sh_mh)
            frrel = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frrel_sl_ml,frrel_sl_mh,frrel_sh_ml,frrel_sh_mh)
!Determin the maximum value of the inonisation binning 
            xbinhi = max(xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh)
            sdmin = min(sdmin_sl_ml,sdmin_sl_mh,sdmin_sh_ml,sdmin_sh_mh)
            sdmax = max(sdmax_sl_ml,sdmax_sl_mh,sdmax_sh_ml,sdmax_sh_mh)
            
         else 

            if(.not. allocated(transe) ) allocate(transe(nex,me,ge,xe))
            if(.not. allocated(transea) ) allocate(transea(nex,me,ge,xe))

!Calculate the Kernel for the given parameters
            status_re_tau = .true.
            call strans(a,h,muobs,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin,sdmax,transe,transea,frobs,frrel,xbinhi,lens)
         endif         
      end if
      
      if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
      if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel
      
!Determine if I need to convolve with the restframe reflection spectrum
      needconv = .false.
      if( needtrans ) needconv = .true.
      do i = 8,12
        if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needconv = .true.
      end do
      if( Cp .ne. Cpsave ) needconv = .true.
      
      if( needconv )then

!Get continuum spectrum
        call getcont_kTbb(nex,earx,Gamma,Afe,kTbb,Ecut_obs,logxi,contx,xillpar)
        if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
        
!Now reflection
        xillpar(7) = -1.0       !reflection fraction of 1        
        reconv = 0.0
        imconv = 0.0
        reconvW1 = 0.0
        imconvW1 = 0.0
        reconvW1a = 0.0
        imconvW1a = 0.0
        
!Loop over ionisatrion, emission angle and Ecut zones
        do xbin = 1,xbinhi   !loop over ionisation zones
          logxir = ximin + (xbin-0.5)*(ximax-ximin)/float(xe)
          xillpar(4) = logxir
          if( xe .eq. 1 ) xillpar(4) = logxi
          do sdbin = 1,ge      !loop over Ecut zones
!Calculate blueshift factor
            gsd = sdmin + (sdbin-0.5) * (sdmax-sdmin)/float(ge)
            !Set local high energy cut-off
            xillpar(3) = gsd * Ecut_s
            kTbb_r     = gsd * kTbb_s
            if( ge .eq. 1 )then
               xillpar(3) = Ecut_s
               kTbb_r     = kTbb_s
            end if
!Loop through emission angles
            do mubin = 1,me      !loop over emission angle zones
              !Calculate input inclination angle
              mue = ( real(mubin) - 0.5 ) / real(me)
              xillpar(6) = acos( mue ) * 180.0 / real(pi)
              if( me .eq. 1 ) xillpar(6) = real( inc )

!Call xillver
              xillpar(1) = Gamma
!              call myxill(earx,nex,xillpar,ifl,Cp,photarx)
              call myxill_T(earx,nex,xillpar,kTbb_r,ifl,photarx)
              
!non linear effects
                 DeltaGamma = 0.01
                 Gamma1 = param(7) -0.5*DeltaGamma
                 Gamma2 = param(7) +0.5*DeltaGamma

                 xillpar(1) = Gamma1
                 call myxill_T(earx,nex,xillpar,kTbb_r,ifl,photarx_1)
                 
                 xillpar(1) = Gamma2
                 call myxill_T(earx,nex,xillpar,kTbb_r,ifl,photarx_2)
                 photarx_delta = (photarx_2 - photarx_1)/DeltaGamma

                 
              if (check .ne. 0) then
!interpolation of the calculated transfer functions 
                 call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_11(:,mubin,sdbin,xbin),transe_12(:,mubin,sdbin,xbin)&
                      ,transe_11a(:,mubin,sdbin,xbin),transe_12a(:,mubin,sdbin,xbin),transe_1,transe_1a)
                 call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_21(:,mubin,sdbin,xbin),transe_22(:,mubin,sdbin,xbin)&
                      ,transe_21a(:,mubin,sdbin,xbin),transe_22a(:,mubin,sdbin,xbin),transe_2,transe_2a)
                 call myinterp_complex2(muobs,mu_lo,mu_hi,nex,transe_1,transe_2,transe_1a,transe_2a,transe_,transe_a)

                 do i=1,nex
                    reline(i) = real(  transe_(i) ) 
                    imline(i) = aimag( transe_(i) )
                    reline_a(i) = real(  transe_a(i) ) 
                    imline_a(i) = aimag( transe_a(i) )
                 enddo

              else
                 
                 do i = 1,nex
                    reline(i) = real(  transe(i,mubin,sdbin,xbin) )
                    imline(i) = aimag( transe(i,mubin,sdbin,xbin) )
                    reline_a(i) = real(  transea(i,mubin,sdbin,xbin) )
                    imline_a(i) = aimag( transea(i,mubin,sdbin,xbin) )
                 end do
              endif

!Convolve with line profile
              call padcnv(1e-7,nex,reline,photarx,reconvmu)
              call padcnv(1e-7,nex,imline,photarx,imconvmu)
              !Add to running sum
              reconv = reconv + reconvmu
              imconv = imconv + imconvmu

              call padcnv(1e-7,nex,reline,photarx_delta,reconvmu)
              call padcnv(1e-7,nex,imline,photarx_delta,imconvmu)
              reconvW1 = reconvW1 + reconvmu
              imconvW1 = imconvW1 + imconvmu

              call padcnv(1e-7,nex,reline_a,photarx,reconvmu)
              call padcnv(1e-7,nex,imline_a,photarx,imconvmu)
              reconvW1a = reconvW1a + reconvmu
              imconvW1a = imconvW1a + imconvmu
              
            end do
          end do
        end do

      end if

      
! Calculate absorption
      call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
      
! Include absorption in the model
      contxabs  = contx  * absorbx
      reconvabs = reconv * absorbx
      imconvabs = imconv * absorbx

      reconvabsW1 = reconvW1 * absorbx
      imconvabsW1 = imconvW1 * absorbx
      reconvabsW1a = reconvW1a * absorbx
      imconvabsW1a = imconvW1a * absorbx

      ! do gbin = 1,nex
      !   ghi = (1+zcos) * 10.0**( (gbin-nex/2)*dloge )
      !   glo = (1+zcos) * 10.0**( (gbin-1-nex/2)*dloge )
      !   write(410,*)6.4*0.5*(glo+ghi),real(transe(gbin,1,1,1))
      ! end do
      ! write(410,*)"no no"
      
! Calculate phiA from instrument response - if this option is set to on      
      call phaseA(nex,earx,contxabs,reconvabs,imconvabs,gso,zcos,Gamma,afac,lens,phiA)
      
!Add on continuum (and include boosting fudge factor) + consider the non-linear effects

      !Stop unphysical parameters for the DC component
      if( fhi .lt. tiny(fhi) .or. flo .lt. tiny(flo) )then
        phiA = 0.0
        g    = 0.0
      end if

      !Calculate cross-spectrum
      do i = 1,nex
        E  = 0.5 * ( earx(i) + earx(i-1) )
        dE = earx(i) - earx(i-1)
        !Direct energy-dependent component
        direct  = contxabs(i) / dE * real(lens) * ( gso / real(1.d0+zcos) )**real(Gamma)
        !Factors
        ldir = direct * log(E*(1+real(zcos))/gso)
        ReW  = afac * reconvabs(i)/dE
        ImW  = afac * imconvabs(i)/dE
        ReW1 = afac * ( reconvabsW1(i) + reconvabsW1a(i) ) / dE
        ImW1 = afac * ( imconvabsW1(i) + imconvabsW1a(i) ) / dE
        !write(80,*)E,E**2*direct,E**2*ReW
        !Real part
        ReGx(i) = cos(phiA)*(direct+ReW) - sin(phiA)*ImW
        ReGx(i) = ReGx(i) + g * ( cos(phiB)*(ldir-ReW1) + sin(phiB)*ImW1 )
        !Imaginary part
        ImGx(i) = cos(phiA)*ImW + sin(phiA)*(direct+ReW)
        ImGx(i) = ImGx(i) + g * ( sin(phiB)*(ldir-ReW1) - cos(phiB)*ImW1 )        
     end do
     
      !Re-bin onto input grid
      call rebinE(earx,ReGx,nex,ear,ReS,ne)
      call rebinE(earx,ImGx,nex,ear,ImS,ne)
      
      !Write output (depends on ReIm parameter)
      if( ReIm .eq. 1 .or. flo .lt. 1e-10 )then   !Real part / DC part
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = ReS(i) * dE
        end do
     else if( ReIm .eq. 2 )then   !Imaginary part
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = ImS(i) * dE
        end do
      else if( ReIm .eq. 3 )then   !Modulus
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
        end do
        write(*,*)"Warning ReIm=3 should not be used for fitting!"
     else if( ReIm .eq. 4 )then   !Time lag (seconds)
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          phase = atan2( ImS(i) , ReS(i) )
          photar(i) = phase / real(2.d0*pi*fc) * dE
       end do
       write(*,*)"Warning ReIm=4 should not be used for fitting!"

    end if

      !Write out options to fit for lags and amplitude
      if( ReIm .gt. 4 )then
        !Fold real and imaginary parts around the telescope response
        call folder(nex,earx,ReGx,ImGx,ne,ear,ReS,ImS)
        if( ReIm .eq. 5 )then   !Modulus
          do i = 1,ne
            E = 0.5 * ( ear(i) + ear(i-1) )
            dE = ear(i) - ear(i-1)
            photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
          end do
        else if( ReIm .eq. 6 )then   !Time lag (seconds)
          do i = 1,ne
            E = 0.5 * ( ear(i) + ear(i-1) )
            dE = ear(i) - ear(i-1)
            phase = atan2( ImS(i) , ReS(i) )
            photar(i) = phase / real(2.d0*pi*fc) * dE
          end do
        end if
      end if
      
      fhisave   = fhi
      flosave   = flo
      nfsave    = nf
      paramsave = param
      Cpsave    = Cp
      ! mesave    = me
      ! gesave    = ge
      ! xesave    = xe
      

    end subroutine genreltransT
!-----------------------------------------------------------------------


      

!-----------------------------------------------------------------------
      subroutine genreltrans(Cp,ear,ne,param,ifl,photar)
      use dyn_gr
      implicit none
      integer nro,nphi,nex,i,nf,ifl,ne,ReIm,nfsave
      integer verbose,mubin,sdbin
      integer me,ge,xe,Cp
      parameter (nex=2**12)
      double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
      double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
      double precision fhisave,flosave,rh,frrel,lens,fc
      real afac,param(19),ear(0:ne),gso,direct,ximin,ximax
      real ldir,ReW,ImW,ReW1,ImW1
      real Afe,Ecut_s,Ecut_obs,logxi,xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
      real reline(nex),imline(nex),photarx(nex),reconv(nex),imconv(nex)
      real reconvmu(nex),imconvmu(nex),mue,sdmin,sdmax,gsd
      real phase,ReS(ne),ImS(ne),photar(ne)
      real paramsave(19),contx(nex),phiA,absorbx(nex),photerx(nex)
      real ReGx(nex),ImGx(nex),Nh
      real contxabs(nex),reconvabs(nex),imconvabs(nex)
      complex,dimension(:,:,:,:),allocatable :: transe,transea
      logical firstcall,needtrans,needconv,needresp
      integer xbin,xbinhi,myenv,Cpsave,gbin
      real logxir,ghi,glo

!variable for the grid reading
      integer :: irec,spin_dim,mu_dim,check
      integer :: s_lo,s_hi,m_lo,m_hi,xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh
      double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      real :: sdmin_sl_ml,sdmax_sl_ml,sdmin_sl_mh,sdmax_sl_mh,sdmin_sh_ml,sdmax_sh_ml&
           ,sdmin_sh_mh,sdmax_sh_mh
      double precision :: frobs_sl_ml,frrel_sl_ml,frobs_sl_mh,frrel_sl_mh&
           ,frobs_sh_ml,frrel_sh_ml,frobs_sh_mh,frrel_sh_mh
      complex :: transe_1(nex),transe_2(nex),transe_(nex)
      complex,dimension(:,:,:,:),allocatable :: transe_11,transe_12,transe_21,transe_22
      
!variable for non linear effects
      complex :: transe_1a(nex),transe_2a(nex),transe_a(nex)
      complex,dimension(:,:,:,:),allocatable :: transe_11a,transe_12a,transe_21a,transe_22a
      real :: photarx_1(nex),photarx_2(nex),photarx_delta(nex),Gamma1,Gamma2,DeltaGamma,phiB,g
      real :: reline_a(nex),imline_a(nex),reconvW1(nex),imconvW1(nex),reconvW1a(nex),imconvW1a(nex)
      real :: reconvabsW1(nex),imconvabsW1(nex),reconvabsW1a(nex),imconvabsW1a(nex)

!variables that can be usefull
!      real :: sum,ReSx(nex),ImSx(nex),ReG(ne),ImG(ne),frac,
!      integer :: nro_grid,nphi_grid,mesave,gesave,xesave,lrec,

!time testing 
!      real :: t0,t1,t2,t3
      
      data firstcall /.true./
      data needresp/.true./
      data Cpsave/2/
      save firstcall,Emax,Emin,dloge,earx
      save lens,xbinhi,contx,needresp,me,ge,xe !,mesave,gesave,xesave
      save paramsave,fhisave,flosave,nfsave,nro,nphi
      save frobs,sdmin,sdmax,frrel,Cpsave !,hsave,rinsave
      save transe,transea,transe_11,transe_12,transe_21,transe_22,transe_11a,transe_12a,transe_21a,transe_22a
      save reconv,imconv,reconvW1,imconvW1,reconvW1a,imconvW1a,check,d,rnmax


      pi = acos(-1.d0)
      ifl = 1
      
! Settings
      dlogf = 0.0073  !This is a resolution parameter (base 10)

! Call environment variables
      verbose = myenv("REV_VERB",0)     !Set verbose level
      
! Initialise
      if( firstcall ) call FNINIT
      call initialiser(firstcall,Emin,Emax,nex,dloge,earx,rnmax,d,needtrans,check&
     ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim,me,ge,xe)


!gridname = "/Users/Gullik/Dropbox/work/reflection_lag/Model_crossen/grid/grid_300_30x30_pem.dat"

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
      Ecut_obs = param(10)
      Nh       = param(11)
      afac     = param(12)
      Mass     = dble( param(13) )
      flo      = dble( param(14) )
      fhi      = dble( param(15) )
      ReIm     = int( param(16) )
      phiA     = param(17)
      phiB     = param(18)
      g        = param(19)

      honr = 0.d0
      muobs = cos( inc * pi / 180.d0 )


!check if the grid values are the same one of the model
            if ( check .ne. 0 .and. honr_grid .ne. honr ) then
            write(*,*) 'grid has a different honr!'
            write(*,*) 'honr of the grid is ', honr_grid
            stop
         endif
      
!Work out how many frequencies to average over
      fc = 0.5d0 * ( flo + fhi )
      nf = ceiling( log10(fhi/flo) / dlogf )
      if( fhi .lt. 1d-10 .or. flo .lt. 1d-10 )then
        fhi = 0.d0
        flo = 0.d0
        nf  = 1
      end if

!Convert frequency bounds from Hz to c/Rg
      fhi = fhi * 4.916d-6 * Mass
      flo = flo * 4.916d-6 * Mass
      
!Set minimum r (ISCO) and convert rin and h to rg
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
      if( verbose .gt. 0 )then
        if( Cp .eq. 0 )then
          write(*,*)"Ecut in source restframe (keV)=",Ecut_s
        else
          write(*,*)"kTe in source restframe (keV)=",Ecut_s
        end if
      end if
      
!Determine if I need to calculate the kernel
      if( .not. needtrans )then
        do i = 1,8
          if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needtrans = .true.
        end do
        if( abs( param(11) - paramsave(11) ) .gt. 1e-7 ) needtrans = .true.
        if( nf .ne. nfsave ) needtrans = .true.
        if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
        if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
        ! if( abs( h - hsave ) .gt. 1e-7 ) needtrans = .true.
        ! if( abs( rin - rinsave ) .gt. 1e-7 ) needtrans = .true.
        ! if( nro .ne. nrosave ) needtrans = .true.
        ! if( nphi .ne. nphisave ) needtrans = .true.
        ! if( me .ne. mesave ) needtrans = .true.
        ! if( ge .ne. gesave ) needtrans = .true.
        ! if( xe .ne. xesave ) needtrans = .true.
      end if

      ximin   = 0.0
      ximax   = 4.7

!      call CPU_TIME(t0)
      
      if( needtrans )then
         if (check .ne. 0) then

!Allocate dinamically the transfer functions
            if(.not. allocated(transe_11) ) allocate(transe_11(nex,me,ge,xe))
            if(.not. allocated(transe_12) ) allocate(transe_12(nex,me,ge,xe))
            if(.not. allocated(transe_21) ) allocate(transe_21(nex,me,ge,xe))
            if(.not. allocated(transe_22) ) allocate(transe_22(nex,me,ge,xe))
            if(.not. allocated(transe_11a) ) allocate(transe_11a(nex,me,ge,xe))
            if(.not. allocated(transe_12a) ) allocate(transe_12a(nex,me,ge,xe))
            if(.not. allocated(transe_21a) ) allocate(transe_21a(nex,me,ge,xe))
            if(.not. allocated(transe_22a) ) allocate(transe_22a(nex,me,ge,xe))
            
!Choose the index of spin and inclination to extract from the grid            
            call chclose(a,spin_start,spin_end,spin_dim,s_lo,s_hi)
            call chclose(muobs,mu_start,mu_end,mu_dim,m_lo,m_hi)
            call ch_ind_val(spin_start,spin_end,spin_dim,s_lo,spin_lo)
            call ch_ind_val(spin_start,spin_end,spin_dim,s_hi,spin_hi)
            call ch_ind_val(mu_start,mu_end,mu_dim,m_lo,mu_lo)
            call ch_ind_val(mu_start,mu_end,mu_dim,m_hi,mu_hi)

!Extraction from the binary grid according to the index and how they have been saved in the grid         
            irec = 3*mu_dim*(s_lo-1) + 3*m_lo - 1 
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            status_re_tau = .false.
!Calculate the Kernel for the given parameters
            call strans(spin_lo,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_ml,sdmax_sl_ml&
                 ,transe_11,transe_11a,frobs_sl_ml,frrel_sl_ml,xbinhi_sl_ml,lens)

!Repeat this for all nthe combinations of spin and inclination (4 times)
            irec = 3*mu_dim*(s_lo-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_lo,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sl_mh,sdmax_sl_mh&
                 ,transe_21,transe_21a,frobs_sl_mh,frrel_sl_mh,xbinhi_sl_mh,lens)

            irec = 3*mu_dim*(s_hi-1) + 3*m_lo - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_ml,sdmax_sh_ml&
                 ,transe_12,transe_12a,frobs_sh_ml,frrel_sh_ml,xbinhi_sh_ml,lens)
            
            irec = 3*mu_dim*(s_hi-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin_sh_mh,sdmax_sh_mh&
                 ,transe_22,transe_22a,frobs_sh_mh,frrel_sh_mh,xbinhi_sh_mh,lens)
   

!Weighted average of the reflection fraction 
            frobs = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frobs_sl_ml,frobs_sl_mh,frobs_sh_ml,frobs_sh_mh)
            frrel = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frrel_sl_ml,frrel_sl_mh,frrel_sh_ml,frrel_sh_mh)
!Determin the maximum value of the inonisation binning 
            xbinhi = max(xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh)
            sdmin = min(sdmin_sl_ml,sdmin_sl_mh,sdmin_sh_ml,sdmin_sh_mh)
            sdmax = max(sdmax_sl_ml,sdmax_sl_mh,sdmax_sh_ml,sdmax_sh_mh)
            
         else 

            if(.not. allocated(transe) ) allocate(transe(nex,me,ge,xe))
            if(.not. allocated(transea) ) allocate(transea(nex,me,ge,xe))

!Calculate the Kernel for the given parameters
            status_re_tau = .true.
            call strans(a,h,muobs,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,nex,dloge,&
                 nf,fhi,flo,ximin,ximax,me,ge,xe,logxi,sdmin,sdmax,transe,transea,frobs,frrel,xbinhi,lens)
         endif         
      end if
      
      if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
      if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel

!Determine if I need to convolve with the restframe reflection spectrum
      needconv = .false.
      if( needtrans ) needconv = .true.
      do i = 8,10
        if( abs( param(i) - paramsave(i) ) .gt. 1e-7 ) needconv = .true.
      end do
      if( Cp .ne. Cpsave ) needconv = .true.
      
      if( needconv )then

!Get continuum spectrum
        call getcont(nex,earx,Gamma,Afe,Ecut_obs,logxi,Cp,contx,xillpar)
        if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass),gso,real(Gamma))
!Now reflection
        xillpar(7) = -1.0       !reflection fraction of 1        
        reconv = 0.0
        imconv = 0.0
        reconvW1 = 0.0
        imconvW1 = 0.0
        reconvW1a = 0.0
        imconvW1a = 0.0
        
!Loop over ionisatrion, emission angle and Ecut zones
        do xbin = 1,xbinhi   !loop over ionisation zones
          logxir = ximin + (xbin-0.5)*(ximax-ximin)/float(xe)
          xillpar(4) = logxir
          if( xe .eq. 1 ) xillpar(4) = logxi
          do sdbin = 1,ge      !loop over Ecut zones
!Calculate blueshift factor
            gsd = sdmin + (sdbin-0.5) * (sdmax-sdmin)/float(ge)
            !Set local high energy cut-off
            xillpar(3) = gsd * Ecut_s
            if( ge .eq. 1 ) xillpar(3) = Ecut_s
!Loop through emission angles
            do mubin = 1,me      !loop over emission angle zones
              !Calculate input inclination angle
              mue = ( real(mubin) - 0.5 ) / real(me)
              xillpar(6) = acos( mue ) * 180.0 / real(pi)
              if( me .eq. 1 ) xillpar(6) = real( inc )

!Call xillver
              xillpar(1) = Gamma
              call myxill(earx,nex,xillpar,ifl,Cp,photarx)

!non linear effects
                 DeltaGamma = 0.01
                 Gamma1 = param(7) -0.5*DeltaGamma
                 Gamma2 = param(7) +0.5*DeltaGamma

                 xillpar(1) = Gamma1
                 call myxill(earx,nex,xillpar,ifl,Cp,photarx_1)
                 xillpar(1) = Gamma2
                 call myxill(earx,nex,xillpar,ifl,Cp,photarx_2)
                 photarx_delta = (photarx_2 - photarx_1)/DeltaGamma

                 
              if (check .ne. 0) then
!interpolation of the calculated transfer functions 
                 call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_11(:,mubin,sdbin,xbin),transe_12(:,mubin,sdbin,xbin)&
                      ,transe_11a(:,mubin,sdbin,xbin),transe_12a(:,mubin,sdbin,xbin),transe_1,transe_1a)
                 call myinterp_complex2(a,spin_lo,spin_hi,nex,transe_21(:,mubin,sdbin,xbin),transe_22(:,mubin,sdbin,xbin)&
                      ,transe_21a(:,mubin,sdbin,xbin),transe_22a(:,mubin,sdbin,xbin),transe_2,transe_2a)
                 call myinterp_complex2(muobs,mu_lo,mu_hi,nex,transe_1,transe_2,transe_1a,transe_2a,transe_,transe_a)

                 do i=1,nex
                    reline(i) = real(  transe_(i) ) 
                    imline(i) = aimag( transe_(i) )
                    reline_a(i) = real(  transe_a(i) ) 
                    imline_a(i) = aimag( transe_a(i) )
                 enddo

              else
                 
                 do i = 1,nex
                    reline(i) = real(  transe(i,mubin,sdbin,xbin) )
                    imline(i) = aimag( transe(i,mubin,sdbin,xbin) )
                    reline_a(i) = real(  transea(i,mubin,sdbin,xbin) )
                    imline_a(i) = aimag( transea(i,mubin,sdbin,xbin) )
                 end do
              endif

!Convolve with line profile
              call padcnv(1e-7,nex,reline,photarx,reconvmu)
              call padcnv(1e-7,nex,imline,photarx,imconvmu)
              !Add to running sum
              reconv = reconv + reconvmu
              imconv = imconv + imconvmu

              call padcnv(1e-7,nex,reline,photarx_delta,reconvmu)
              call padcnv(1e-7,nex,imline,photarx_delta,imconvmu)
              reconvW1 = reconvW1 + reconvmu
              imconvW1 = imconvW1 + imconvmu

              call padcnv(1e-7,nex,reline_a,photarx,reconvmu)
              call padcnv(1e-7,nex,imline_a,photarx,imconvmu)
              reconvW1a = reconvW1a + reconvmu
              imconvW1a = imconvW1a + imconvmu
              
            end do
          end do
        end do

      end if

      
! Calculate absorption
      call tbabs(earx,nex,nh,Ifl,absorbx,photerx)
      
! Include absorption in the model
      contxabs  = contx  * absorbx
      reconvabs = reconv * absorbx
      imconvabs = imconv * absorbx

      reconvabsW1 = reconvW1 * absorbx
      imconvabsW1 = imconvW1 * absorbx
      reconvabsW1a = reconvW1a * absorbx
      imconvabsW1a = imconvW1a * absorbx

      ! do gbin = 1,nex
      !   ghi = (1+zcos) * 10.0**( (gbin-nex/2)*dloge )
      !   glo = (1+zcos) * 10.0**( (gbin-1-nex/2)*dloge )
      !   write(410,*)6.4*0.5*(glo+ghi),real(transe(gbin,1,1,1))
      ! end do
      ! write(410,*)"no no"
      
! Calculate phiA from instrument response - if this option is set to on      
      call phaseA(nex,earx,contxabs,reconvabs,imconvabs,gso,zcos,Gamma,afac,lens,phiA)
      
!Add on continuum (and include boosting fudge factor) + consider the non-linear effects

      !Stop unphysical parameters for the DC component
      if( fhi .lt. tiny(fhi) .or. flo .lt. tiny(flo) )then
        phiA = 0.0
        g    = 0.0
      end if

      !Calculate cross-spectrum
      do i = 1,nex
        E  = 0.5 * ( earx(i) + earx(i-1) )
        dE = earx(i) - earx(i-1)
        !Direct energy-dependent component
!        direct  = contxabs(i) / dE * real(lens) * ( gso / real(1.d0+zcos) )**real(2.d0+Gamma)   
        direct  = contxabs(i) / dE * real(lens) * ( gso / real(1.d0+zcos) )**real(Gamma)   
        !Factors
        ldir = direct * log(E*(1+real(zcos))/gso)
        ReW  = afac * reconvabs(i)/dE
        ImW  = afac * imconvabs(i)/dE
        ReW1 = afac * ( reconvabsW1(i) + reconvabsW1a(i) ) / dE
        ImW1 = afac * ( imconvabsW1(i) + imconvabsW1a(i) ) / dE
        !Real part
        ReGx(i) = cos(phiA)*(direct+ReW) - sin(phiA)*ImW
        ReGx(i) = ReGx(i) + g * ( cos(phiB)*(ldir-ReW1) + sin(phiB)*ImW1 )
        !Imaginary part
        ImGx(i) = cos(phiA)*ImW + sin(phiA)*(direct+ReW)
        ImGx(i) = ImGx(i) + g * ( sin(phiB)*(ldir-ReW1) - cos(phiB)*ImW1 )        
     end do

        ! ReGx(i) = cos(phiA) * (direct + afac*reconvabs(i)/dE) - sin(phiA)*afac*imconvabs(i)/dE + &
        ! g * (cos(phiB) * ( log(E*(1+real(zcos))/gso)*direct - afac*reconvabsW1a(i)/dE - afac*reconvabsW1(i)/dE)+&
        ! sin(phiB)* (imconvabsW1a(i)/dE + imconvabsW1(i)/dE) )
        ! ImGx(i) = sin(phiA) * (direct + afac*reconvabs(i)/dE) + cos(phiA)*afac*imconvabs(i)/dE + &
        ! g * (sin(phiB) * ( log(E*(1+real(zcos))/gso)*direct - afac*reconvabsW1a(i)/dE - afac*reconvabsW1(i)/dE)-&
        ! cos(phiB) * (afac*imconvabsW1a(i)/dE + afac*imconvabsW1(i)/dE) )
     
      !Re-bin onto input grid
      call rebinE(earx,ReGx,nex,ear,ReS,ne)
      call rebinE(earx,ImGx,nex,ear,ImS,ne)
      
      !Write output (depends on ReIm parameter)
      if( ReIm .eq. 1 .or. flo .lt. 1e-10 )then   !Real part / DC part
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = ReS(i) * dE
        end do
     else if( ReIm .eq. 2 )then   !Imaginary part
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = ImS(i) * dE
        end do
      else if( ReIm .eq. 3 )then   !Modulus
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
        end do
        write(*,*)"Warning ReIm=3 should not be used for fitting!"
     else if( ReIm .eq. 4 )then   !Time lag (seconds)
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          phase = atan2( ImS(i) , ReS(i) )
          photar(i) = phase / real(2.d0*pi*fc) * dE
       end do
       write(*,*)"Warning ReIm=4 should not be used for fitting!"

    end if

      !Write out options to fit for lags and amplitude
      if( ReIm .gt. 4 )then
        !Fold real and imaginary parts around the telescope response
        call folder(nex,earx,ReGx,ImGx,ne,ear,ReS,ImS)
        if( ReIm .eq. 5 )then   !Modulus
          do i = 1,ne
            E = 0.5 * ( ear(i) + ear(i-1) )
            dE = ear(i) - ear(i-1)
            photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) * dE
          end do
        else if( ReIm .eq. 6 )then   !Time lag (seconds)
          do i = 1,ne
            E = 0.5 * ( ear(i) + ear(i-1) )
            dE = ear(i) - ear(i-1)
            phase = atan2( ImS(i) , ReS(i) )
            photar(i) = phase / real(2.d0*pi*fc) * dE
          end do
        end if
      end if
      
      fhisave   = fhi
      flosave   = flo
      nfsave    = nf
      paramsave = param
      Cpsave    = Cp
      ! mesave    = me
      ! gesave    = ge
      ! xesave    = xe
      

    end subroutine genreltrans
!-----------------------------------------------------------------------



