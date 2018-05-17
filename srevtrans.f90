! gfortran -L/usr/local/opt/cfitsio/ -lcfitsio srevtrans.f90 /Users/adamingram/Dropbox/relxill/relbase.f90 /Users/adamingram/Dropbox/relxill/relxill.f90
! gfortran -L/usr/local/opt/cfitsio/ -lcfitsio srevtrans.f90 relbase.f90 relxill.f90
!
! gfortran -L/Data/adamingram/heasoft/heasoft-6.22.1/heacore/cfitsio/ -lcfitsio srevtrans.f90 relbase.f90 relxill.f90
!
! make -f revmakefile run   or make -f desktop_revmakefile run
! ./exe
!
! To do:
! x Bring includes into main body of code
! x Create lmodel.dat & make into xspec code
! x Do firstcalls in the main body of the code
! x Do firstcalls in the kernal calculation
! x Put in environment variables: verbose and resolution - just verbose, resolution is hardwired
!
! This calculates a relativistic transfer function for an on-axis lamppost source
! including as many effects as possible:
! 1) The shifting of cut-off energy for different radii and for the observer
! 2) The adjustment of the continuum flux for g^2
! 3) The angular dependence of the reflection spectrum (viewing angle)
! 4) If possible, the dependence on incident angle (how well can this be mimicked by tweaking the ionization?)

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
      subroutine genreltrans(Cp,ear,ne,param,ifl,photar)
      use dyn_gr
      implicit none
      integer nro,nphi,ndelta,nex,i,nf,ifl,ne,ReIm,nfsave
      integer nrosave,nphisave,verbose,mubin,sdbin
      integer me,ge,xe,Cp
!      parameter (ndelta=1000,nex=2**12,mex=1,gex=1,xex=1)
      parameter (ndelta=1000,nex=2**12)
      double precision a,h,Gamma,inc,pi,rout,rmin,disco,muobs,rin
      double precision Mass,flo,fhi,dlogf,dgsofac,zcos,frobs,honr,rnmax,d
      double precision fhisave,flosave,rh,hsave,rinsave,frrel
      real afac,fc,param(19),ear(0:ne),gso,direct
      real Afe,Ecut_s,Ecut_obs,logxi,xillpar(7),E,dE,earx(0:nex),Emax,Emin,dloge
      real reline(nex),imline(nex),photarx(nex),reconv(nex),imconv(nex)
      real reconvmu(nex),imconvmu(nex),mue,sdmin,sdmax,gsd,ximin,ximax
      real phase,t0,t1,ReSx(nex),ImSx(nex),ReS(ne),ImS(ne),photar(ne)
      real paramsave(19),contx(nex),frac,phiA
      real ReGx(nex),ImGx(nex),sum
      complex,dimension(:,:,:,:),allocatable :: transe,transea
      logical firstcall,needtrans,needconv,needresp
      integer xbin,xbinhi,myenv,Cpsave
      real logxir

!variable for the grid reading
      integer :: lrec,irec,nphi_grid,nro_grid,spin_dim,mu_dim,check
      integer :: s_lo,s_hi,m_lo,m_hi,xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh
      double precision :: honr_grid,spin_lo,spin_hi,mu_lo,mu_hi,spin_start,spin_end,mu_start,mu_end,ave_weight2D
      real :: ximin_sl_ml,ximax_sl_ml,ximin_sl_mh,ximax_sl_mh,ximin_sh_ml,ximax_sh_ml&
           ,ximin_sh_mh,ximax_sh_mh
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

      real :: reline_reb(ne)
      
      data firstcall /.true./
      data nrosave,nphisave /0,0/
      data needresp/.true./
      data Cpsave/2/
      save firstcall,Emax,Emin,dloge,earx
      save paramsave,fhisave,flosave,nfsave,nrosave,nphisave
      save frobs,hsave,rinsave,sdmin,sdmax,frrel,Cpsave
      save transe,transea,transe_11,transe_12,transe_21,transe_22,transe_11a,transe_12a,transe_21a,transe_22a
      save reconv,imconv,reconvW1,imconvW1,reconvW1a,imconvW1a,check
      pi = acos(-1.d0)
      ifl = 1
      
! Settings
      nro   = 300    !resolution variables - these could be made parameters
      nphi  = 300    !  "
      ! if( nro .gt. nmax ) nro = nmax
      ! if( nphi .gt. nmax ) nphi = nmax
      dlogf = 0.0073  !This is a resolution parameter (base 10)

! Call environment variables
      verbose = myenv("REV_VERB",0)     !Set verbose level
      me      = myenv("MU_ZONES",1)     !Set number of mu_e zones used
      ge      = myenv("ECUT_ZONES",1)   !Set number of Ecut zones used
      xe      = myenv("ION_ZONES",100)    !Set number of ionisation zones used
      
! ! Make sure they haven't exceeded their maximum allowed values
      ! call sizecheck(me,mex)
      ! call sizecheck(ge,gex)
      ! call sizecheck(xe,xex)
      
! Initialise
      call initialiser(firstcall,Emin,Emax,nex,dloge,earx,rnmax,d,needtrans,check&
     ,nphi,nro,honr_grid,spin_start,spin_end,mu_start,mu_end,spin_dim,mu_dim)

!check if the grid values are the same one of the model
      
         if ( check .ne. 0 .and. honr_grid .ne. honr ) then
            write(*,*) 'grid has a different honr!'
            write(*,*) 'honr of the grid is ', honr_grid
            stop
         endif

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
      honr     = dble( param(11) )
      afac     = param(12)
      Mass     = dble( param(13) )
      flo      = dble( param(14) )
      fhi      = dble( param(15) )
      ReIm     = int( param(16) )
      phiA     = param(17)
      phiB     = param(18)
      g        = param(19)
      
      muobs = cos( inc * pi / 180.d0 )

      
!Work out how many frequencies to average over
      fc = 0.5 * ( flo + fhi )
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
      Ecut_s = (1.0+zcos) * Ecut_obs / gso
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
        if( abs( h - hsave ) .gt. 1e-7 ) needtrans = .true.
        if( abs( rin - rinsave ) .gt. 1e-7 ) needtrans = .true.
        if( nro .ne. nrosave ) needtrans = .true.
        if( nphi .ne. nphisave ) needtrans = .true.
        if( nf .ne. nfsave ) needtrans = .true.
        if( abs( fhi - fhisave ) .gt. 1e-7 ) needtrans = .true.
        if( abs( flo - flosave ) .gt. 1e-7 ) needtrans = .true.
      end if

      needtrans = .true.
      
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
            
!            write(*,*) 'GRID extraction'
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
            call strans(spin_lo,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ndelta,nex,dloge,&
                 nf,fhi,flo,me,ge,xe,logxi,sdmin_sl_ml,sdmax_sl_ml,ximin_sl_ml,ximax_sl_ml&
                 ,transe_11,transe_11a,frobs_sl_ml,frrel_sl_ml,xbinhi_sl_ml)

!Repeat this for all nthe combinations of spin and inclination (4 times)
            irec = 3*mu_dim*(s_lo-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_lo,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ndelta,nex,dloge,&
                 nf,fhi,flo,me,ge,xe,logxi,sdmin_sl_mh,sdmax_sl_mh,ximin_sl_mh,ximax_sl_mh&
                 ,transe_21,transe_21a,frobs_sl_mh,frrel_sl_mh,xbinhi_sl_mh)

            irec = 3*mu_dim*(s_hi-1) + 3*m_lo - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_lo,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ndelta,nex,dloge,&
                 nf,fhi,flo,me,ge,xe,logxi,sdmin_sh_ml,sdmax_sh_ml,ximin_sh_ml,ximax_sh_ml&
                 ,transe_12,transe_12a,frobs_sh_ml,frrel_sh_ml,xbinhi_sh_ml)
            
            irec = 3*mu_dim*(s_hi-1) + 3*m_hi - 1
            read(98,rec=irec) re1
            read(98,rec=irec+1) taudo1
            read(98,rec=irec+2) pem1
            call strans(spin_hi,h,mu_hi,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ndelta,nex,dloge,&
                 nf,fhi,flo,me,ge,xe,logxi,sdmin_sh_mh,sdmax_sh_mh,ximin_sh_mh,ximax_sh_mh&
                 ,transe_22,transe_22a,frobs_sh_mh,frrel_sh_mh,xbinhi_sh_mh)
   
            ! write(*,*)spin_lo,h,mu_hi,gamma,rin,rout,honr,zcos,nro,nphi,ndelta,nex,dloge &
            !      ,nf,fhi,flo,nmu,refl_sh_ml

!Weighted average of the reflection fraction 
            frobs = ave_weight2D(a,spin_lo,spin_hi,muobs,mu_lo,mu_hi,frobs_sl_ml,frobs_sl_mh,frobs_sh_ml,frobs_sh_mh)
!Determin the maximum value of the inonisation binning 
            xbinhi = max(xbinhi_sl_ml,xbinhi_sl_mh,xbinhi_sh_ml,xbinhi_sh_mh)

         else 
!            write(*,*) 'transfer calculation'
            if(.not. allocated(transe) ) allocate(transe(nex,me,ge,xe))
            if(.not. allocated(transea) ) allocate(transea(nex,me,ge,xe))

!Calculate the Kernel for the given parameters
            status_re_tau = .true.
            call strans(a,h,muobs,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ndelta,nex,dloge,&
                 nf,fhi,flo,me,ge,xe,logxi,sdmin,sdmax,ximin,ximax,transe,transea,frobs,frrel,xbinhi)
         endif         
      end if
      
      if( verbose .gt. 0 ) write(*,*)"Observer's reflection fraction=",afac*frobs
!      if( verbose .gt. 0 ) write(*,*)"Relxill reflection fraction=",frrel
      
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
        if( verbose .gt. 0 ) call sourcelum(nex,earx,contx,real(mass))
!Now reflection
        xillpar(7) = -1.0       !reflection fraction of 1        
        reconv = 0.0
        imconv = 0.0
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
              xillpar(6) = acos( mue ) * 180.0 / pi
              if( me .eq. 1 ) xillpar(6) = real( inc )

!Call xillver
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
                 ! write(*,*) 'interpolation for the GRID'
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
     
      !Re-bin onto input grid
      call rebinE(earx,reline,nex,ear,reline_reb,ne)
      
! Calculate phiA from instrument response - if this option is set to on      
      call phaseA(nex,earx,contx,reconv,imconv,gso,zcos,Gamma,afac,phiA)
              
      !Add on continuum (and include boosting fudge factor)
!       do i = 1,nex
      !         E = 0.5 * ( earx(i) + earx(i-1) )
!         dE = earx(i) - earx(i-1)
!         direct  = contx(i) / dE * ( gso / (1.0+zcos) )**(2+Gamma)
!         ReSx(i) = direct + afac * reconv(i) / dE
!         ImSx(i) = afac * imconv(i) / dE
!         ReGx(i) = cos(phiA) * ReSx(i) - sin(phiA) * ImSx(i)
!         ImGx(i) = sin(phiA) * ReSx(i) + cos(phiA) * ImSx(i)
! !        write(300,*)E,dE,E**2*ReSx(i),E**2*direct,E**2*afac*reconv(i)/dE
!       end do

!Add on continuum (and include boosting fudge factor) + consider the non-linear effects
      do i = 1,nex
        E = 0.5 * ( earx(i) + earx(i-1) )
        dE = earx(i) - earx(i-1)
        direct  = contx(i) / dE * ( gso / (1.0+zcos) )**(2+Gamma)

        ReGx(i) = cos(phiA) * (direct + afac*reconv(i)/dE) - sin(phiA)*afac*imconv(i)/dE + &
        g * (cos(phiB) * ( log(E*(1+zcos)/gso)*direct - afac*reconvW1a(i)/dE - afac*reconvW1(i)/dE)+&
        sin(phiB)* (imconvW1a(i)/dE + imconvW1(i)/dE) )
        ImGx(i) = sin(phiA) * (direct + afac*reconv(i)/dE) + cos(phiA)*afac*imconv(i)/dE + &
        g * (sin(phiB) * ( log(E*(1+zcos)/gso)*direct - afac*reconvW1a(i)/dE - afac*reconvW1(i)/dE)-&
        cos(phiB) * (afac*imconvW1a(i)/dE + afac*imconvW1(i)/dE) )
        
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
        write(*,*)"Warning, only fit to data with ReIm=1 & 2"
      else if( ReIm .eq. 4 )then   !Phase lag (cycles)
        do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          phase = atan2( ImS(i) , ReS(i) )
          photar(i) = phase / (2.0*pi) * dE
        end do
        write(*,*)"Warning, only fit to data with ReIm=1 & 2"
      else                         !Time lag (seconds)
         do i = 1,ne
          E = 0.5 * ( ear(i) + ear(i-1) )
          dE = ear(i) - ear(i-1)
          phase = atan2( ImS(i) , ReS(i) )
          photar(i) = phase / (2.0*pi*fc) * dE
       end do
        write(*,*)"Warning, only fit to data with ReIm=1 & 2"
      end if

      nrosave   = nro
      nphisave  = nphi
      fhisave   = fhi
      flosave   = flo
      nfsave    = nf
      rinsave   = rin
      hsave     = h
      paramsave = param
      Cpsave    = Cp

      end subroutine genreltrans
!-----------------------------------------------------------------------




      


! !-----------------------------------------------------------------------
!       function na(r,rin)
! ! Number density for zone a of a Shakura & Sunyaev disk (eqn 2.11 in SS73)
!       implicit none
!       double precision na,r,rin,x,f
!       x = r / rin
!       f = 1.0 - x**(-0.5)
!       na = r**1.5 / f**2
!       return
!       end function na
! !-----------------------------------------------------------------------

