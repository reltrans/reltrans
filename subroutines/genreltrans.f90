!-----------------------------------------------------------------------
subroutine genreltrans(Cp, dset, nlp, ear, ne, param, ifl, photar)
! All reltrans flavours are calculated in this subroutine.
! Cp and dset are the settings:
! |Cp|=1 means use cut-off power-law, |Cp|=2 means use nthcomp
! Cp>1 means there is a density parameter, Cp<1 means density is hardwired  
! dset=0 means ionisation is a parameter, dset=1 means ionization is calculated
! from distance. What to do about ION_ZONES=1 in the distance model?

! The parameter array has 27 parameters. No one model actually has 27
! parameters. In each model, some of these parameters are hardwired, but
! the parameters must be sorted into the param(1:27) array for this subroutine.
  
!    Arg:
! 
!  Internal variables:
!         constants:
!         pi: greek pi
!         rnmax: maximum radius to consider GR effects
!         nphi, rno: resolution variables, number of pixels on the observer's camera(b and phib)
!         Emax, Emin: minimum and maximum range of the internal energy grid which is different than the xspec one
!         dlogf: resolution parameter of the frequency grid
!         dyn:   limit to check the saved values
!         ionvar: sets the ionisation variation (1 = w/ ion var; 0 = w/o ion var)
  
    use dyn_gr
    use conv_mod
    implicit none
    !Constants
    integer         , parameter :: nphi = 200, nro = 200!, ionvar! = 1 
    real            , parameter :: Emin = 1e-2, Emax = 3e3, dyn = 1e-7
    double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0, &
                                   dlogf = 0.09 !This is a resolution parameter (base 10)       
    !Args:
    integer, intent(inout) :: ifl
    integer, intent(in)    :: Cp, dset, ne, nlp
    real   , intent(inout) :: param(32)
    real   , intent(out)   :: photar(ne)  
    !Variables of the subroutine
    !initializer
    integer          :: verbose, me, xe, m, ionvar, refvar
    logical          :: firstcall, needtrans, needconv, test
    double precision :: d
    !Parameters of the model:
    double precision :: h(nlp), a, inc, rin, rout, zcos, Gamma, honr, muobs 
    real             :: logxi, Afe, lognep, Ecut_obs, Ecut_s, Dkpc, Anorm, beta_p
    real             :: Nh, boost, Mass, floHz, fhiHz, DelA, DelAB(nlp), g(nlp)
    integer          :: ReIm, resp_matr
    double precision :: qboost,b1,b2, eta, eta_0
    !internal frequency grid
    integer          :: nf 
    real             :: f, fac
    double precision :: fc, flo, fhi
    ! internal energy grid (nex) and output/xspec (ne) energy grid
    real             :: E, dE, dloge
    real             :: earx(0:nex)   
    real             :: ear(0:ne)
    ! internal frequency grid, for when we do lag/frequency spectra
    integer           :: fbinx 
    real, allocatable :: fix(:)
    !relativistic parameters and limit on rin and h
    double precision :: rmin, rh 
    double precision :: gso(nlp),tauso(nlp),cosdelta_obs(nlp),height(nlp),contx_int(nlp)
    !lens needs to be allocatable to save it. 
    double precision, allocatable :: lens(:),frobs(:),frrel(:)
    !TRANSFER FUNCTIONS and Cross spectrum dynamic allocation + variables
   ! complex, dimension(:,:,:,:,:), allocatable :: transe, transea
    complex, dimension(:,:,:,:,:), allocatable :: ker_W0,ker_W1,ker_W2,ker_W3
    real   , dimension(:,:,:)    , allocatable :: ReW0,ImW0,ReW1,ImW1
    real   , dimension(:,:,:)    , allocatable :: ReW2,ImW2,ReW3,ImW3
    real   , dimension(:,:)      , allocatable :: ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG                                                
    !double precision :: frobs(nlp), frrel(nlp)  !reflection fraction variables (verbose)
    !Radial and angle profile 
    integer                       :: mubin, rbin, ibin
    double precision, allocatable :: logxir(:),gsdr(:), logner(:)
    real    :: contx(nex,nlp)
    real    :: mue, logxi0, reline_w0(nlp,nex), imline_w0(nlp,nex), photarx(nex), photerx(nex)
    real    :: absorbx(nex), ImGbar(nex), ReGbar(nex)
    real    :: ReGx(nex),ImGx(nex),ReS(ne),ImS(ne)
    !variable for non linear effects
    integer ::  DC, ionvariation
    real    :: photarx_1(nex), photarx_2(nex), photarx_delta(nex), photarx_dlogxi(nex)
    real    :: reline_w1(nlp,nex),imline_w1(nlp,nex),reline_w2(nlp,nex),imline_w2(nlp,nex)
    real    :: reline_w3(nlp,nex),imline_w3(nlp,nex)
    real    :: dlogxi1, dlogxi2, Gamma1, Gamma2, DeltaGamma  
    !SAVE 
    integer          :: nfsave, Cpsave
    real             :: paramsave(32)
    double precision :: fhisave, flosave
    !Functions
    integer          :: i, j
    double precision :: disco, dgsofac
    ! New  
    double precision :: fcons,get_fcons,contx_temp!,ell13pt6,lacc,get_lacc,
    real             :: Gamma0,logne,Ecut0,thetae,logxiin
    integer          :: Cp_cont
    real time_start,time_end        !runtime stuff
 
    data firstcall /.true./
    data Cpsave/2/
    data nfsave /-1/  
    !Save the first call variables
    save firstcall, dloge, earx, me, xe, d, verbose
    save paramsave, fhisave, flosave, nfsave, refvar, ionvar
    save frobs, frrel, Cpsave, needtrans, lens
    save ker_W0, ker_W1, ker_W2, ker_W3, logxir, gsdr, logner
    save ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3
    save ReSraw,ImSraw,ReSrawa,ImSrawa,ReGrawa,ImGrawa,ReG,ImG

    ifl = 1
    ! Initialise some parameters 
    call initialiser(firstcall,Emin,Emax,dloge,earx,rnmax,d,needtrans,me,xe,refvar,ionvar,verbose, test)

    if (test) then 
       call FNINIT
    endif
    
    !Allocate dynamically the array to calculate the trasfer function 
    if (.not. allocated(re1)) allocate(re1(nphi,nro))
    if (.not. allocated(taudo1)) allocate(taudo1(nphi,nro))
    if (.not. allocated(pem1)) allocate(pem1(nphi,nro))
    
    !Note: the two different calls are because for the double lP we set the temperature from the coronal frame(s), but for the single
    !LP we use the temperature in the observer frame
    if (nlp .eq. 1) then
        call set_param(dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut_obs,&
                       eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                       g,Anorm,resp_matr,refvar,verbose)        
    else 
        call set_param(dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut_s,&
                       eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                       g,Anorm,resp_matr,refvar,verbose) 
    end if 
    
    muobs = cos( inc * pi / 180.d0 )

    !this needs to go in a subroutine - model_mode or something
    !rework this logic so that low frequencies always result in time independent spectrum, not just for reim<7
    if( ReIm .eq. 7 .and. fhiHz .gt. tiny(fhiHz) .and. floHz .gt. tiny(floHz)) then
        !set up frequency grid in Hz if using lag frequency mode, depending on whether we're looking at AGN or XRBs
        if( Mass .gt. 1000 ) then
            floHz = 1.e-5
            fhiHz = 5.e-2 
        else
            floHz = 0.07
            fhiHz = 700.
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
        !Note that the frequency grid is using a higher resolution since it's what we care about in this mode 
        !TBD: find a way to reduce energy resolution since we don't need it, otherwise the runtime is insanely slow  
        nf = ceiling( log10(fhiHz/floHz) / 0.01 )
        allocate(fix(0:nf))
        do fbinx = 0, nf 
            fix(fbinx) = floHz * (fhiHz / floHz)**( real(fbinx) / real(nf) )
        end do
    else 
        !if doing lag-energy spectra, just work out how many frequencies to average over 
        fc = 0.5d0 * ( floHz + fhiHz )
        nf = ceiling( log10(fhiHz/floHz) / dlogf )
        if( fhiHz .lt. tiny(fhiHz) .or. floHz .lt. tiny(floHz) )then
            fhiHz = 0.d0
            floHz = 0.d0
            nf    = 1
        end if
        !Convert frequency bounds from Hz to c/Rg (now being more accurate with constants)
        fhi   = dble(fhiHz) * 4.92695275718945d-06 * Mass
        flo   = dble(floHz) * 4.92695275718945d-06 * Mass
    end if   

    !Decide if this is the DC component/time averaged spectrum or not
    if( flo .lt. tiny(flo) .or. fhi .lt. tiny(fhi) )then
        DC     = 1
        g      = 0.0
        DelAB  = 0.0
        DelA   = 0.0
        ReIm   = 1
        eta    = eta_0
        beta_p = 1. !this is an ugly hack for the double LP model, to calculate the time-averaged spectrum
    else
        DC     = 0
        boost  = abs(boost)
    end if
    !this could go into a subroutine -- just put it in set_params?
    !Set minimum r (ISCO) and convert rin and h to rg
    if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
    rmin   = disco( a )
    if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
    rh     = 1.d0+sqrt(1.d0-a**2)
    if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
    if( rin .lt. rmin )then
        write(*,*)"Warning! rin<ISCO! Set to ISCO"
        rin = rmin
    end if
    do m=1,nlp 
        if( h(m) .lt. 0.d0 ) h(m) = abs(h(m)) * rh
        if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h(m)
        if( h(m) .lt. 1.5d0*rh )then
            write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
            h(m) = 1.5d0 * rh
        end if 
    end do

    !Determine if I need to calculate the kernel 
    call need_check(Cp,Cpsave,param,paramsave,fhi,flo,fhisave,flosave,nf,nfsave,needtrans,needconv)

    ! Allocate arrays that depend on frequency
    if( nf .ne. nfsave )then
        if( allocated(ker_W0 ) ) deallocate(ker_W0 )
        if( allocated(ker_W1) ) deallocate(ker_W1 )
        if( allocated(ker_W2 ) ) deallocate(ker_W2 )
        if( allocated(ker_W3 ) ) deallocate(ker_W3 )
        allocate( ker_W0(nlp,nex,nf,me,xe) )
        allocate( ker_W1(nlp,nex,nf,me,xe) )
        allocate( ker_W2(nlp,nex,nf,me,xe) )
        allocate( ker_W3(nlp,nex,nf,me,xe) )
        if( allocated(ReW0) ) deallocate(ReW0)
        if( allocated(ImW0) ) deallocate(ImW0)
        if( allocated(ReW1) ) deallocate(ReW1)
        if( allocated(ImW1) ) deallocate(ImW1)
        if( allocated(ReW2) ) deallocate(ReW2)
        if( allocated(ImW2) ) deallocate(ImW2)
        if( allocated(ReW3) ) deallocate(ReW3)
        if( allocated(ImW3) ) deallocate(ImW3)
        allocate( ReW0(nlp,nex,nf) )
        allocate( ImW0(nlp,nex,nf) )
        allocate( ReW1(nlp,nex,nf) )
        allocate( ImW1(nlp,nex,nf) )
        allocate( ReW2(nlp,nex,nf) )
        allocate( ImW2(nlp,nex,nf) )
        allocate( ReW3(nlp,nex,nf) )
        allocate( ImW3(nlp,nex,nf) )
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
  
    !allocate lensing/reflection fraction arrays if necessary
    if( needtrans ) then
        if( allocated(lens) ) deallocate( lens )
        allocate (lens(nlp))
        if( allocated(frobs) ) deallocate( frobs )
        allocate (frobs(nlp))
        if( allocated(frrel) ) deallocate( frrel )
        allocate (frrel(nlp))
    end if

    !set up the continuum spectrum plus relative quantities (cutoff energies, lensing/gfactors, luminosity, etc)
    call init_cont(nlp,a,h,zcos,Ecut_s,Ecut_obs,logxi, lognep, gso,&
                   muobs,lens,tauso,cosdelta_obs,Cp_cont,Cp,fcons,Gamma,&
                   Dkpc,Mass,earx,Emin,Emax,contx,dlogE,verbose,dset,Anorm,contx_int,eta)    

    if (verbose .gt. 2) call CPU_TIME (time_start)
    if( needtrans )then
        !Allocate arrays for kernels   
        if (allocated (logxir)) deallocate (logxir)
        allocate (logxir(xe))
        if (allocated (gsdr)) deallocate (gsdr)
        allocate (gsdr(xe))
        if (allocated (logner)) deallocate (logner)
        allocate (logner(xe))
        !Calculate the Kernel for the given parameters
        status_re_tau = .true.
        call rtrans(verbose,dset,nlp,a,h,gso,muobs,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,eta_0,&
                    fcons,contx_int,tauso,lens,cosdelta_obs,nro,nphi,nex,dloge,nf,fhi,flo,me,xe,logxi,lognep,&
                    ker_W0,ker_W1,ker_W2,ker_W3,frobs,frrel,logxir,gsdr,logner)         
    end if
    if( verbose .gt. 2 ) then
        call CPU_TIME (time_end)
        print *, 'Transfer function runtime: ', time_end - time_start, ' seconds'
    end if
  
    !do this for each lamp post, then find some sort of weird average?
    if( verbose .gt. 0) write(*,*)"Observer's reflection fraction for each source:",boost*frobs
    if( verbose .gt. 0) write(*,*)"Relxill reflection fraction for each source:",frrel    
    
    if( verbose .gt. 2) call CPU_TIME (time_start)  
    if( needconv )then
        !needtrans = .false.
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
        !Get logxi values corresponding to Gamma1 and Gamma2
        call xilimits(nex,earx,nlp,contx,DeltaGamma,real(gso),real(lens),real(zcos),dlogxi1,dlogxi2)
        !Set the ion-variation to 1, there is an if inside the radial loop to check if either the ionvar is 0 or the logxi is 0 to
        !set ionvariation to 0  it is important that ionvariation is different than ionvar because ionvar  is used also later in
        !the rawS subroutine to calculate the cross-spectrum
        ionvariation = 1
        !Loop over radius, emission angle and frequency
        do rbin = 1, xe  !Loop over radial zones
            !Set parameters with radial dependence
            Gamma0 = real(Gamma)
            logne  = logner(rbin)
            Ecut0  = real( gsdr(rbin) ) * Ecut_s
            logxi0 = real( logxir(rbin) )    
            if( xe .eq. 1 )then
                Ecut0  = Ecut_s
                logne  = lognep
                logxi0 = logxi
            end if
            !Avoid negative values of the ionisation parameter 
            if (logxi0 .eq. 0.0 .or. ionvar .eq. 0) then
                ionvariation = 0.0
            end if
            do mubin = 1, me      !loop over emission angle zones
                !Calculate input emission angle
                mue    = ( real(mubin) - 0.5 ) / real(me)
                thetae = acos( mue ) * 180.0 / real(pi)
                if( me .eq. 1 ) thetae = real(inc)
                !Call restframe reflection model
                call rest_frame(earx,nex,Gamma0,Afe,logne,Ecut0,logxi0,thetae,Cp,photarx)
                !NON LINEAR EFFECTS
                if (DC .eq. 0) then 
                    !Gamma variations
                    logxiin = logxi0 + ionvariation * dlogxi1
                    call rest_frame(earx,nex,Gamma1,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_1)
                    logxiin = logxi0 + ionvariation * dlogxi2
                    call rest_frame(earx,nex,Gamma2,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_2)
                    photarx_delta = (photarx_2 - photarx_1)/(Gamma2-Gamma1)
                    !xi variations
                    logxiin = logxi0 + ionvariation * dlogxi1
                    call rest_frame(earx,nex,Gamma0,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_1)
                    logxiin = logxi0 + ionvariation * dlogxi2
                    call rest_frame(earx,nex,Gamma0,Afe,logne,Ecut0,logxiin,thetae,Cp,photarx_2)
                    photarx_dlogxi = 0.434294481 * (photarx_2 - photarx_1) / (dlogxi2-dlogxi1) !pre-factor is 1/ln10           
                end if
                !Loop through frequencies and lamp posts
                do j = 1,nf
                    do i = 1,nex
                        do m=1,nlp
                            reline_w0(m,i) = real( ker_W0(m,i,j,mubin,rbin) )
                            imline_w0(m,i) = aimag( ker_W0(m,i,j,mubin,rbin) )
                            reline_w1(m,i) = real( ker_W1(m,i,j,mubin,rbin) )
                            imline_w1(m,i) = aimag( ker_W1(m,i,j,mubin,rbin) )
                            reline_w2(m,i) = real( ker_W2(m,i,j,mubin,rbin) )
                            imline_w2(m,i) = aimag( ker_W2(m,i,j,mubin,rbin) )
                            reline_w3(m,i) = real( ker_W3(m,i,j,mubin,rbin) )
                            imline_w3(m,i) = aimag( ker_W3(m,i,j,mubin,rbin) )
                        end do  
                    end do
                    !always: convolution for reverberation/DC spectrum
                    !TBD: add flag here to do this convolution if no reflection time, or different convolution with complex
                    !xillver if tref > 0 or something.

                    if (test) then 
                       call conv_one_FFT(dyn,photarx,reline_w0,imline_w0,ReW0(:,:,j),ImW0(:,:,j),DC,nlp)
                       if(DC .eq. 0 .and. refvar .eq. 1) then
                          call conv_one_FFT(dyn,photarx,reline_w1,imline_w1,ReW1(:,:,j),ImW1(:,:,j),DC,nlp)
                          call conv_one_FFT(dyn,photarx_delta,reline_w2,imline_w2,ReW2(:,:,j),ImW2(:,:,j),DC,nlp)
                       end if
                       if(DC .eq. 0 .and. ionvar .eq. 1) then
                          call conv_one_FFT(dyn,photarx_dlogxi,reline_w3,imline_w3,ReW3(:,:,j),ImW3(:,:,j),DC,nlp)
                       end if
                    else
                       call conv_one_FFTw(dyn,photarx,reline_w0,imline_w0,ReW0(:,:,j),ImW0(:,:,j),DC,nlp)
                       if(DC .eq. 0 .and. refvar .eq. 1) then
                          call conv_one_FFTw(dyn,photarx,reline_w1,imline_w1,ReW1(:,:,j),ImW1(:,:,j),DC,nlp)
                          call conv_one_FFTw(dyn,photarx_delta,reline_w2,imline_w2,ReW2(:,:,j),ImW2(:,:,j),DC,nlp)
                       end if
                       if(DC .eq. 0 .and. ionvar .eq. 1) then
                          call conv_one_FFTw(dyn,photarx_dlogxi,reline_w3,imline_w3,ReW3(:,:,j),ImW3(:,:,j),DC,nlp)
                       end if
                    endif
                    !old call: always convolve every single transfer function in one go
                    !call conv_all_FFTw(dyn,photarx,photarx_delta,photarx_dlogxi,reline_w0,imline_w0,reline_w1,imline_w1,&
                    !     reline_w2,imline_w2,reline_w3,imline_w3,ReW0(:,:,j),ImW0(:,:,j),ReW1(:,:,j),ImW1(:,:,j),&
                    !     ReW2(:,:,j),ImW2(:,:,j),ReW3(:,:,j),ImW3(:,:,j),DC,nlp)                        
                end do
            end do
        end do
    end if
    if( verbose .gt. 2 ) then
        call CPU_TIME (time_end)
        print *, 'Convolutions runtime: ', time_end - time_start, ' seconds' 
    endif

    ! do i = 1, nex
    !    write(32,*) (earx(i) + earx(i-1))*0.5, &
    !          contx(i,1)/(earx(i) - earx(i-1)) *  ((earx(i) + earx(i-1))*0.5)**2 
    ! enddo
       
    ! Calculate absorption 
    call tbabs(earx,nex,nh,Ifl,absorbx,photerx)

    !TBD coherence check - if zero coherence between lamp posts, call a different subroutine 
    if( ReIm .eq. 7 ) then
        !tbd - implement zero cohernece in lag_freq
        if(nlp .gt. 1 .and. beta_p .eq. 0. ) then
            call lag_freq_nocoh(nex,earx,nf,fix,real(flo),real(fhi),Emin,Emax,nlp,contx,absorbx,real(tauso),real(gso),&
                                ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),boost,&
                                g,DelAB,ionvar,ReGbar,ImGbar)
        else
            call lag_freq(nex,earx,nf,fix,real(flo),real(fhi),Emin,Emax,nlp,contx,absorbx,real(tauso),real(gso),&
                          ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),beta_p,&
                          boost,g,DelAB,ionvar,ReGbar,ImGbar)        
        end if
    else if (nlp .gt. 1 .and. beta_p .eq. 0.) then
        call rawG(nex,earx,nf,real(flo),real(fhi),nlp,contx,absorbx,real(tauso),real(gso),ReW0,ImW0,&
                  ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),boost,ReIm,g,DelAB,&
                  ionvar,DC,resp_matr,ReGrawa,ImGrawa)                
    else
        !Calculate raw FT of the full spectrum without absorption
        call rawS(nex,earx,nf,real(flo),real(fhi),nlp,contx,real(tauso),real(gso),ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                  real(h),real(zcos),real(Gamma),real(eta),beta_p,boost,g,DelAB,ionvar,DC,ReSraw,ImSraw)
        !Include absorption in the model
        do j = 1, nf
            do i = 1, nex
                ReSrawa(i,j) = ReSraw(i,j) * absorbx(i)
                ImSrawa(i,j) = ImSraw(i,j) * absorbx(i)
            end do
        end do        
    end if

    if( DC .eq. 1 )then
        !Norm is applied internally for DC/time averaged spectrum component of dset=1
        !No need for the immaginary part in DC
        do i = 1, nex
            ReGbar(i) = (Anorm/real(1.+eta)) * ReSrawa(i,1)
        end do
    else if (ReIm .eq. 7) then     
        !if calculating the lag-frequency spectrum, just rebin the arrays 
        call rebinE(fix, ReGbar, nf, ear, ReS, ne)
        call rebinE(fix, ImGbar, nf, ear, ImS, ne)     
    else 
        !In this case, calculate the lag-energy spectrum
        !Calculate raw cross-spectrum from Sraw(E,\nu) and the reference band parameters
        !note: this must be done by rawG for two incoherent lamp posts, hence the skip below
        if(nlp .eq. 1 .or. beta_p .ne. 0.) then
            if (ReIm .gt. 0.0) then
                call propercross(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa, resp_matr)
            else
                call propercross_NOmatrix(nex, nf, earx, ReSrawa, ImSrawa, ReGrawa, ImGrawa)
            endif
        end if
        !Apply phase correction parameter to the cross-spectral model (for bad calibration)
        !this is where coherence = 0 or = 1 cases merge back 
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
        !This means that norm for the AC components in the dset=1 model is power in squared fractional rms format
        !note: the factor eta is to have the same normalization as the single LP model, it's 100% arbitrary
        ReGbar = ReGbar * fac * (Anorm/real(1.+eta))**2  
        ImGbar = ImGbar * fac * (Anorm/real(1.+eta))**2  
    end if
    
    !Write output depending on ReIm parameter
    if( ReIm .eq. 7 ) then
        do i=1,ne 
            dE = ear(i) - ear(i-1) 
            photar(i) = atan2(ImS(i),ReS(i))/(pi*(ear(i) + ear(i-1)))*dE
        end do
    else if( abs(ReIm) .le. 4 )then
        call crebin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS) !S is in photar form
        ! do i = 1, ne
        !    write(98,*) (ear(i) + ear(i-1))*0.5, ReS(i)/(ear(i) - ear(i-1)) *  ((ear(i) + ear(i-1))*0.5)**2
        ! enddo
        
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
        call cfoldandbin(nex,earx,ReGbar,ImGbar,ne,ear,ReS,ImS,resp_matr) !S is count rate
        if( abs(ReIm) .eq. 5 )then        !Modulus
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                photar(i) = sqrt( ReS(i)**2 + ImS(i)**2 ) 
            end do
        else if( abs(ReIm) .eq. 6 )then   !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                photar(i) = atan2( ImS(i) , ReS(i) ) / ( 2.0*pi*fc ) * dE
            end do
        end if
    end if

    if (verbose .gt. 1 .and. abs(ReIm) .gt. 0 .and. ReIm .lt. 7) then
        if( DC .eq. 0 .and. beta_p .eq. 0) then
           call write_components(ne,ear,nex,earx,nf,real(flo),real(fhi),nlp,contx,absorbx,real(tauso),real(gso),&
                                  ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,real(h),real(zcos),real(Gamma),real(eta),&
                                  beta_p,boost,floHz,fhiHz,ReIm,DelA,DelAB,g,ionvar,resp_matr)
        !catch case here for coherence = 0 or 1
        end if                
        !this writes the full model as returned to Xspec 
        !note that xspec gets output in e.g. lags*dE, and we want just the lags, so a factor dE needs to be included
        !add writing of components for lag frequency spectrum
        open (unit = 14, file = 'Output/Total.dat', status='replace', action = 'write')     
        do i = 1,ne 
            dE = ear(i) - ear(i-1)
            write (14,*) (ear(i)+ear(i-1))/2., photar(i)/dE      
        end do 
        close(14)  
        !print continuum for both single and multiple LPs REDO THIS 
        open (unit = 24, file = 'Output/Continuum_spec.dat', status='replace', action = 'write')
        do i=1,nex
            dE = earx(i) - earx(i-1)
            if( nlp .eq. 1 ) then
                contx_temp = contx(i,1)/dE
            else
                contx_temp = 0.
                do m=1,nlp 
                    contx_temp = contx_temp + contx(i,m)
                end do
                contx_temp =  contx_temp/((1.+eta)*dE)      
            end if
            write (24,*) (earx(i)+earx(i-1))/2., contx_temp
        end do
        close(24)
    else if (ReIm .eq. 7) then
       open (unit = 14, file = 'Output/Total.dat', status='replace', action = 'write')     
        do i = 1,ne 
            dE = ear(i) - ear(i-1)
            write (14,*)  (ear(i)+ear(i-1))/2., photar(i)/dE        
        end do 
        close(14)
    endif 

    fhisave   = fhi
    flosave   = flo
    nfsave    = nf
    paramsave = param
    Cpsave    = Cp
  
end subroutine genreltrans
!-----------------------------------------------------------------------
