!-----------------------------------------------------------------------
subroutine rtrans(verbose,dset,nlp,spin,h,mu0,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,eta_0,&
                  fcons,nro,nphi,ne,dloge,nf,fhi,flo,me,xe,ker_W0,ker_W1,ker_W2,ker_W3,frobs,frrel)
    ! Code to calculate the transfer function for an accretion disk.
    ! This code first does full GR ray tracing for a camera with impact parameters < bmax
    ! It then also does straight line ray tracing for impact parameters >bmax
    ! It adds both up to produce a transfer function for a disk extending from rin to rout
    ! INPUT
    ! verbose               Decides whether to print radial scalings to file or not
    ! dset                  dset=1 means calculate ionization from distance, dset=0 means ignore distance
    ! nlp                   number of lamp post height considered
    ! spin,h,mu0,Gamma      Physical parameters (spin, source height(S), cos(inclination), photon index)
    ! nlp                   Number of lampposts considered (for now either 1 or 2)
    ! rin,rout,honr         Physical parameters (disk inner radius, outer radius & scaleheight)
    ! d,rnmax               Physical parameters (distance of the source, max radius for which GR ray tracing is used)
    ! zcos                  Cosmological redshift
    ! b1                    Linear coefficient of angular emissivity function
    ! b2                    Quadratic coefficient of angular emissivity function
    ! qboost                Asymmetry parameter of angular emissivity function
    ! fcons                 Used to calculate ionization from distance
    ! contx_int             Integral of the continuum flux over energy; needed to calculate the radial ionisation profile
    !                       in the double lamppost case
    ! nro,nphi              Number of pixels on the observer's camera (b and phib)
    ! ne, dloge             Number of energy bins and maximum energy (compatible with FFT convolution)
    ! nf,fhi,flo            nf = Number of logarithmic frequency bins used, range= flo to fhi
    ! me                    Number of mue bins
    ! xe                    Number of logr bins: bins 1:xe-1 are logarithmically spaced, bin xe is everything else
    ! OUTPUT
    ! ker_W0(nlp,ne,nf,me,xe)  Transfer function W0 - linear transfer function
    ! ker_W1(nlp,ne,nf,me,xe)  Transfer function W1 - one aspect of photon index variations
    ! ker_W2(nlp,ne,nf,me,xe)  Transfer function W2 - other aspect of photon index variations
    ! ker_W3(nlp,ne,nf,me,xe)  Transfer function W3 - ionization variations
    ! frobs                 Observer's reflection fraction
    ! frrel                 Reflection fraction defined by relxilllp
    use dyn_gr
    use blcoordinate
    use radial_grids
    use gr_continuum
    implicit none
    integer nro,nphi,ne,nf,me,xe,dset,nlp
    double precision spin,h(nlp),mu0,Gamma,rin,rout,zcos,fhi,flo,honr
    double precision b1,b2,qboost
    double precision fcons,cosdout(nlp)
    real dloge, lognep
    complex cexp!,transe(0:nlp,ne,nf,me,xe),transe_a(nlp,ne,nf,me,xe)

    integer i,j,odisc,n,gbin,rbin,mubin,l,m,k,nl
    double precision domega(nro),d,taudo,g,dlgfacthick,dFe(nlp),newtex,contx_int(nlp)
    !double precision rlp_column(ndelta),dcosdr_column(ndelta),tlp_column(ndelta),cosd_column(ndelta)
    double precision alpha,beta,cos0,sin0,phie,re,gsd(nlp)
    double precision tau(nlp),tausd,emissivity(nlp),cosfac,dglpfacthick,dareafac
    integer kk,fbin,get_index
    double precision rmin,disco,rfunc,scal,velocity(3),mudisk,sysfref
    double precision rnmax,rnmin,rn(nro),phin,mueff,dlogr,interper
    double precision fi(nf),dgsofac,sindisk,mue,demang,frobs(nlp),cosdin,frrel(nlp)
    double precision pnormer,mus,ptf,pfunc_raw,ang_fac
    integer nron,nphin,nrosav,nphisav,verbose
    double precision spinsav,musav,routsav,mudsav,rnn(nro),domegan(nro)
    integer get_env_int
    double precision lximax
    double precision eta_0
    logical dotrace

    !new stuff - move back above once it's implemented properly    
    complex ker_W0(nlp,ne,nf,me,xe),ker_W1(nlp,ne,nf,me,xe),ker_W2(nlp,ne,nf,me,xe),ker_W3(nlp,ne,nf,me,xe)
    real emisfac,thetafac(nlp),kfac,normfac
    
    !arrays to save the transfer function
    integer, parameter :: nt = 2**9
    integer            :: tbin
    double precision   :: tmin, tmax, sumresp, tar(0:nt), dlogt, dg, E
    double precision, allocatable :: resp(:,:)

    data nrosav,nphisav,spinsav,musav /0,0,2.d0,2.d0/
    save nrosav,nphisav,spinsav,musav,routsav,mudsav
       
    ! Settings/initialization
    nron     = 100
    nphin    = 100
    rmin     = disco( spin )
    scal     = 1.d0
    velocity = 0.d0
    mudisk   = honr / sqrt( honr**2 + 1.d0  )
    sindisk  = sqrt( 1.d0 - mudisk**2 )
    dfer_arr = 0.
    ker_W0 = 0.
    ker_W1 = 0.
    ker_W2 = 0.
    ker_W3 = 0.
    
    !set up saving the impulse response function if user desieres
!note: the ideal parameters to plot the transfer function are nro~=7000,nphi~=7000,nt~=2e9,nex~=2e10
    if (verbose .gt. 1) then    
        ! Create time grid in units of Rg
        tmin = 1.0
        tmax = 2.0e3
        dlogt = log10( tmax/tmin ) / float(nt)
        do i = 0,nt
            tar(i) = tmin * 10.0**( i * dlogt )
        end do
        ! Create energy grid optimised for plotting the transfer function (linear)
        dg = 2.0 / float(ne)
        !allocate and initialize impulse response function    
        if (.not. allocated(resp)) allocate(resp(ne, nt))
        resp = 0.0
        !add files to be printed here
        !open (unit = 104, file = 'Output/Impulse_2dImpulse.dat', status='replace', action = 'write')
        open (unit = 103, file = 'Output/Impulse_1dImpulseVsTime.dat', status='replace', action = 'write')
        open (unit = 105, file = 'Output/Impulse_1dImpulseVsEnergy.dat', status='replace', action = 'write')
        !note: integrated1 is a fucking bad name
        open (unit = 200, file = 'Output/Impulse_Integrated1.dat', status='replace', action = 'write')
        !open (unit = 201, file = 'Output/Impulse_Integrated2.dat', status='replace', action = 'write')
        !open (unit = 202, file = 'Output/Impulse_Integrated3.dat', status='replace', action = 'write')
    endif 

    !get the GR ray-tracing CONTINUUM parameters which are stored in the module gr_continuum
    if (nlp .eq. 1) then
       gso(1) = real( dgsofac(spin,h(1)) )
       call getlens(spin,h(1),mu0,lens(1),tauso(1),cosdelta_obs(1))
       if( tauso(1) .ne. tauso(1) ) stop "tauso is NaN"
    else
       !here the observed cutoffs are set from the temperature in the source frame   
       do m = 1, nlp
          gso(m) = real( dgsofac(spin,h(m)) )
          call getlens(spin,h(m),mu0,lens(m),tauso(m),cosdelta_obs(m))
          if( tauso(m) .ne. tauso(m) ) stop "tauso is NaN"
       enddo
    endif
    
    ! Set up observer's camera ( alpha = rn sin(phin), beta = mueff rn cos(phin) )
    ! to do full GR ray tracing with      
    mueff  = max( mu0 , 0.3d0 )
    rnmin  = rfunc(spin,mu0)
    !Grid to do in full GR
    call getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
    !Grid for Newtonian approximation
    call getrgrid(rnmax,rout,mueff,nron,nphin,rnn,domegan)

    ! Trace rays in full GR for the small camera (ie with relativistic effects) from the osberver to the disk,
    !which is why it doesnt depend on h
    if(status_re_tau) then !Only if the geodesics grid isn't loaded
        dotrace = .false.
        if( abs(spinsav-spin)  .gt. tiny(spin)   ) dotrace = .true.
        if( abs(musav-mu0)     .gt. tiny(mu0)    ) dotrace = .true.
        if( abs(routsav-rout)  .gt. tiny(rout)   ) dotrace = .true.
        if( abs(mudsav-mudisk) .gt. tiny(mudisk) ) dotrace = .true.         
        if( dotrace )then
            call GRtrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
            spinsav = spin
            musav   = mu0
            routsav = rout
            mudsav  = mudisk
        end if
    end if

    ! Set frequency array
    do fbin = 1,nf
        fi(fbin) = flo * (fhi/flo)**((float(fbin)-0.5d0)/dble(nf))
    end do
    if( fhi .lt. tiny(fhi) ) fi(1) = 0.0d0

    !initialize radius grid, angles, and transfer functions
    dlogr    = log10(rnmax/rin) / real(xe-1)
    cos0     = mu0
    sin0     = sqrt(1.0-cos0**2)
    frobs    = 0.0 !Initialised observer's reflection fraction

    ! Calculate dcos/dr and time lags vs r for the lamppost model
    call getdcos(spin,h,mudisk,ndelta,nlp,rout,npts,rlp,dcosdr,tlp,cosd,cosdout) 

    !set continuum normalisations depending on model flavour 
    if( dset .eq. 0 )then
        pnorm = 1.d0 / ( 4.d0 * pi )
    else
        pnorm = pnormer(b1,b2,qboost)  
    end if

    !loop over all photon directions (l), disk radii (i), disk azimuth (j), and calculate the contribution to the
    !transfer function/convolution kernel in energy (gbin), frequency (fbin), emission angle (mubin), disk radial 
    !bin (rbin) from the m-th/nl-th lamp post
    
    ! Construct the transfer function by summing over all pixels
    odisc    = 1       !flag to ensure the chosen disk radius is between rin and rout
    i        = nro + 1
    do while( odisc .eq. 1 .and. i .gt. 1 )             !main loops of the subroutine: first is for GR
        i = i - 1                                       !i counts over the camera until it reaches the disk inner radius
        odisc = 0
        do j = 1,NPHI                                   !azimuth over BH on the disk
            phin  = (j-0.5) * 2.d0 * pi / dble(nphi) 
            alpha = rn(i) * sin(phin)
            beta  = -rn(i) * cos(phin) * mueff
            !If the ray hits the disk, calculate flux and time lag
            if( pem1(j,i) .gt. 0.0d0 )then
                re    = re1(j,i)
                if( re .gt. rin .and. re .lt. rout )then
                    odisc = 1  
                    do m=1,nlp                           
                        taudo = taudo1(j,i)           
                        g = dlgfacthick(spin,mu0,alpha,re,mudisk) !disk to observer g factor
                        !Find the rlp bin that corresponds to re
                        kk = get_index(rlp(:,m),ndelta,re,rmin,npts(m))
                        !Interpolate (or extrapolate) the time function
                        tausd = interper(rlp(:,m),tlp(:,m),ndelta,re,kk)
                        tau(m) = (1.d0+zcos)*(tausd+taudo-tauso(1)) !Time lag between direct and reflected photons
                        !Interpolate |dcos\delta/dr| function                  
                        cosfac = interper(rlp(:,m),dcosdr(:,m),ndelta,re,kk)
                        mus = interper(rlp(:,m),cosd(:,m),ndelta,re,kk)
                        !Extrapolate to Newtonian if need be
                        if( kk .eq. npts(m) ) then
                            cosfac = newtex(rlp(:,m),dcosdr(:,m),ndelta,re,h(m),honr,kk)
                            mus = newtex(rlp(:,m),cosd(:,m),ndelta,re,h(m),honr,kk)
                        end if
                        !Calculate angular emissivity
                        ptf = pnorm * pfunc_raw(-mus,b1,b2,qboost)
                        !Calculate flux from pixel
                        gsd(m) = dglpfacthick(re,spin,h(m),mudisk) !source to disk g factor
                        emissivity(m) = gsd(m)**Gamma * 2.d0 * pi * ptf
                        emissivity(m) = emissivity(m) * cosfac / dareafac(re,spin)                  
                        dFe(m) = emissivity(m) * g**3 * domega(i) / (1.d0+zcos)**3
                        !calculate extra factors that go into the transfer functions for double lps
                        if (nlp .gt. 1) then
                            thetafac(m) = emissivity(m)*gso(m)**(Gamma-2.)*gsd(m)**(2.-Gamma)                      
                        else !single lamp post case, double check this later
                            thetafac(m) = 1.                            
                        endif                        
                        !Add to reflection fraction
                        frobs(m) = frobs(m) + 2.0*g**3*gsd(m)*cosfac/dareafac(re,spin)*domega(i)    
                        ! write(*,*) g, gsd(m), cosfac, dareafac(re, spin), domega(i)
                    end do
                   
                    !tbd: put a second for loop over lps here, now that both emissivity/dfe/tau are known
                    do nl=1,nlp 
                        !Work out energy bin
                        gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
                        gbin = MAX( 1    , gbin  )
                        gbin = MIN( gbin , ne    )         
                        !Work out radial bin
                        rbin = ceiling( log10(re/rin) / dlogr )
                        rbin = MAX( rbin , 1  )
                        rbin = MIN( rbin , xe )
                        !Add to the radial dependence of the transfer function TBD MAKE SURE THIS IS RIGHT
                        dfer_arr(rbin) = dfer_arr(rbin) + dFe(nl)                 
                        !Calculate emission angle and work out which mue bin to add to
                        mue   = demang(spin,mu0,re,alpha,beta)
                        mubin = ceiling( mue * dble(me) )
                        !calculate the extra factors for w2/3
                        !if (nl .eq. 1 .and. nlp .gt. 1) then
                        !    emisfac = emissivity(1)+eta_0*emissivity(2)                          
                        if (nlp .gt. 1) then
                            emisfac = (emissivity(1)+eta_0*emissivity(2))/(1.+eta_0)
                            kfac = (emissivity(1)+eta_0*emissivity(2))/(thetafac(1)+eta_0*thetafac(2)) 
                        else
                            emisfac = emissivity(1)
                            kfac = emissivity(1)
                            !single lamp post case, double check this later
                        endif     
                        !this is just to make the formatting below less ugly     
                        normfac = real(g**3*domega(i)/(1.d0+zcos)**3)                 
                        !Add to the transfer function integral
                        do fbin = 1,nf
                            cexp = cmplx(cos(real(2.d0*pi*tau(nl)*fi(fbin))),sin(real(2.d0*pi*tau(nl)*fi(fbin))))
                            ker_W0(nl,gbin,fbin,mubin,rbin) = ker_W0(nl,gbin,fbin,mubin,rbin) + real(dFe(nl))*cexp
                            ker_W1(nl,gbin,fbin,mubin,rbin) = ker_W1(nl,gbin,fbin,mubin,rbin) + &
                                                              real(log(gsd(nl)))*real(dFe(nl))*cexp  
                            !tbd redo these transfer functions                             
                            ker_W2(nl,gbin,fbin,mubin,rbin) = ker_W2(nl,gbin,fbin,mubin,rbin) + &
                                                              emisfac*normfac*cexp
                            ker_W3(nl,gbin,fbin,mubin,rbin) = ker_W3(nl,gbin,fbin,mubin,rbin) + &
                                                              kfac*thetafac(nl)*normfac*cexp 
                        end do
                        !if large verbose, start saving the impulse response function to file 
                        if( verbose .gt. 1 ) then
                            !find the appropriate energy and time bins
                            gbin = ceiling(g/dg) 
                            gbin = MAX( 1    , gbin  )
                            gbin = MIN( gbin , ne    )
                            tbin = ceiling( log10( tau(nl) / tar(0) ) / dlogt )
                            !write(102,*)re,tau,log10( tau / tar(0) ) / dlogt
                            tbin = MAX( 1    , tbin )
                            tbin = MIN( tbin , nt   )
                            ! kernel of the impulse response function              
                            resp(gbin,tbin) = resp(gbin,tbin) + dFe(nl)  
                        end if 
                    end do                    
                end if
            end if                
        end do
    end do

    ! Now trace rays for that bigger camera (obviously a lot easier because it's Newtonian)
    do i = 1,nron
        do j = 1,nphin
            phin  = (j-0.5) * 2.d0 * pi / dble(nphin) 
            alpha = rnn(i) * sin(phin)
            beta  = -rnn(i) * cos(phin) * mueff
            call drandphithick(alpha,beta,mu0,mudisk,re,phie)
            !If the ray hits the disk, calculate flux and time lag
            if( re .gt. rin .and. re .lt. rout )then
                do m=1,nlp
                    g = dlgfacthick( spin,mu0,alpha,re,mudisk )
                    !Find the rlp bin that corresponds to re
                    kk = get_index(rlp(:,m),ndelta,re,rmin,npts(m))
                    !Time lag
                    tau(m) = sqrt(re**2+(h(m)-honr*re)**2) - re*(sin0*sindisk*cos(phie)+mu0*mudisk ) + h(1)*mu0
                    tau(m) = (1.d0+zcos)*tau(m)
                    !Interpolate |dcos\delta/dr| function
                    cosfac = interper(rlp(:,m),dcosdr(:,m),ndelta,re,kk)
                    mus = interper(rlp(:,m),cosd(:,m),ndelta,re,kk)
                    !Extrapolate to Newtonian if needs be
                    if( kk .eq. npts(m) )then
                        cosfac = newtex(rlp(:,m),dcosdr(:,m),ndelta,re,h(m),honr,kk)
                        mus = newtex(rlp(:,m),cosd(:,m),ndelta,re,h(m),honr,kk)
                    end if
                    !Calculate angular emissivity
                    ptf = pnorm * pfunc_raw(-mus,b1,b2,qboost)
                    !Calculate flux from pixel
                    gsd(m) = dglpfacthick(re,spin,h(m),mudisk)
                    emissivity(m) = gsd(m)**Gamma * 2.d0 * pi * ptf
                    emissivity(m) = emissivity(m) * cosfac / dareafac(re,spin)
                    dFe(m) = emissivity(m) * g**3 * domegan(i) / (1.d0+zcos)**3
                    if (nlp .gt. 1) then
                        thetafac(m) = emissivity(m)*gso(m)**(Gamma-2.)*gsd(m)**(2.-Gamma)                      
                    else !single lamp post case, double check this later
                        thetafac(m) = 1.                      
                    endif
                    !Add to reflection fraction
                    frobs(m) = frobs(m) + 2.0*g**3*gsd(m)*cosfac/dareafac(re,spin)*domegan(i)
                end do 
                
                do nl=1,nlp
                    !Work out energy bin
                    gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
                    gbin = MAX( 1    , gbin  )
                    gbin = MIN( gbin , ne    )
                    !Work out radial bin
                    rbin = ceiling( log10(re/rin) / dlogr )
                    rbin = MAX( rbin , 1  )
                    rbin = MIN( rbin , xe )
                    !Add to the radial dependence of the transfer function
                    dfer_arr(rbin) = dfer_arr(rbin) + dFe(nl)                     
                    !Calculate emission angle and work out which mue bin to add to
                    mue = demang(spin,mu0,re,alpha,beta)
                    mubin = ceiling( mue * dble(me) )
                    !calculate the extra factors for w2/3
                    if (nlp .gt. 1) then
                        emisfac = (emissivity(1)+eta_0*emissivity(2))/(1.+eta_0)
                        kfac = (emissivity(1)+eta_0*emissivity(2))/(thetafac(1)+eta_0*thetafac(2)) 
                    else
                        emisfac = emissivity(1)
                        kfac = emissivity(1)
                        !single lamp post case, double check this later
                    endif  
                    !this is just to make the formatting below less ugly     
                    normfac = real(g**3*domegan(i)/(1.d0+zcos)**3)                 
                    !Add to the transfer function integral
                    do fbin = 1,nf
                        cexp = cmplx(cos(real(2.d0*pi*tau(nl)*fi(fbin))),sin(real(2.d0*pi*tau(nl)*fi(fbin))))
                        ker_W0(nl,gbin,fbin,mubin,rbin) = ker_W0(nl,gbin,fbin,mubin,rbin) + real(dFe(nl))*cexp
                        ker_W1(nl,gbin,fbin,mubin,rbin) = ker_W1(nl,gbin,fbin,mubin,rbin) + &
                                                          real(log(gsd(nl)))*real(dFe(nl))*cexp  
                        !tbd redo these transfer functions                             
                        ker_W2(nl,gbin,fbin,mubin,rbin) = ker_W2(nl,gbin,fbin,mubin,rbin) + &
                                                          emisfac*normfac*cexp
                        ker_W3(nl,gbin,fbin,mubin,rbin) = ker_W3(nl,gbin,fbin,mubin,rbin) + &
                                                          kfac*thetafac(nl)*normfac*cexp 
                    end do          
                    !if large verbose, start saving the impulse response function to file 
                    if( verbose .gt. 1 ) then
                        !find the appropriate energy and time bins
                        gbin = ceiling(g/dg) 
                        gbin = MAX( 1    , gbin  )
                        gbin = MIN( gbin , ne    )
                        tbin = ceiling( log10( tau(nl) / tar(0) ) / dlogt )
                        !write(102,*)re,tau,log10( tau / tar(0) ) / dlogt
                        tbin = MAX( 1    , tbin )
                        tbin = MIN( tbin , nt   )
                        ! kernel of the impulse response function              
                        resp(gbin,tbin) = resp(gbin,tbin) + dFe(nl)  
                    end if 
                end do
            end if
        end do
    end do
    
    do m=1,nlp 
        ! Calculate 4pi p(theta0,phi0) = ang_fac
        ang_fac = 4.d0 * pi * pnorm * pfunc_raw(-cosdelta_obs(m),b1,b2,qboost)
        ! Adjust the lensing factor (easiest way to keep track)
        lens(m) = lens(m) * ang_fac                 
        ! Calculate the relxill reflection fraction for one columncosdout
        frrel(m) = sysfref(rin,rlp(:,m),cosd(:,m),ndelta,cosdout(m))    
        !Finish calculation of observer's reflection fraction
        frobs(m) = frobs(m) / dgsofac(spin,h(m)) / lens(m)
    end do

    
    !TBD DOUBLE CHECK WTF IS BEING PRINTED TO FILE HERE I MEAN SERIOUSLY
    !finish saving the impulse response function to file
    if( verbose .gt. 1 ) then
        ! Deal with edge effects
        do tbin = 1,nt
            resp(1,tbin)  = 0.0
            resp(ne,tbin) = 0.0
        end do
        do gbin = 1,ne
            resp(gbin,1)  = 0.0
            resp(gbin,nt) = 0.0
        end do          
        do tbin = 1,nt
            sumresp = 0.0
            do gbin = 1,ne
                sumresp = sumresp + resp(gbin,tbin)
                E = gbin*dg  !10**( float(gbin-ne/2) * dloge )
                !write(104,*)0.5*(tar(tbin)+tar(tbin-1)),E,resp(gbin,tbin)
            end do
            write(103,*)0.5*(tar(tbin)+tar(tbin-1)),sumresp
            !write(201, *) sumresp
        end do

        do gbin = 1,ne
            sumresp = 0.0
            do tbin = 1,nt
                sumresp = sumresp + resp(gbin,tbin)
            end do
            write(105,*)gbin*dg,sumresp
            !write(202, *) sumresp
        end do
        
        do gbin = 1,ne
            write(200,*) resp(gbin, :)
        enddo
    end if    

    ! !calculate the ionization/density/gsd radial profiles 
    ! if( dset .eq. 0 .or. size(h) .eq. 2) then
    !     call radfunctions_dens(verbose,xe,rin,rnmax,eta_0,dble(rlxi),dble(lognep),spin,h,Gamma,honr,rlp&
    !                            &,dcosdr,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner,dfer_arr)
    ! else
    !     call radfuncs_dist(xe,rin,rnmax,b1,b2,qboost,fcons,&
    !                        & dble(lognep),spin,h(1),honr,rlp,dcosdr,cosd,ndelta,rmin,npts(1),&
    !                        & logxir,gsdr,logner,pnorm)
    ! end if
    ! !Outputs: logxir(1:xe),gsdr(1:xe), logner(1:xe)

    if (verbose .gt. 1) then
        !close(102)
        close(103)
        !close(104)
        close(105)
        close(200)
        !close(201)
        !close(202)
    endif 

    return
end subroutine rtrans
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
! Extrapolates using Newtonian value
  implicit none
  integer ndelta,kk
  double precision newtex,rlp(ndelta),dcosdr(ndelta),re,h,honr,cosfac
  newtex = dcosdr(kk) *  ( (h-honr*rlp(kk))**2 + rlp(kk)**2 )**1.5 / rlp(kk)
  newtex = newtex * re / ( (h-honr*re     )**2 + re**2      )**1.5
  return
end function newtex  
!-----------------------------------------------------------------------
