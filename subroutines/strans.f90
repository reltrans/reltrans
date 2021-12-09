!-----------------------------------------------------------------------
subroutine rtrans(verbose,dset,nlp,spin,h,mu0,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,lumratio,&
                  fcons,contx_int,tauso,lens,cosdelta_obs,nro,nphi,ne,dloge,nf,fhi,flo,me,xe,rlxi,lognep,&
                  transe,transe_a,frobs,frrel,logxir,gsdr,logner)
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
    ! rlxi                  Maximum value of logxi(r) in the disc
    ! OUTPUT
    ! transe(ne,nf,me,xe)   Transfer function as a function of energy, frequency, emission ngle and radius
    ! transe_a(ne,nf,me,xe) Second transfer function as a function of energy, frequency, emission ngle and radius
    !                            This is for the non-linear effects
    ! frobs                 Observer's reflection fraction
    ! frrel                 Reflection fraction defined by relxilllp
    ! lens                  Lensing factor for direct emission * 4pi p(theta0,phi0)
    ! logxir(xe),gsdr(xe)   logxi (ionization parameter) and gsd (source to disc blueshift) as a function of radius
    use dyn_gr
    use blcoordinate
    implicit none
    integer nro,nphi,ndelta,ne,nf,me,xe,dset,nlp
    double precision spin,h(nlp),mu0,Gamma,rin,rout,zcos,fhi,flo,honr
    double precision b1,b2,qboost
    double precision fcons,cosdout,logxir(xe),gsdr(xe),logner(xe),dfer_arr(xe)
    real rlxi, dloge, lognep
    complex cexp,transe(ne,nf,me,xe),transe_a(ne,nf,me,xe)
    integer i,npts(nlp),j,odisc,n,gbin,rbin,mubin,l,m,k
    parameter (ndelta=1000)
    double precision domega(nro),d,taudo,g,dlgfacthick,dFe,newtex
    double precision tauso(nlp),lens(nlp),cosdelta_obs(nlp),contx_int(nlp)
    double precision rlp_column(ndelta),dcosdr_column(ndelta),tlp_column(ndelta),cosd_column(ndelta)
    double precision rlp(ndelta,nlp),dcosdr(ndelta,nlp),tlp(ndelta,nlp),cosd(ndelta,nlp)
    double precision alpha,beta,cos0,sin0,phie,re,gsd
    double precision tau,tausd,emissivity,cosfac,dglpfacthick,dareafac
    integer kk,fbin,get_index
    double precision rmin,disco,rfunc,scal,velocity(3),mudisk,sysfref
    double precision rnmax,rnmin,rn(nro),phin,mueff,dlogr,interper
    double precision fi(nf),dgsofac,sindisk,mue,demang,frobs,cosdin,frrel
    double precision pnorm,pnormer,mus,ptf,pfunc_raw,ang_fac
    integer nron,nphin,nrosav,nphisav,verbose
    double precision spinsav,musav,routsav,mudsav,rnn(nro),domegan(nro)
    integer myenv
    double precision lximax
    double precision lumratio 
    logical dotrace
    
    !arrays to save the transfer function
    integer, parameter :: nt = 2**9
    integer            :: tbin
    double precision   :: tmin, tmax, sumresp, tar(0:nt), dlogt, dg, E
    double precision, allocatable :: resp(:,:)

    data nrosav,nphisav,spinsav,musav /0,0,2.d0,2.d0/
    save nrosav,nphisav,spinsav,musav,routsav,mudsav
       
    ! Settings 
    nron     = 100
    nphin    = 100
    rmin     = disco( spin )
    scal     = 1.d0
    velocity = 0.d0
    mudisk   = honr / sqrt( honr**2 + 1.d0  )
    sindisk  = sqrt( 1.d0 - mudisk**2 )
    dfer_arr = 0.
    
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
        !open (unit = 102, file = 'Output/Impulse_ReTau.dat', status='replace', action = 'write')
        !open (unit = 104, file = 'Output/Impulse_2dImpulse.dat', status='replace', action = 'write')
        open (unit = 103, file = 'Output/Impulse_1dImpulseVsTime.dat', status='replace', action = 'write')
        open (unit = 105, file = 'Output/Impulse_1dImpulseVsEnergy.dat', status='replace', action = 'write')
        !note: integrated1 is a fucking bad name
        open (unit = 200, file = 'Output/Impulse_Integrated1.dat', status='replace', action = 'write')
        !open (unit = 201, file = 'Output/Impulse_Integrated2.dat', status='replace', action = 'write')
        !open (unit = 202, file = 'Output/Impulse_Integrated3.dat', status='replace', action = 'write')
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
    transe   = 0.0 !Initialised transfer function
    transe_a = 0.0 !Initialised transfer function
    frobs    = 0.0 !Initialised observer's reflection fraction

    ! Calculate dcos/dr and time lags vs r for the lamppost model
    call getdcos(spin,h,mudisk,ndelta,nlp,rout,npts,rlp,dcosdr,tlp,cosd,cosdout) 

    !set continuum normalisations depending on model flavour 
    if( dset .eq. 0 )then
        pnorm = 1.d0 / ( 4.d0 * pi )
    else
        pnorm = pnormer(b1,b2,qboost)  
    end if

    !loop over all lmaposts (m), photon directions (l), disk radii (i), disk azimuth (j), and calculate the contribution to the
    !transfer function/convolution kernel in energy (gbin), frequency (fbin), emission angle (mubin), disk radialb in (rbin)
    do m=1,nlp
        !set normalisation of second LP wrt to first
        if (m .gt. 1) pnorm = lumratio*pnorm
        !get appropriate arrays for rlp/tlp/dcosdr/cosd                 
        do l=1,ndelta
            rlp_column(l)=rlp(l,m)
            tlp_column(l)=tlp(l,m)
            dcosdr_column(l)=dcosdr(l,m)
            cosd_column(l)=cosd(l,m)    
        end do     
        ! Calculate 4pi p(theta0,phi0) = ang_fac
        ang_fac = 4.d0 * pi * pnorm * pfunc_raw(-cosdelta_obs(m),b1,b2,qboost)
        !to do luminosity ratio between LPs: include factor in front of pnorm here.
        ! Adjust the lensing factor (easiest way to keep track)
        lens(m) = lens(m) * ang_fac                 
        ! Calculate the relxill reflection fraction for one column...need to do this for both
        frrel = sysfref(rin,rlp_column,cosd_column,ndelta,cosdout)    
        ! Construct the transfer function by summing over all pixels
        odisc    = 1                                        !flag to ensure the chosen disk radius is between rin and rout
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
                        taudo = taudo1(j,i)           
                        g     = dlgfacthick( spin,mu0,alpha,re,mudisk )     !this is disk to observer g factor                                               
                        !Find the rlp bin that corresponds to re
                        kk = get_index(rlp_column,ndelta,re,rmin,npts(m))
                        !Interpolate (or extrapolate) the time function
                        tausd = interper(rlp_column,tlp_column,ndelta,re,kk)
                        tau   = (1.d0+zcos) * (tausd+taudo-tauso(m))           !This is the time lag between direct and reflected photons
                        !Interpolate |dcos\delta/dr| function                  
                        cosfac = interper(rlp_column,dcosdr_column,ndelta,re,kk)
                        mus    = interper(rlp_column,cosd_column,ndelta,re,kk)
                        !Extrapolate to Newtonian if needs be
                        if( kk .eq. npts(m) )then
                            cosfac = newtex(rlp_column,dcosdr_column,ndelta,re,h(m),honr,kk)
                            mus    = newtex(rlp_column,cosd_column,ndelta,re,h(m),honr,kk)
                        end if
                        !Calculate angular emissivity
                        ptf        = pnorm * pfunc_raw(-mus,b1,b2,qboost)
                        !Calculate flux from pixel
                        gsd        = dglpfacthick(re,spin,h(m),mudisk)  !source to disk g factor
                        emissivity = gsd**Gamma * 2.d0 * pi * ptf
                        emissivity = emissivity * cosfac / dareafac(re,spin)                  
                        dFe        = emissivity * g**3 * domega(i) / (1.d0+zcos)**3
                        !Add to reflection fraction
                        frobs      = frobs + 2.0 * g**3 * gsd * cosfac/dareafac(re,spin) * domega(i)
                        !Work out energy bin
                        gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
                        gbin = MAX( 1    , gbin  )
                        gbin = MIN( gbin , ne    )                        
                        !Work out radial bin
                        rbin = ceiling( log10(re/rin) / dlogr )
                        rbin = MAX( rbin , 1  )
                        rbin = MIN( rbin , xe )
                        !Add to the radial dependence of the transfer function
                        dfer_arr(rbin) = dfer_arr(rbin) + dFe                 
                        !Calculate emission angle and work out which mue bin to add to
                        mue   = demang(spin,mu0,re,alpha,beta)
                        mubin = ceiling( mue * dble(me) )
                        !Add to the transfer function integral
                        do fbin = 1,nf
                            cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                            transe(gbin,fbin,mubin,rbin)   = transe(gbin,fbin,mubin,rbin)                  + real(dFe) * cexp
                            transe_a(gbin,fbin,mubin,rbin) = transe_a(gbin,fbin,mubin,rbin) + real(log(g)) * real(dFe) * cexp
                        end do
                        !if large verbose, start saving the impulse response function to file 
                        if( verbose .gt. 1 ) then
                            !find the appropriate energy and time bins
                            gbin = ceiling(g/dg) 
                            gbin = MAX( 1    , gbin  )
                            gbin = MIN( gbin , ne    )
                            tbin = ceiling( log10( tau / tar(0) ) / dlogt )
                            write(102,*)re,tau,log10( tau / tar(0) ) / dlogt
                            tbin = MAX( 1    , tbin )
                            tbin = MIN( tbin , nt   )
                            ! kernel of the impulse response function              
                            resp(gbin,tbin) = resp(gbin,tbin) + dFe  
                        end if
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
                    g = dlgfacthick( spin,mu0,alpha,re,mudisk )
                    !Find the rlp bin that corresponds to re
                    kk = get_index(rlp_column,ndelta,re,rmin,npts(m))
                    !Time lag
                    tau = sqrt(re**2+(h(m)-honr*re)**2) - re*(sin0*sindisk*cos(phie)+mu0*mudisk ) + h(m)*mu0
                    tau = (1.d0+zcos) * tau
                    !Interpolate |dcos\delta/dr| function
                    cosfac = interper(rlp_column,dcosdr_column,ndelta,re,kk)
                    mus    = interper(rlp_column,cosd_column,ndelta,re,kk)
                    !Extrapolate to Newtonian if needs be
                    if( kk .eq. npts(m) )then
                        cosfac = newtex(rlp_column,dcosdr_column,ndelta,re,h(m),honr,kk)
                        mus    = newtex(rlp_column,cosd_column,ndelta,re,h(m),honr,kk)
                    end if
                    !Calculate angular emissivity
                    ptf        = pnorm * pfunc_raw(-mus,b1,b2,qboost)
                    !Calculate flux from pixel
                    gsd        = dglpfacthick(re,spin,h(m),mudisk)
                    emissivity = gsd**Gamma * 2.d0 * pi * ptf
                    emissivity = emissivity * cosfac / dareafac(re,spin)
                    dFe        = emissivity * g**3 * domegan(i) / (1.d0+zcos)**3
                    !Add to reflection fraction
                    frobs      = frobs + 2.0 * g**3 * gsd * cosfac/dareafac(re,spin) * domega(i)
                    !Work out energy bin
                    gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
                    gbin = MAX( 1    , gbin  )
                    gbin = MIN( gbin , ne    )
                    !Work out radial bin
                    rbin = ceiling( log10(re/rin) / dlogr )
                    rbin = MAX( rbin , 1  )
                    rbin = MIN( rbin , xe )
                    !Add to the radial dependence of the transfer function
                    dfer_arr(rbin) = dfer_arr(rbin) + dFe                     
                    !Calculate emission angle and work out which mue bin to add to
                    mue   = demang(spin,mu0,re,alpha,beta)
                    mubin = ceiling( mue * dble(me) )
                    !Add to the transfer function integral
                    do fbin = 1,nf
                        cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                        transe(gbin,fbin,mubin,rbin)   = transe(gbin,fbin,mubin,rbin)                  + real(dFe) * cexp
                        transe_a(gbin,fbin,mubin,rbin) = transe_a(gbin,fbin,mubin,rbin) + real(log(g)) * real(dFe) * cexp
                    end do
                    !if large verbose, start saving the impulse response function to file 
                    if( verbose .gt. 1 ) then
                        !find the appropriate energy and time bins
                        gbin = ceiling(g/dg) 
                        gbin = MAX( 1    , gbin  )
                        gbin = MIN( gbin , ne    )
                        tbin = ceiling( log10( tau / tar(0) ) / dlogt )
                        write(102,*)re,tau,log10( tau / tar(0) ) / dlogt
                        tbin = MAX( 1    , tbin )
                        tbin = MIN( tbin , nt   )
                        ! kernel of the impulse response function              
                        resp(gbin,tbin) = resp(gbin,tbin) + dFe  
                    end if
                end if
            end do
        end do
        !Finish calculation of observer's reflection fraction - FIGURE OUT A WAY TO DO THIS FOR ALL LPS
        frobs = frobs / dgsofac(spin,h(m)) / lens(m)
    end do 

    
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

    !calculate the ionization/density/gsd radial profiles 
    if( dset .eq. 0 .or. size(h) .eq. 2) then
        call radfunctions_dens(verbose,xe,rin,rnmax,lumratio,dble(rlxi),dble(lognep),spin,h,Gamma,honr,rlp&
                               &,dcosdr,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner,dfer_arr)
    else
        call radfuncs_dist(xe,rin,rnmax,b1,b2,qboost,fcons,&
                           & dble(lognep),spin,h(1),honr,rlp,dcosdr,cosd,ndelta,rmin,npts(1),&
                           & logxir,gsdr,logner,pnorm)
    end if
    !Outputs: logxir(1:xe),gsdr(1:xe), logner(1:xe)

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
function sysfref(rin,rlp,cosd,ndelta,cosdout)
! Calculates the relxill definition of reflection fraction        
! In: rin,rlp,ndelta,cosd,cosdout
! Out: stsfref
  implicit none
  integer ndelta
  double precision sysfref,rin,rlp(ndelta),cosd(ndelta),cosdout
  integer n
  double precision cosdin
  n = 1
  do while( rin .gt. rlp(n) )
    n = n + 1
  end do
  !Inperpolate
  if( n .eq. 1 )then
    !Need to extrapolate
    cosdin = (cosd(n+1)-cosd(n))*(rin-rlp(n))/(rlp(n+1)-rlp(n))
    cosdin = cosdin + cosd(n)
  else
    !Can actually interpolate
    cosdin = (cosd(n)-cosd(n-1))*(rin-rlp(n-1))/(rlp(n)-rlp(n-1))
    cosdin = cosdin + cosd(n-1)
  end if
  sysfref = ( cosdin - cosdout ) / ( 1.0 + cosdout )
  return
end function sysfref  
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine radfuncs_dist(xe, rin, rnmax, b1, b2, qboost, fcons,&
     & lognep, spin, h, honr, rlp, dcosdr, cosd, ndelta, rmin, npts,&
     & logxieff, gsdr, logner, pnorm)
! Calculates ionization parameter, g-factor and disc density as a
! function of disc radius. This subroutine is only called for rtdist
!
! In  : xe,rin,rnmax,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
! Out :
! logxir(1:xe) -- Effective ionisation as a function of r
! gsdr(1:xe)   -- Source-to-disc blueshift as a function of r
! logner(1:xe) -- Log10 of electron density as a function of r
!note: this does not work with multiple lamposts for now
  implicit none
  integer         , intent(IN)   :: xe, ndelta, npts 
  double precision, intent(IN)   :: rin, rmin, rnmax, b1, b2, qboost
  double precision, intent(IN)   :: fcons, lognep, spin, h, honr
  double precision, intent(IN)   :: rlp(ndelta), dcosdr(ndelta), cosd(ndelta)
  double precision, intent(INOUT):: logxieff(xe), gsdr(xe), logner(xe)
  integer          :: i, kk, get_index, myenv, adensity, verbose
  double precision :: pnorm,re,re1(xe),zA_logne,cosfac,mus,interper,newtex,mudisk
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: ptf,pfunc_raw,gsd,dglpfacthick,eps_bol,Fx(xe),logxir(xe),mui,dinang
  double precision :: pnormer,dareafac,lximax
! Decide on zone a density profile or constant density profile
  adensity = myenv("A_DENSITY",1)
  adensity = min( adensity , 1 )
  adensity = max( adensity , 0 )
! Set disk opening angle
  mudisk   = honr / sqrt( honr**2 + 1.d0  )
! Now loop through xe radial bins
  do i = 1,xe
     !Radius
     re     = (rnmax/rin)**(real(i-1) / real(xe))
     re     = re + (rnmax/rin)**(real(i) / real(xe))
     re     = re * rin * 0.5
     re1(i) = re
     !Density
     logner(i) = zA_logne(re,rin,lognep)     
     !Interpolate functions from rpl grid
     kk     = get_index(rlp,ndelta,re,rmin,npts)
     cosfac = interper(rlp,dcosdr,ndelta,re,kk)
     mus    = interper(rlp,cosd  ,ndelta,re,kk)
     if( kk .eq. npts )then !Extrapolation onto Newtonian grid
        cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
        mus    = newtex(rlp,cosd  ,ndelta,re,h,honr,kk)
     end if
     !Calculate 13.6 eV - 13.6 keV illuminating flux
     ptf     = pnorm * pfunc_raw(-mus,b1,b2,qboost)
     gsd     = dglpfacthick(re,spin,h,mudisk)
     !gsd     = dglpfac(re,spin,h)
     gsdr(i) = gsd
     eps_bol = gsd**2 * 2.0 * pi * ptf
     eps_bol = eps_bol * cosfac / dareafac(re,spin)
     Fx(i)   = fcons * 4.0 * pi * eps_bol
     !Calculate logxi(r)
     logxir(i) = log10( 4.0 * pi * Fx(i) ) - logner(i)
     !Now adjust to effective ionization parameter
     mui         = dinang(spin, re, h, mus)
     logxieff(i) = logxir(i) - 0.1505 - log10(mui)     
!     write(188,*)re,logxir(i),logxieff(i)
     
  end do
!  write(188,*)"no no"
  
!check max and min for both ionisation and density
  logxieff = max( logxieff , 0.d0  )
  logxieff = min( logxieff , 4.7d0 )
  ! logner   = max( logner , 15.d0  )
  ! logner   = min( logner , 22.d0 )
  !...no need to enforce limits on logne since this is done in myreflect()
  !This is needed because reflionx has a different maximum to xillverDCp

  verbose = myenv("REV_VERB",0)
  if( verbose .gt. 2 )then
     !Write out logxir for plots
     lximax = -huge(lximax)
     do i = 1,xe
        write(188,*)re1(i),Fx(i),logxir(i)
        lximax = max( lximax , logxieff(i) )
     end do
     write(188,*)"no no"
     write(*,*)"MAX LOGXIeff = ",lximax
  end if
  
  return
end subroutine radfuncs_dist  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
function pnormer(b1,b2,boost)
  implicit none
  double precision pnormer,b1,b2,boost,pi
  integer i,imax
  parameter (imax=1000)
  double precision integral,mu,pfunc_raw
  pi = acos(-1.d0)
  integral = 0.0
  do i = 1,imax
     mu       = -1.0 + 2.0*dble(i-1)/dble(imax)
     integral = integral + pfunc_raw(mu,b1,b2,boost)
  end do
  integral = integral * 2.0 / dble(imax)
  pnormer  = 1.0 / ( 2.0*pi*integral )
  return
end function pnormer
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function pfunc_raw(mu,b1,b2,boost)
  implicit none
  double precision pfunc_raw,mu,b1,b2,boost
  double precision calB,norm,pm,mup,p
  calB = 1.0/boost
  if( mu .le. 0.0 )then
     norm = 1.0
  else
     norm = calB
  end if
  pm = mu/sign(mu,1.d0)
  if( mu .eq. 0.d0 ) pm = 1.0
  mup = pm * ( norm**2 * (1.0/mu**2-1.0) + 1.0 )**(-0.5)
  p   = 1.0 + (b1+abs(b2))*abs(mup) + b2*mup**2
  p   = p * sqrt( 1.0 + mup**2 * (norm**2-1.0) )
  pfunc_raw = p  
  return
end function pfunc_raw  
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine radfunctions_dens(verbose,xe,rin,rnmax,lumratio,logxip,lognep,spin,h,Gamma,honr,rlp,dcosdr&
     &,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner,dfer_arr)
    ! In  : xe,rin,rnmax,lumratio,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
    ! Out : logxir(1:xe), gsdr(1:xe), logner(1:xe)
    implicit none
    integer         , intent(IN)   :: xe, ndelta, nlp, npts(nlp)
    double precision, intent(IN)   :: rin,rmin,rnmax,lumratio,logxip,lognep,spin,h(nlp),honr,Gamma,dfer_arr(xe)
    real                           :: gso(nlp)
    double precision, intent(IN)   :: rlp(ndelta,nlp), dcosdr(ndelta,nlp), cosd(ndelta,nlp), contx_int(nlp)
    double precision :: rlp_column(ndelta),dcosdr_column(ndelta),cosd_column(ndelta), dgsofac
    double precision, intent(INOUT):: logxir(xe), gsdr(xe), logner(xe)
    integer          :: i, kk, get_index, myenv, adensity, l, m, verbose
    double precision :: rp, logxinorm, lognenorm,  mus, interper, newtex, mui, dinang, gsd(nlp), dglpfacthick
    double precision :: xi_lp(xe,nlp), logxi_lp(xe,nlp), logxip_lp(nlp), xitot, xiraw, mylogne, mudisk, gsd_temp
    double precision, allocatable :: rad(:)

    !Decide on zone a density profile or constant density profile
    adensity = myenv("A_DENSITY",1)
    adensity = min( adensity , 1 )
    adensity = max( adensity , 0 )
    ! Set disk opening angle
    mudisk   = honr / sqrt( honr**2 + 1.d0  )
    
    allocate(rad(xe))
    !Now calculate logxi itself
    ! The loop calculates the raw xi and raw n_e.
    ! This means they are without normalization: only to find the maximum and the minimum. Remember that the max of the ionisation is not the same as the minumim in the density because the flux depends on r
    !The loops calculates also the correction factor mui

    !TBD: include luminosity ratio between LPs 
    do i = 1, xe        
        rad(i) = (rnmax/rin)**(real(i-1) / real(xe))
        rad(i) = rad(i) + (rnmax/rin)**(real(i) / real(xe))
        rad(i) = rad(i) * rin * 0.5
        !Initialize total ionization tracker
        xitot = 0. 
        gsd_temp = 0.
        !Now calculate the raw density (this matters only for high dens model reltransD)
        logner(i) = adensity * mylogne(rad(i), rin)
        do m=1,nlp
            do l=1,ndelta
                rlp_column(l) = rlp(l,m)
                dcosdr_column(l) = dcosdr(l,m)
                cosd_column(l) = cosd(l,m)
            end do    
            gso(m) = real( dgsofac(spin,h(m)) )     
            xi_lp(i,m) = xiraw(rad(i),spin,h(m),honr,rlp_column,dcosdr_column,ndelta,rmin,npts(m),mudisk,gsd(m))            
            if (m .eq. 2) xi_lp(i,m) = lumratio*xi_lp(i,m)
            !Calculate the incident angle for this bin
            kk = get_index(rlp_column, ndelta, rad(i), rmin, npts(m))
            mus = interper(rlp_column, cosd_column, ndelta, rad(i), kk)
            if( kk .eq. npts(m) ) mus = newtex(rlp_column, cosd_column, ndelta, rad(i), h(m), honr, kk)
            mui = dinang(spin, rad(i), h(m), mus)
            !Correction to account for the radial dependence of incident angle, and for the g factors
            xi_lp(i,m) = xi_lp(i,m)/(sqrt(2.)*mui)*contx_int(m)*(gso(m))**(Gamma-2)  
            xitot = xitot + xi_lp(i,m)
            gsd_temp = gsd_temp + gsd(m)*xi_lp(i,m)
        end do 
        !This and the line above calculate the gsd factor along the disk, averaging over the flux the disk sees from each LP 
        gsdr(i) = gsd_temp/xitot
        logxir(i) = log10(xitot) - logner(i)     
    end do
    !After the loop calculate the max and the min - ionization renormalized wrt to the first LP
    logxinorm = maxval(logxir)
    lognenorm = minval(logner)
    logxir = logxir - (logxinorm - logxip) 
    logner = logner - (lognenorm - lognep)    
    
    do m=1,nlp 
        do i=1,xe
            logxi_lp(i,m) = log10(xi_lp(i,m)) - logner(i) - lognenorm - logxinorm + lognep + logxip        
        end do
        logxip_lp(m) = max(maxval(logxi_lp(:,m)),0.)
    end do
    
    !Write radii, ionisation (for both and each LP), gamma factors, and log(xi(r))+log(ne(r)) (which is nearly the same as
    !epsilon(r) for identical coronal spectrra and gamma=2) to file. 
    !note 1) we need to do this before the ionisation array is set to have a minimum of 0, in order
    !to recover the correct scaling of the emissivity at large radii
    !2) in order to correctly compare the dfer_arr array with the single LP case, it has to be renormalized by (1+lumratio)
    if( verbose .gt. 1 ) then
        print*, "Peak ionisations from each LP: first " , logxip_lp(1), " second ", logxip_lp(2)
        open (unit = 27, file = 'Output/RadialScalings.dat', status='replace', action = 'write')
            do i = 1, xe
                write(27,*) rad(i), logxir(i), gsdr(i), logxir(i)+logner(i), logxi_lp(i,1), logxi_lp(i,2), dfer_arr(i) 
            end do 
        close(27)    
    end if
    
    !check max and min for ionisation 
    logxir = max( logxir , 0.d0  )
    logxir = min( logxir , 4.7d0 )
    
    deallocate(rad)

    return
end subroutine radfunctions_dens
!-----------------------------------------------------------------------

  
!-----------------------------------------------------------------------
function logxiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,mudisk,gsd)
! In: re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts
! Out: logxiraw,gsd
  implicit none
  integer ndelta,npts,kk,get_index
  double precision re,spin,h,honr,rlp(ndelta),dcosdr(ndelta),rmin,gsd
  double precision cosfac,interper,logxiraw,dareafac,dglpfacthick,newtex
  double precision mudisk
  !Calculate source to disc blueshift at this radius
  gsd = dglpfacthick(re,spin,h,mudisk)
  !gsd = dglpfac(re,spin,h)
  !Find the rlp bin that corresponds to re
  kk = get_index(rlp,ndelta,re,rmin,npts)
  !Interpolate to get |d\cos\delta/dr| at r=re
  cosfac = interper(rlp,dcosdr,ndelta,re,kk)
  !Extrapolate to Newtonian if needs be
  if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
  !Now can do the calculation
  logxiraw = log10( gsd**2 * cosfac / dareafac(re,spin) )
  return
end function logxiraw  
!-----------------------------------------------------------------------

function xiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,mudisk,gsd)
    ! In: re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts
    ! Out: logxiraw,gsd
    implicit none
    integer ndelta,npts,kk,get_index
    double precision re,spin,h,honr,rlp(ndelta),dcosdr(ndelta),rmin,gsd
    double precision cosfac,interper,xiraw,dareafac,dglpfacthick,newtex
    double precision mudisk
    !Calculate source to disc blueshift at this radius
    gsd = dglpfacthick(re,spin,h,mudisk)
    !Find the rlp bin that corresponds to re
    kk = get_index(rlp,ndelta,re,rmin,npts)
    !Interpolate to get |d\cos\delta/dr| at r=re
    cosfac = interper(rlp,dcosdr,ndelta,re,kk)
    !Extrapolate to Newtonian if needs be
    if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
    !Now can do the calculation
    xiraw = gsd**2 * cosfac / dareafac(re,spin) 
    return
end function xiraw  
!-----------------------------------------------------------------------
  

!-----------------------------------------------------------------------
function zA_logne(r,rin,lognep)
! log(ne), where ne is the density.
! This function is normalised to have a maximum of lognep.
! We can therefore calculate the ionization parameter by taking
! 4 \pi * Fx / nemin * zA_one_on_ne
  implicit none
  double precision zA_logne,r,rin,lognep,rp
  rp       = 25./9. * rin
  zA_logne = lognep + 1.5*log10(r/rp)
  zA_logne = zA_logne + 2.0*log10( 1.0 - sqrt( rin / rp ) )
  zA_logne = zA_logne - 2.0*log10( 1.0 - sqrt( rin / r  ) )
  return
end function zA_logne
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function mylogne(r,rin)
! Calculates log10(ne). Don't let r = rin
  implicit none
  double precision mylogne,r,rin
  mylogne = 1.5 * log10(r) - 2.0 * log10( 1.0 - sqrt(rin/r) )
  return
end function mylogne 
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function get_index(rlp,ndelta,re,rmin,npts)
  implicit none
  integer get_index,ndelta,npts,kk
  double precision rlp(ndelta),re,rmin
  kk = 2
  do while( ( rlp(kk) .le. re .or. rlp(kk-1) .lt. rmin ) .and. kk .lt. npts )
     kk = kk + 1
  end do
  get_index = kk
  return
end function get_index
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function interper(rlp,ylp,ndelta,re,kk)
! Interpolates the array ylp between ylp(kk-1) and ylp(kk)
! ylp is a function of rlp and rlp(kk-1) .le. re .le. rlp(kk)
  implicit none
  integer ndelta,kk
  double precision interper,rlp(ndelta),ylp(ndelta),re
  interper = (ylp(kk)-ylp(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
  interper = interper + ylp(kk)
  return
end function interper
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
