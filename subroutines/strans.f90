!-----------------------------------------------------------------------
      subroutine strans(spin,h,mu0,Gamma,rin,rout,honr,zcos,nro,nphi,ndelta,ne,dloge,&
           ear,nf,fhi,flo,mex,gex,xex,me,ge,xe,rlxi,sdmin,sdmax,ximin,ximax,transe,frobs,frrel,xbinhi,lens)
! Code to calculate the transfer function for an accretion disk.
! This code first does full GR ray tracing for a camera with impact parameters < bmax
! It then also does straight line ray tracing for impact parameters >bmax
! It adds both up to produce a transfer function for a disk extending from rin to rout
! INPUT
! spin,h,mu0,Gamma      Physical parameters (spin, source height, cos(inclination), photon index)
! rin,rout,honr         Physical parameters (disk inner radius, outer radius & scaleheight)
! zcos                  Cosmological redshift
! nro,nphi              Number of pixels on the observer's camera (b and phib)
! ndelta                Number of \delta values considered
! ne, dloge             Number of energy bins and maximum energy (compatible with FFT convolution)
! ear(0:ne)             Energy grid
! nf,fhi,flo            nf = Number of logarithmic frequency bins averaged over, range= flo to fhi
! mex
! gex
! xex
! rlxi
! OUTPUT
! sdmin,sdmax
! ximin ximax           Minimum and maximum values of logxi_eff (depends on emissivity, density and incidence angle)
! transe(ne,mex,gex,xex)    Transfer function as a function of energy and emission angle for the given frequency range
! frobs                 Observer's reflection fraction
! frrel                 Reflection fraction defined by relxilllp
! xbinhi                Highest xi bin that has an entry 
! lens                  Lensing factor for direct emission
        use blcoordinate
      implicit none
      integer nro,nphi,ndelta,ne,nf,mex,gex,sdbin,xex,me,ge,xe
      double precision spin,h,mu0,Gamma,rin,rout,zcos,fhi,flo,honr,cosdout
      real rlxi,ximin,ximax
      real ear(0:ne),dloge,sdmin,sdmax
      complex cexp,transe(ne,mex,gex,xex)
      integer i,npts,j,k,l,odisc,jj,nmax,n,xbin
      parameter (nmax=1000)
      double precision domega(nro),d
      double precision tauso,rlp(ndelta),dcosdr(ndelta),tlp(ndelta),cosd(ndelta)
      double precision alpha,beta,cos0,sin0,phi0,phie,re,gsd
      double precision taudo,dtaudo,g,dlgfac,dFe,lens
      double precision tau,tausd,emissivity,cosfac,dglpfac,dareafac,line
      integer gbin,kk,fbin,ebin,ifl,stat,myenv
      double precision f,rmin,disco,rfunc
      double precision amp,phase,sum,scal,velocity(3),f1234(4),lambda,q
      double precision pem,mudisk,mucros,sigmacros
      double precision rnmax,rnmin,rn(nro),phin,mueff
      double precision fi(nf),dgsofac,sindisk,mue,demang,frobs,cosdin,frrel
      double precision pem1(nmax,nmax),re1(nmax,nmax),taudo1(nmax,nmax)
      integer nron,nphin,nrosav,nphisav,mubin,xbinhi,adensity
      double precision spinsav,musav,routsav,mudsav,t1,t0,rnn(nro),domegan(nro)
      double precision mus,mui,dinang,xir,logxieff,logxir,xinorm,logxip,logxihi
      logical dotrace
      character (len=1) A_DENSITY
      data nrosav,nphisav,spinsav,musav /0,0,2.d0,2.d0/
      save nrosav,nphisav,spinsav,musav,routsav,mudsav,pem1,taudo1,re1
      
! Settings 
      nron     = 100
      nphin    = 100
      ifl      = 1
      rmin     = disco( spin )
      scal     = 1.d0
      velocity = 0.d0
      mudisk   = honr / sqrt( honr**2 + 1.d0  )
      sindisk  = sqrt( 1.d0 - mudisk**2 )
      
! Set up observer's camera ( alpha = rn sin(phin), beta = mueff rn cos(phin) )
! to do full GR ray tracing with      
      mueff  = max( mu0 , 0.3d0 )
      rnmin  = rfunc(spin,mu0)
      rnmax  = 300d0
      !Grid to do in full GR
      call getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
      !Grid for Newtonian approximation
      call getrgrid(rnmax,rout,mueff,nron,nphin,rnn,domegan)
      
! Set frequency array
      do fbin = 1,nf
        fi(fbin) = flo * (fhi/flo)**((float(fbin)-0.5d0)/dble(nf))
      end do
      if( fhi .lt. 1d-10 ) fi(1) = 0.0d0
      
! Set sensible distance for observer from the BH
      d = max( 1.0d4 , 2.0d2 * rnmax**2 )

! Set up g_{sd} array
      sdmax = real( dglpfac(rin ,spin,h) )
      sdmin = real( dglpfac(rout,spin,h) )

! Calculate lensing factor "lens" and source to observer time lag "tauso"
      call getlens(spin,h,mu0,lens,tauso)
      if( tauso .ne. tauso ) stop "tauso is NaN"
      
! Calculate dcos/dr and time lags vs r for the lamppost model
      call getdcos(spin,h,mudisk,ndelta,rout,npts,rlp,dcosdr,tlp,cosd,cosdout)

! Calculate cosdin in order to calculate the relxill reflection fraction
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
      frrel = ( cosdin - cosdout ) / ( 1.0 + cosdout )
      
! Now trace rays in full GR down to mu=0 for impact parameters with b<bmax
      cos0  = mu0
      sin0  = sqrt(1.0-cos0**2)
      transe  = 0.0 !Initialised transfer function
      frobs   = 0.0
      
! Trace rays in full GR for the small camera
! to convert alpha and beta to r and tau_do (don't care about phi)
! First work out if we even need to call this
      dotrace = .false.
      if( nro .ne. nrosav ) dotrace = .true.
      if( nphi .ne. nphisav ) dotrace = .true.
      if( abs(spinsav-spin) .gt. 1d-6 ) dotrace = .true.
      if( abs(musav-mu0) .gt. 1d-6 ) dotrace = .true.
      if( abs(routsav-rout) .gt. 1d-6 ) dotrace = .true.
      if( abs(mudsav-mudisk) .gt. 1d-6 ) dotrace = .true.
      if( dotrace )then
        call GRtrace(nmax,nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d,pem1,taudo1,re1)
        nrosav  = nro
        nphisav = nphi
        spinsav = spin
        musav   = mu0
        routsav = rout
        mudsav  = mudisk
      end if
      !Could replace this with reading in from a grid   ***Could paralellize the GRtrace subroutine***

! Decide on zone a density profile or constant density profile
      adensity = myenv("A_DENSITY",1)
      
! Set up logxi part of the calculation
      
      logxip = dble(rlxi)
      ximin   = 0.0
      ximax   = 4.7
      if( adensity .eq. 1 )then
        re      = (11.d0/9.d0)**2 * rin  !Max \xi for Fx~r^{-3} and n=zone a SS73 (eqn 2.11)
      else
        re = rin
      end if
! => \xi \propto r^{-9/2} [ 1 - sqrt(rin/r) ]^2 => rp = (11/9)**2 rin
      kk      = 2
      do while( ( rlp(kk) .le. re .or. rlp(kk-1) .lt. rmin ) .and. kk .lt. npts )
        kk = kk + 1
      end do
      !Interpolate |dcos\delta/dr| function
      cosfac = (dcosdr(kk)-dcosdr(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
      cosfac = cosfac + dcosdr(kk)
      if( adensity .eq. 1 )then
        xinorm = dglpfac(re,spin,h)**4 * cosfac / dareafac(re,spin)
        xinorm = xinorm * re**(-1.5) * ( 1.0 - sqrt(rin/re) )**2
      else
        xinorm = dglpfac(re,spin,h)**4 * cosfac / dareafac(re,spin)
      end if
      
! Construct the transfer function by summing over all pixels ***Could paralellize this loop***
      logxihi = 0.0
      odisc = 1
      i = nro + 1
      do while( odisc .eq. 1 .and. i .gt. 1 )
        i = i - 1
        odisc = 0
        do j = 1,NPHI
          phin  = (j-0.5) * 2.d0 * pi / dble(nphi) 
          alpha = rn(i) * sin(phin)
          beta  = -rn(i) * cos(phin) * mueff
          !If the ray hits the disk, calculate flux and time lag
          if( pem1(j,i) .gt. 0.0d0 )then
            re    = re1(j,i)
            if( re .gt. rin .and. re .lt. rout )then
              odisc = 1
              taudo = taudo1(j,i)
              g     = dlgfac( spin,mu0,alpha,re )
              !Calculate tausd by interpolating tlp
              kk = 2
              do while( ( rlp(kk) .le. re .or. rlp(kk-1) .lt. rmin ) .and. kk .lt. npts )
                kk = kk + 1
              end do
              !Interpolate (or extrapolate) the time function
              tausd = (tlp(kk)-tlp(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
              tausd = tausd + tlp(kk)
              tau   = tausd+taudo-tauso !This is the time lag between direct and reflected photons 
              !Interpolate |dcos\delta/dr| function
              cosfac = (dcosdr(kk)-dcosdr(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
              cosfac = cosfac + dcosdr(kk)              
              if( kk .eq. npts )then !i.e. extrapolate to Newtonian if needs be
                cosfac = dcosdr(kk) *  ( (h-honr*rlp(kk))**2 + rlp(kk)**2 )**1.5 / rlp(kk)
                cosfac = cosfac * re / ( (h-honr*re     )**2 + re**2      )**1.5
              end if             
              !Add flux from pixel to transfer function
              gsd        = dglpfac(re,spin,h)
              emissivity = gsd**(2.d0+Gamma)
              emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
              dFe        = emissivity * g**3 * domega(i) / (1.d0+zcos)**3
              frobs      = frobs + 2.0 * g**3 * gsd**3 * cosfac/dareafac(re,spin) * domega(i)
              gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
              gbin = MAX( 1    , gbin  )
              gbin = MIN( gbin , ne    )              
              !Calculate emission angle and work out which mue bin to add to
              mue   = demang(spin,mu0,re,alpha,beta)
              mubin = ceiling( mue * dble(me) )
              !Calculate incidence angle
              mus = (cosd(kk)-cosd(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
              mus = mus + cosd(kk)
              if( kk .eq. npts )then !i.e. extrapolate to Newtonian if needs be
                mus = h/sqrt(h**2+re**2) - h/sqrt(h**2+rout**2) + cosdout
              end if              
              mui = dinang(spin,re,h,mus)
              !write(285,*)re,mue,mui
              !Calculate ionisation parameter and effective ionisation parameter
              xir = gsd**4 * cosfac / dareafac(re,spin)
              xir = xir / xinorm
              if( adensity .eq. 1 ) xir = xir * re**(-1.5) * ( 1.0 - sqrt(rin/re) )**2
              logxir   = logxip + log10(xir)
              logxieff = logxir - 0.1505 - log10(mui)
              logxihi  = max( logxihi , logxieff )
              !Calculate which logxieff bin to add to
              xbin = min( xe , ceiling( (logxieff-ximin)/(ximax-ximin) * float(xe) ) )
              xbin = max( 1 , xbin )
              !Calculate which gsd bin to add to
              sdbin = min( ge , ceiling( (gsd-sdmin)/(sdmax-sdmin) * float(ge) ) )
              sdbin = max( 1 , sdbin )
              !write(124,*)re,mue
              !Sum up over frequency range (if flo=fhi=0, this is DC component)
              do fbin = 1,nf
                cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                transe(gbin,mubin,sdbin,xbin) = transe(gbin,mubin,sdbin,xbin) + real(dFe) * cexp
              end do
            end if
          end if
        end do
      end do
      !write(124,*)"no no"
      !write(285,*)"no no"
      
      xbinhi = min( xe , ceiling( (logxihi-ximin)/(ximax-ximin) * float(xe) ) )
      if( xe .eq. 1 ) xbinhi = 1
      
! Now trace rays for that bigger camera (obviously a lot easier)
      cos0 = mu0
      sin0 = sqrt( 1.0 - mu0**2 )
      do i = 1,nron
        do j = 1,nphin
          phin  = (j-0.5) * 2.d0 * pi / dble(nphin) 
          alpha = rnn(i) * sin(phin)
          beta  = -rnn(i) * cos(phin) * mueff
          call drandphithick(alpha,beta,mu0,mudisk,re,phie)
          !If the ray hits the disk, calculate flux and time lag
          if( re .gt. rin .and. re .lt. rout )then
            g = dlgfac( spin,mu0,alpha,re )
            !Time lag
            tau = sqrt(re**2+(h-honr*re)**2) - re*(sin0*sindisk*cos(phie)+mu0*mudisk ) + h*mu0
            !Interpolate (or most likely extrapolate) |dcos\delta/dr| function
            kk = 2
            do while( ( rlp(kk) .le. re .or. rlp(kk-1) .lt. rmin ) .and. kk .lt. npts )
              kk = kk + 1
            end do
            cosfac = (dcosdr(kk)-dcosdr(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
            cosfac = cosfac + dcosdr(kk)
            if( kk .eq. npts )then !i.e. extrapolate to Newtonian if needs be
              cosfac = dcosdr(kk) *  ( (h-honr*rlp(kk))**2 + rlp(kk)**2 )**1.5 / rlp(kk)
              cosfac = cosfac * re / ( (h-honr*re     )**2 + re**2      )**1.5
            end if
            !Add flux from pixel to transfer function
            gsd        = dglpfac(re,spin,h)
            emissivity = gsd**(2.d0+Gamma)
            emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
            dFe        = emissivity * g**3 * domegan(i) / (1.d0+zcos)**3
            frobs      = frobs + 2.0 * g**3 * gsd**3 * cosfac/dareafac(re,spin) * domega(i)
            gbin = ceiling( log10( g/(1.d0+zcos) ) / dloge ) + ne / 2
            gbin = MAX( 1    , gbin  )
            gbin = MIN( gbin , ne    )
            !Calculate emission angle and work out which mue bin to add to
            mue   = demang(spin,mu0,re,alpha,beta)
            mubin = ceiling( mue * dble(me) )
            !Calculate incidence angle
            mus = (cosd(kk)-cosd(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
            mus = mus + cosd(kk)
            if( kk .eq. npts )then !i.e. extrapolate to Newtonian if needs be
              mus = h/sqrt(h**2+re**2) - h/sqrt(h**2+rout**2) + cosdout
            end if              
            mui = dinang(spin,re,h,mus)
            !write(285,*)re,mue,mui
            !Calculate ionisation parameter and effective ionisation parameter
            xir = gsd**4 * cosfac / dareafac(re,spin)
            xir = xir / xinorm
            if( adensity .eq. 1 ) xir = xir * re**(-1.5) * ( 1.0 - sqrt(rin/re) )**2
            logxir   = logxip + log10(xir)
            logxieff = logxir - 0.1505 - log10(mui)              
            !Calculate which logxieff bin to add to
            xbin = min( xe , ceiling( (logxieff-ximin)/(ximax-ximin) * float(xe) ) )
            xbin = max( 1 , xbin )
            !Calculate which gsd bin to add to
            sdbin = min( ge , ceiling( (gsd-sdmin)/(sdmax-sdmin) * float(ge) ) )
            sdbin = max( 1 , sdbin )            
            !write(124,*)re,mue
            !Sum up over frequency range (if flo=fhi=0, this is DC component)
            do fbin = 1,nf
              cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
              transe(gbin,mubin,sdbin,xbin) = transe(gbin,mubin,sdbin,xbin) + real(dFe) * cexp
            end do
          end if
        end do
      end do

      !Deal with edge effects
      do xbin = 1,xe
        do sdbin = 1,ge
          do mubin = 1,me
            transe(1,mubin,sdbin,xbin)  = 0.0
            transe(1,mubin,sdbin,xbin) = 0.0
          end do
        end do
      end do
      !Normalise
      transe = transe / real(nf)

      !Finish calculation of the reflection fraction
      frobs = frobs / (dgsofac(spin,h))**3 / lens
        
      return
    end subroutine strans
!-----------------------------------------------------------------------
