!-----------------------------------------------------------------------
subroutine rtrans(spin,h,mu0,Gamma,rin,rout,honr,d,rnmax,zcos,nro,nphi,ne,dloge,&
           nf,fhi,flo,me,xe,rlxi, lognep, transe,transe_a,frobs,frrel,lens&
           &,logxir,gsdr, logner)
! Code to calculate the transfer function for an accretion disk.
! This code first does full GR ray tracing for a camera with impact parameters < bmax
! It then also does straight line ray tracing for impact parameters >bmax
! It adds both up to produce a transfer function for a disk extending from rin to rout
! INPUT
! spin,h,mu0,Gamma      Physical parameters (spin, source height, cos(inclination), photon index)
! rin,rout,honr         Physical parameters (disk inner radius, outer radius & scaleheight)
! d,rnmax               Physical parameters (distance of the source, max radius for which GR ray tracing is used)
! zcos                  Cosmological redshift
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
! lens                  Lensing factor for direct emission
! logxir(xe),gsdr(xe)   logxi (ionization parameter) and gsd (source to disc blueshift) as a function of radius
  use dyn_gr
  use blcoordinate
  implicit none
  integer nro,nphi,ndelta,ne,nf,me,xe
  double precision spin,h,mu0,Gamma,rin,rout,zcos,fhi,flo,honr&
       &,cosdout,logxir(xe),gsdr(xe), logner(xe)
  real rlxi, dloge, lognep
  complex cexp,transe(ne,nf,me,xe),transe_a(ne,nf,me,xe)
  integer i,npts,j,odisc,n,gbin,rbin,mubin
  parameter (ndelta=1000)
  double precision domega(nro),d,taudo,g,dlgfac,dFe,lens,newtex
  double precision tauso,rlp(ndelta),dcosdr(ndelta),tlp(ndelta),cosd(ndelta)
  double precision alpha,beta,cos0,sin0,phie,re,gsd
  double precision tau,tausd,emissivity,cosfac,dglpfac,dareafac
  integer kk,fbin,get_index
  double precision rmin,disco,rfunc,scal,velocity(3),mudisk,sysfref
  double precision rnmax,rnmin,rn(nro),phin,mueff,dlogr,interper
  double precision fi(nf),dgsofac,sindisk,mue,demang,frobs,cosdin,frrel
  integer nron,nphin,nrosav,nphisav,verbose
  double precision spinsav,musav,routsav,mudsav,rnn(nro),domegan(nro)
  logical dotrace
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
      
! Set up observer's camera ( alpha = rn sin(phin), beta = mueff rn cos(phin) )
! to do full GR ray tracing with      
  mueff  = max( mu0 , 0.3d0 )
  rnmin  = rfunc(spin,mu0)
  !Grid to do in full GR
  call getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
  !Grid for Newtonian approximation
  call getrgrid(rnmax,rout,mueff,nron,nphin,rnn,domegan)
      
! Set frequency array
  do fbin = 1,nf
    fi(fbin) = flo * (fhi/flo)**((float(fbin)-0.5d0)/dble(nf))
  end do
  if( fhi .lt. tiny(fhi) ) fi(1) = 0.0d0

! Calculate lensing factor "lens" and source to observer time lag "tauso"
  call getlens(spin,h,mu0,lens,tauso)
  if( tauso .ne. tauso ) stop "tauso is NaN"
      
! Calculate dcos/dr and time lags vs r for the lamppost model
  call getdcos(spin,h,mudisk,ndelta,rout,npts,rlp,dcosdr,tlp,cosd,cosdout)

  ! Set up grids for ionization and g-factor as a function of radius
  call radfunctions(xe,rin,rnmax,dble(rlxi),dble(lognep), spin,h,honr,rlp,dcosdr&
       &,cosd,ndelta,rmin,npts,logxir,gsdr, logner)
  !Outputs: logxir(1:xe),gsdr(1:xe), logner(1:xe)
 
! Calculate the relxill reflection fraction
  frrel = sysfref(rin,rlp,cosd,ndelta,cosdout)      
      
! Trace rays in full GR for the small camera
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

! Construct the transfer function by summing over all pixels
  dlogr    = log10(rnmax/rin) / real(xe-1)
  cos0     = mu0
  sin0     = sqrt(1.0-cos0**2)
  transe   = 0.0 !Initialised transfer function
  transe_a = 0.0 !Initialised transfer function
  frobs    = 0.0 !Initialised observer's reflection fraction
  odisc    = 1
  i        = nro + 1
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
              !Find the rlp bin that corresponds to re
              kk = get_index(rlp,ndelta,re,rmin,npts)
              !Interpolate (or extrapolate) the time function
              tausd = interper(rlp,tlp,ndelta,re,kk)
              tau   = (1.d0+zcos) * (tausd+taudo-tauso) !This is the time lag between direct and reflected photons
              !Interpolate |dcos\delta/dr| function
              cosfac = interper(rlp,dcosdr,ndelta,re,kk)
              !Extrapolate to Newtonian if needs be
              if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
              !Calculate flux from pixel
              gsd        = dglpfac(re,spin,h)
              emissivity = gsd**Gamma
              emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
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
              !Calculate emission angle and work out which mue bin to add to
              mue   = demang(spin,mu0,re,alpha,beta)
              mubin = ceiling( mue * dble(me) )
              !Add to the transfer function integral
              do fbin = 1,nf
                 cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                 transe(gbin,fbin,mubin,rbin)   = transe(gbin,fbin,mubin,rbin)                  + real(dFe) * cexp
                 transe_a(gbin,fbin,mubin,rbin) = transe_a(gbin,fbin,mubin,rbin) + real(log(g)) * real(dFe) * cexp
              end do
           end if
        end if
     end do
  end do

  ! verbose = myenv("REV_VERB",0)     !Set verbose level
  ! if( verbose .gt. 0 )then
  !    write(*,*)"[gsd^{2-Gamma} epsilon(r)]max=",rfacmax
  ! end if

! Now trace rays for that bigger camera (obviously a lot easier)
  do i = 1,nron
     do j = 1,nphin
        phin  = (j-0.5) * 2.d0 * pi / dble(nphin) 
        alpha = rnn(i) * sin(phin)
        beta  = -rnn(i) * cos(phin) * mueff
        call drandphithick(alpha,beta,mu0,mudisk,re,phie)
        !If the ray hits the disk, calculate flux and time lag
        if( re .gt. rin .and. re .lt. rout )then
           g = dlgfac( spin,mu0,alpha,re )
           !Find the rlp bin that corresponds to re
           kk = get_index(rlp,ndelta,re,rmin,npts)
           !Time lag
           tau = sqrt(re**2+(h-honr*re)**2) - re*(sin0*sindisk*cos(phie)+mu0*mudisk ) + h*mu0
           tau = (1.d0+zcos) * tau
           !Interpolate |dcos\delta/dr| function
           cosfac = interper(rlp,dcosdr,ndelta,re,kk)
           !Extrapolate to Newtonian if needs be (most likely in this loop)
           if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
           !Calculate flux from pixel
           gsd        = dglpfac(re,spin,h)
           emissivity = gsd**Gamma
           emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
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
           !Calculate emission angle and work out which mue bin to add to
           mue   = demang(spin,mu0,re,alpha,beta)
           mubin = ceiling( mue * dble(me) )
           !Add to the transfer function integral
           do fbin = 1,nf
              cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
              transe(gbin,fbin,mubin,rbin)   = transe(gbin,fbin,mubin,rbin)                  + real(dFe) * cexp
              transe_a(gbin,fbin,mubin,rbin) = transe_a(gbin,fbin,mubin,rbin) + real(log(g)) * real(dFe) * cexp
           end do
        end if
     end do
  end do

  ! !Deal with edge effects
  ! do rbin = 1,xe
  !    do mubin = 1,me
  !       do fbin = 1,nf
  !          transe(1,fbin,mubin,rbin)    = 0.0
  !          transe(ne,fbin,mubin,rbin)   = 0.0
  !          transe_a(1,fbin,mubin,rbin)  = 0.0
  !          transe_a(ne,fbin,mubin,rbin) = 0.0
  !       end do
  !    end do
  ! end do
      
  !Finish calculation of observer's reflection fraction
  frobs = frobs / dgsofac(spin,h) / lens
    
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
subroutine radfunctions(xe,rin,rnmax,logxip, lognep, spin,h,honr,rlp,dcosdr&
     &,cosd,ndelta,rmin,npts,logxir,gsdr, logner)
! In  : xe,rin,rnmax,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
! Out : logxir(1:xe),gsdr(1:xe)
  implicit none
  integer xe,ndelta,npts,adensity
  double precision rin,rnmax,logxip,lognep,spin,h,honr,rlp(ndelta),dcosdr(ndelta)
  double precision cosd(ndelta),rmin,logxir(xe),gsdr(xe), logner(xe)
  integer i,kk,get_index,myenv
  double precision rp,logxinorm,logxiraw,mylogne,re,mus,interper,newtex
  double precision mui,dinang,gsd,dglpfac
  !Decide on zone a density profile or constant density profile
  adensity = myenv("A_DENSITY",1)
  adensity = min( adensity , 1 )
  adensity = max( adensity , 0 )

  !Calculate normalisation of logxi(r) function
  ! The last numeric factor is to avoid to have rp = rin in case of A_DENSITY=0. This would cause a problem in the function mylogne (it'd go to infinity)
  rp = rin * ( 11.0 / 9.0 )**(2*adensity)  + 0.0000001d0
  
  logxinorm = logxiraw(rp,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,gsd)
  logxinorm = logxinorm - adensity * mylogne(rp,rin)
  
  !Now calculate logxi itself
  rp = rin * ( 11.0 / 9.0 )**2
  do i = 1, xe - 1 
     !Calculate the radius for this bin
     if( i .eq. 1 )then
        re = 10.0**( 0.5 * ( log10(rin) + log10(rp) ) )
     else
        re =     (rnmax/rp)**(real(i-1)  /real(xe-1))
        re = re + (rnmax/rp)**(real(i-2)/real(xe-1))
        re = re * rp / 2.0        
     end if
!Calculate the incident angle for this bin
     kk = get_index(rlp,ndelta,re,rmin,npts)
     mus = interper(rlp,cosd,ndelta,re,kk)
     if( kk .eq. npts ) mus = newtex(rlp,cosd,ndelta,re,h,honr,kk)
     mui = dinang(spin,re,h,mus)
!Now logxi(r)
     logxir(i) = logxiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,gsd)
     logxir(i) = logxir(i) - adensity * mylogne(re,rin)
     logxir(i) = logxir(i) - logxinorm + logxip
     logxir(i) = logxir(i) - 0.1505 - log10(mui)
     logxir(i) = max( logxir(i) , 0.d0  )
     logxir(i) = min( logxir(i) , 4.7d0 )

!Also save gsd(r)
     gsdr(i) = gsd

     logner(i) = lognep + adensity * mylogne(re, rin)
! Check if the density is in the limits 
     logner(i) = max( logner(i) , 15.d0  )
     logner(i) = min( logner(i) , 19.d0 )
     write(*,*) 'logner ', logner(i)
  end do
  logxir(xe) = 0.0
  gsdr(xe)   = dglpfac(1000.d0,spin,h)
!The last bin has the maximum density available a part in the case A_DENSITY=0 in that case it should have the density assigned by the user (parameter in the model)
  logner(xe) = lognep + adensity * 4.d0
  return
end subroutine radfunctions
!-----------------------------------------------------------------------


  
!-----------------------------------------------------------------------
function logxiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,gsd)
! In: re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts
! Out: logxiraw,gsd
  implicit none
  integer ndelta,npts,kk,get_index
  double precision re,spin,h,honr,rlp(ndelta),dcosdr(ndelta),rmin,gsd
  double precision cosfac,interper,logxiraw,dareafac,dglpfac,newtex
  !Calculate source to disc blueshift at this radius
  gsd = dglpfac(re,spin,h)
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







!***************************************************



!-----------------------------------------------------------------------
      subroutine strans(spin,h,mu0,Gamma,rin,rout,rnmax,d,honr,zcos,nro,nphi,ne,dloge,&
           nf,fhi,flo,ximin,ximax,me,ge,xe,rlxi,sdmin,sdmax,transe,transe_a,frobs,frrel,xbinhi,lens)
! Code to calculate the transfer function for an accretion disk.
! This code first does full GR ray tracing for a camera with impact parameters < bmax
! It then also does straight line ray tracing for impact parameters >bmax
! It adds both up to produce a transfer function for a disk extending from rin to rout
! INPUT
! spin,h,mu0,Gamma      Physical parameters (spin, source height, cos(inclination), photon index)
! rin,rout,honr         Physical parameters (disk inner radius, outer radius & scaleheight)
! d,rnmax               Physical parameters (distance of the source, )
! zcos                  Cosmological redshift
! nro,nphi              Number of pixels on the observer's camera (b and phib)
! ndelta                Number of \delta values considered
! ne, dloge             Number of energy bins and maximum energy (compatible with FFT convolution)
! ear(0:ne)             Energy grid
! nf,fhi,flo            nf = Number of logarithmic frequency bins averaged over, range= flo to fhi
! ximin ximax           Minimum and maximum values of logxi_eff (depends on emissivity, density and incidence angle)
! me
! ge
! xe
! rlxi
! OUTPUT
! sdmin,sdmax
! transe(nf,ne,mex,gex,xex)    Transfer function as a function of energy and emission angle for the given frequency range
! transe_a(nf,ne,mex,gex,xex)  Second transfer function as a function of energy and emission angle for the given frequency range
!                            This is for the non-linear effects
! frobs                 Observer's reflection fraction
! frrel                 Reflection fraction defined by relxilllp
! xbinhi                Highest xi bin that has an entry 
! lens                  Lensing factor for direct emission
        use dyn_gr
        use blcoordinate
      implicit none
      integer nro,nphi,ndelta,ne,nf,sdbin,me,ge,xe
      double precision spin,h,mu0,Gamma,rin,rout,zcos,fhi,flo,honr,cosdout
      real rlxi,ximin,ximax
      real dloge,sdmin,sdmax
      complex cexp,transe(ne,nf,me,ge,xe),transe_a(ne,nf,me,ge,xe)
      integer i,npts,j,odisc,n,xbin
      parameter (ndelta=1000)
      double precision domega(nro),d
      double precision tauso,rlp(ndelta),dcosdr(ndelta),tlp(ndelta),cosd(ndelta)
      double precision alpha,beta,cos0,sin0,phie,re,gsd
      double precision taudo,g,dlgfac,dFe,lens
      double precision tau,tausd,emissivity,cosfac,dglpfac,dareafac
      integer gbin,kk,fbin,ifl,myenv
      double precision rmin,disco,rfunc
      double precision scal,velocity(3)
      double precision mudisk
      double precision rnmax,rnmin,rn(nro),phin,mueff
      double precision fi(nf),dgsofac,sindisk,mue,demang,frobs,cosdin,frrel
      integer nron,nphin,nrosav,nphisav,mubin,xbinhi,adensity,verbose
      double precision spinsav,musav,routsav,mudsav,rnn(nro),domegan(nro)
      double precision mus,mui,dinang,xir,logxieff,logxir,xinorm,logxip,logxihi,rfacmax
      logical dotrace
!      real :: t0,t1,t2,t3,t4
      !      character (len=1) A_DENSITY
      !      double precision amp,dtaudo,f,ebin,f1234(4),lambda,line,mucros,pem,phase,q,phi0,sigmacros,sum
      !      integer jj,k,l,stat
      data nrosav,nphisav,spinsav,musav /0,0,2.d0,2.d0/
      save nrosav,nphisav,spinsav,musav,routsav,mudsav
      
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
      !Grid to do in full GR
      call getrgrid(rnmin,rnmax,mueff,nro,nphi,rn,domega)
      !Grid for Newtonian approximation
      call getrgrid(rnmax,rout,mueff,nron,nphin,rnn,domegan)
      
! Set frequency array
      do fbin = 1,nf
        fi(fbin) = flo * (fhi/flo)**((float(fbin)-0.5d0)/dble(nf))
      end do
      if( fhi .lt. 1d-10 ) fi(1) = 0.0d0
      

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

      ! write(*,*)"---------------------"
      ! write(*,*)"Reflection fraction calculations"
      ! write(*,*)"cosdin=",cosdin
      ! write(*,*)"cosdout=",cosdout
      ! write(*,*)"---------------------"    
      
      
! Now trace rays in full GR down to mu=0 for impact parameters with b<bmax
      cos0  = mu0
      sin0  = sqrt(1.0-cos0**2)
      transe  = 0.0 !Initialised transfer function
      transe_a  = 0.0 !Initialised transfer function
      frobs   = 0.0

      
!Check if we need to calculate re and tau and pem or we took them from a grid      
      if(status_re_tau) then
         
! Trace rays in full GR for the small camera
! to convert alpha and beta to r and tau_do (don't care about phi)
! First work out if we even need to call this
         dotrace = .false.
         ! if( nro .ne. nrosav ) dotrace = .true.
         ! if( nphi .ne. nphisav ) dotrace = .true.
         if( abs(spinsav-spin) .gt. 1d-6 ) dotrace = .true.
         if( abs(musav-mu0) .gt. 1d-6 ) dotrace = .true.
         if( abs(routsav-rout) .gt. 1d-6 ) dotrace = .true.
         if( abs(mudsav-mudisk) .gt. 1d-6 ) dotrace = .true.         
         if( dotrace )then
            call GRtrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
            spinsav = spin
            musav   = mu0
            routsav = rout
            mudsav  = mudisk
         end if
      end if
      

! Decide on zone a density profile or constant density profile
      adensity = myenv("A_DENSITY",1)

! Set up logxi part of the calculation
      logxip = dble(rlxi)
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
        xinorm = dglpfac(re,spin,h)**2 * cosfac / dareafac(re,spin)
        xinorm = xinorm * re**(-1.5) * ( 1.0 - sqrt(rin/re) )**2
      else
        xinorm = dglpfac(re,spin,h)**2 * cosfac / dareafac(re,spin)
      end if


!      call CPU_TIME(t0)
!      write(*,*) 'NPHI', NPHI
!      write(*,*) 'nf',nf
!      t4 = 0.0
! Construct the transfer function by summing over all pixels ***Could paralellize this loop***

      rfacmax = -huge(rfacmax)
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
              tau   = (1.d0+zcos) * (tausd+taudo-tauso) !This is the time lag between direct and reflected photons
              !Interpolate |dcos\delta/dr| function
              cosfac = (dcosdr(kk)-dcosdr(kk-1))*(re-rlp(kk))/(rlp(kk)-rlp(kk-1))
              cosfac = cosfac + dcosdr(kk)              
              if( kk .eq. npts )then !i.e. extrapolate to Newtonian if needs be
                cosfac = dcosdr(kk) *  ( (h-honr*rlp(kk))**2 + rlp(kk)**2 )**1.5 / rlp(kk)
                cosfac = cosfac * re / ( (h-honr*re     )**2 + re**2      )**1.5
              end if             
              !Add flux from pixel to transfer function
              gsd        = dglpfac(re,spin,h)
!              emissivity = gsd**(2.d0+Gamma)
              emissivity = gsd**(Gamma)
              emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
              dFe        = emissivity * g**3 * domega(i) / (1.d0+zcos)**3
              frobs      = frobs + 2.0 * g**3 * gsd * cosfac/dareafac(re,spin) * domega(i)
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
              xir = gsd**2 * cosfac / dareafac(re,spin)
              rfacmax = max( rfacmax , 0.5*xir )
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
!              call cpu_time(t2)
              do fbin = 1,nf
                 cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                 transe(gbin,fbin,mubin,sdbin,xbin) = real(dFe) * cexp
                 transe_a(gbin,fbin,mubin,sdbin,xbin) = real(log(g)) * real(dFe) * cexp
              end do
              ! call cpu_time(t3)
              ! t4 = t4 + (t3-t2)
              
            end if
          end if
        end do
      end do

      verbose = myenv("REV_VERB",0)     !Set verbose level
      if( verbose .gt. 0 )then
         write(*,*)"[gsd^{2-Gamma} epsilon(r)]max=",rfacmax
         write(*,*)"A_DENSITY=",adensity
      end if
      
      xbinhi = min( xe , ceiling( (logxihi-ximin)/(ximax-ximin) * float(xe) ) )
      if( xe .eq. 1 ) xbinhi = 1
      
!      call CPU_TIME(t0)     
      ! write(*,*) 'nron', nron
      ! write(*,*) 'nphin',nphin
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
            tau = (1.d0+zcos) * tau
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
!            emissivity = gsd**(2.d0+Gamma)
            emissivity = gsd**(Gamma)
            emissivity = 0.5 * emissivity * cosfac / dareafac(re,spin)
            dFe        = emissivity * g**3 * domegan(i) / (1.d0+zcos)**3
            frobs      = frobs + 2.0 * g**3 * gsd * cosfac/dareafac(re,spin) * domega(i)
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
            xir = gsd**2 * cosfac / dareafac(re,spin)
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
               transe(gbin,fbin,mubin,sdbin,xbin)   = real(dFe) * cexp
               transe_a(gbin,fbin,mubin,sdbin,xbin) = real(log(g)) * real(dFe) * cexp
            end do
          end if
        end do
     end do


!           call CPU_TIME(t1)
!      write(*,*)"Newton kernel CPU time=",t1-t0

      ! !Deal with edge effects
      ! do xbin = 1,xe
      !   do sdbin = 1,ge
      !     do mubin = 1,me
      !       transe(1,mubin,sdbin,xbin)  = 0.0
      !       transe(ne,mubin,sdbin,xbin) = 0.0
      !       transe_a(1,mubin,sdbin,xbin)  = 0.0
      !       transe_a(ne,mubin,sdbin,xbin) = 0.0
      !     end do
      !   end do
      ! end do
      ! !Normalise
      ! transe = transe / real(nf)
      ! transe_a = transe_a / real(nf)
      
      !Finish calculation of the reflection fraction
      frobs = frobs / dgsofac(spin,h) / lens
    
      return
    end subroutine strans
!-----------------------------------------------------------------------
