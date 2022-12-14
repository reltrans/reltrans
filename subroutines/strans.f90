!-----------------------------------------------------------------------
subroutine rtrans(spin,h,mu0,Gamma,rin,rout,honr,d,rnmax,zcos,b1,b2,qboost,&
     &dset,fcons,nro,nphi,ne,dloge,nf,fhi,flo,me,xe,rlxi,lognep,&
     &transe,transe_a,frobs,frrel,lens,logxir,gsdr,logner)
! Code to calculate the transfer function for an accretion disk.
! This code first does full GR ray tracing for a camera with impact parameters < bmax
! It then also does straight line ray tracing for impact parameters >bmax
! It adds both up to produce a transfer function for a disk extending from rin to rout
! INPUT
! spin,h,mu0,Gamma      Physical parameters (spin, source height, cos(inclination), photon index)
! rin,rout,honr         Physical parameters (disk inner radius, outer radius & scaleheight)
! d,rnmax               Physical parameters (distance of the source, max radius for which GR ray tracing is used)
! zcos                  Cosmological redshift
! b1                    Linear coefficient of angular emissivity function
! b2                    Quadratic coefficient of angular emissivity function
! qboost                Asymmetry parameter of angular emissivity function
! dset                  dset=1 means calculate ionization from distance, dset=0 means ignore distance
! fcons                 Used to calculate ionization from distance
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
  integer nro,nphi,ndelta,ne,nf,me,xe,dset
  double precision spin,h,mu0,Gamma,rin,rout,zcos,fhi,flo,honr
  double precision b1,b2,qboost
  double precision fcons,cosdout,logxir(xe),gsdr(xe),logner(xe)
  real rlxi, dloge, lognep
  complex cexp,transe(ne,nf,me,xe),transe_a(ne,nf,me,xe)
  integer i,npts,j,odisc,n,gbin,rbin,mubin
  parameter (ndelta=1000)
  double precision domega(nro),d,taudo,g,dlgfacthick,dFe,lens,newtex
  double precision tauso,rlp(ndelta),dcosdr(ndelta),tlp(ndelta),cosd(ndelta)
  double precision alpha,beta,cos0,sin0,phie,re,gsd
  double precision tau,tausd,emissivity,cosfac,dglpfacthick,dareafac
  integer kk,fbin,get_index
  double precision rmin,disco,rfunc,scal,velocity(3),mudisk,sysfref
  double precision rnmax,rnmin,rn(nro),phin,mueff,dlogr,interper
  double precision fi(nf),dgsofac,sindisk,mue,demang,frobs,cosdin,frrel
  double precision pnorm,mus,ptf,pfunc_raw,cosdelta_obs,ang_fac
  integer nron,nphin,nrosav,nphisav,verbose
  double precision spinsav,musav,routsav,mudsav,rnn(nro),domegan(nro)
  integer myenv
  double precision lximax
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
  call getlens(spin,h,mu0,lens,tauso,cosdelta_obs)
  if( tauso .ne. tauso ) stop "tauso is NaN"
      
! Calculate dcos/dr and time lags vs r for the lamppost model
  call getdcos(spin,h,mudisk,ndelta,rout,npts,rlp,dcosdr,tlp,cosd,cosdout)

  ! Set up grids for ionization and g-factor as a function of radius
  if( dset .eq. 0 )then
     call radfunctions_dens(xe,rin,rnmax,dble(rlxi),dble(lognep), spin,h,honr,rlp,dcosdr&
          &,cosd,ndelta,rmin,npts,logxir,gsdr, logner)
     pnorm = 1.d0 / ( 4.d0 * pi )
  else
     call radfuncs_dist(xe, rin, rnmax, b1, b2, qboost, fcons,&
     & dble(lognep), spin, h, honr, rlp, dcosdr, cosd, ndelta, rmin, npts,&
     & logxir, gsdr, logner, pnorm)
  end if
  !Outputs: logxir(1:xe),gsdr(1:xe), logner(1:xe)
  
! Calculate 4pi p(theta0,phi0) = ang_fac
  ang_fac = 4.d0 * pi * pnorm * pfunc_raw(-cosdelta_obs,b1,b2,qboost)
! Adjust the lensing factor (easiest way to keep track)
  lens = lens * ang_fac
     
! Calculate the relxill reflection fraction
  frrel = sysfref(rin,rlp,cosd,ndelta,cosdout)      
      
! Trace rays in full GR for the small camera (ie with relativistic effects) from the osberver to the disk, which is why it doesnt depend on h
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
  do while( odisc .eq. 1 .and. i .gt. 1 )         !main loops of the subroutine: first is for GR
     i = i - 1                                      !i counts over the camera until it reaches the disk inner radius
     odisc = 0
     do j = 1,NPHI                              !azimuth over BH on the disk
        phin  = (j-0.5) * 2.d0 * pi / dble(nphi) 
        alpha = rn(i) * sin(phin)
        beta  = -rn(i) * cos(phin) * mueff
        !If the ray hits the disk, calculate flux and time lag
        if( pem1(j,i) .gt. 0.0d0 )then
           re    = re1(j,i)
           if( re .gt. rin .and. re .lt. rout )then              
              odisc = 1                             
              taudo = taudo1(j,i)           !these calculate g factors, emissivity, etc....and will change with a new lamppost
              g     = dlgfacthick( spin,mu0,alpha,re,mudisk ) !ideally tbd so that we eventually can do this for N lamposts, not 2
              !g = dlgfac( spin,mu0,alpha,re )                  !which basically means generalize interper and dlgfacthick
              !Find the rlp bin that corresponds to re
              kk = get_index(rlp,ndelta,re,rmin,npts)
              !Interpolate (or extrapolate) the time function
              tausd = interper(rlp,tlp,ndelta,re,kk)
              tau   = (1.d0+zcos) * (tausd+taudo-tauso) !This is the time lag between direct and reflected photons
              !Interpolate |dcos\delta/dr| function
              cosfac = interper(rlp,dcosdr,ndelta,re,kk)
              mus    = interper(rlp,cosd,ndelta,re,kk)
              !Extrapolate to Newtonian if needs be
              if( kk .eq. npts )then
                 cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
                 mus    = newtex(rlp,cosd,ndelta,re,h,honr,kk)
              end if
              !Calculate angular emissivity
              ptf        = pnorm * pfunc_raw(-mus,b1,b2,qboost)
              !Calculate flux from pixel
              gsd        = dglpfacthick(re,spin,h,mudisk)
              !gsd        = dglpfac(re,spin,h)
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
              !Calculate emission angle and work out which mue bin to add to
              mue   = demang(spin,mu0,re,alpha,beta)
              mubin = ceiling( mue * dble(me) )
              !Add to the transfer function integral
              !gbin: energy; fbin: Fourier frequency bin; mubin: emission angle bin; rbin: radial bin
              do fbin = 1,nf
                 cexp = cmplx( cos(real(2.d0*pi*tau*fi(fbin))) , sin(real(2.d0*pi*tau*fi(fbin))) )
                 transe(gbin,fbin,mubin,rbin)   = transe(gbin,fbin,mubin,rbin)                  + real(dFe) * cexp
                 transe_a(gbin,fbin,mubin,rbin) = transe_a(gbin,fbin,mubin,rbin) + real(log(g)) * real(dFe) * cexp
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
           g = dlgfacthick( spin,mu0,alpha,re,mudisk )
           !g = dlgfac( spin,mu0,alpha,re )
           !Find the rlp bin that corresponds to re
           kk = get_index(rlp,ndelta,re,rmin,npts)
           !Time lag
           tau = sqrt(re**2+(h-honr*re)**2) - re*(sin0*sindisk*cos(phie)+mu0*mudisk ) + h*mu0
           tau = (1.d0+zcos) * tau
           !Interpolate |dcos\delta/dr| function
           cosfac = interper(rlp,dcosdr,ndelta,re,kk)
           mus    = interper(rlp,cosd,ndelta,re,kk)
           !Extrapolate to Newtonian if needs be
           if( kk .eq. npts )then
              cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
              mus    = newtex(rlp,cosd,ndelta,re,h,honr,kk)
           end if
           !Calculate angular emissivity
           ptf        = pnorm * pfunc_raw(-mus,b1,b2,qboost)
           !Calculate flux from pixel
           gsd        = dglpfacthick(re,spin,h,mudisk)
           !gsd        = dglpfac(re,spin,h)
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
  use env_variables
  implicit none
  integer         , intent(IN)   :: xe, ndelta, npts 
  double precision, intent(IN)   :: rin, rmin, rnmax, b1, b2, qboost
  double precision, intent(IN)   :: fcons, lognep, spin, h, honr
  double precision, intent(IN)   :: rlp(ndelta), dcosdr(ndelta), cosd(ndelta)
  double precision, intent(INOUT):: logxieff(xe), gsdr(xe), logner(xe)
  integer          :: i, kk, get_index, myenv, verbose
  double precision :: pnorm,re,re1(xe),zA_logne,cosfac,mus,interper,newtex,mudisk
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: ptf,pfunc_raw,gsd,dglpfacthick,eps_bol,Fx(xe),logxir(xe),mui,dinang
  double precision :: pnormer,dareafac,lximax
! Set disk opening angle
  mudisk   = honr / sqrt( honr**2 + 1.d0  )
! Normalise the angular emissivity profile
  pnorm = pnormer(b1,b2,qboost)  
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
  if( verbose .gt. 10 )then
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
subroutine radfunctions_dens(xe, rin, rnmax, logxip, lognep, spin, h, honr, rlp, dcosdr&
     &, cosd, ndelta, rmin, npts, logxir, gsdr, logner)
! In  : xe,rin,rnmax,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
  ! Out : logxir(1:xe), gsdr(1:xe), logner(1:xe)
  use env_variables
  implicit none
  integer         , intent(IN)   :: xe, ndelta, npts 
  double precision, intent(IN)   :: rin, rmin, rnmax, logxip, lognep, spin, h, honr, rlp(ndelta), dcosdr(ndelta), cosd(ndelta)
  double precision, intent(INOUT):: logxir(xe), gsdr(xe), logner(xe)

  integer          :: i, kk, get_index, myenv
  double precision :: rp, logxinorm, lognenorm,  mus, interper, newtex, mui, dinang, gsd, dglpfacthick
  double precision :: logxiraw, mylogne,mudisk
  double precision, allocatable :: rad(:)
! Set disk opening angle
  mudisk   = honr / sqrt( honr**2 + 1.d0  )
  
  allocate(rad(xe))
  !Now calculate logxi itself

  
  !radius calculation 
  ! do i = 1, xe
  !    rad(i) = (rnmax/rin)**(real(i-1) / real(xe))
  !    rad(i) = rad(i) + (rnmax/rin)**(real(i) / real(xe))
  !    rad(i) = rad(i) * rin * 0.5
  !    ! write(*,*) i, rad(i)
  ! enddo
  
  ! The loop calculates the raw xi and raw n_e.
  ! This means they are without normalization: only to find the maximum and the minimum. Remember that the max of the ionisation is not the same as the minumim in the density because the flux depends on r
  !The loops calculates also the correction factor mui 
  do i = 1, xe 
     rad(i) = (rnmax/rin)**(real(i-1) / real(xe))
     rad(i) = rad(i) + (rnmax/rin)**(real(i) / real(xe))
     rad(i) = rad(i) * rin * 0.5

!Now calculate the raw density (this matters only for high dens model reltransD)
     logner(i) = adensity * mylogne(rad(i), rin)
     
!Now logxi(r)
     logxir(i) = logxiraw(rad(i),spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,mudisk,gsd)
     logxir(i) = logxir(i) - logner(i)

!Calculate the incident angle for this bin
     kk = get_index(rlp, ndelta, rad(i), rmin, npts)
     mus = interper(rlp, cosd, ndelta, rad(i), kk)
     if( kk .eq. npts ) mus = newtex(rlp, cosd, ndelta, rad(i), h, honr, kk)
     mui = dinang(spin, rad(i), h, mus)
!Correction to account for the radial dependence of incident angle
     logxir(i) = logxir(i) - 0.1505 - log10(mui)

     !Also save gsd(r)
     gsdr(i) = gsd
     ! write(*,*) 'logxir, logner', rad(i), logxir(i), logner(i)
  end do

  !After the loop calculate the max and the min
  logxinorm = maxval(logxir)
  lognenorm = minval(logner)
  ! write(*,*)'-----------------------'
  ! write(*,*) logxinorm, lognenorm  
  ! write(*,*)'-----------------------'
  logxir = logxir - (logxinorm - logxip)
  logner = logner - (lognenorm - lognep)
!check max and min for both ionisation and density
  logxir = max( logxir , 0.d0  )
  logxir = min( logxir , 4.7d0 )
  ! logner = max( logner , 15.d0  )
  ! logner = min( logner , 22.d0 )
  !...no need to enforce limits on logne since this is done in myreflect()
  !This is needed because reflionx has a different maximum to xillverDCp
  
!   do i = 1, xe      
!      logxir(i) = logxir(i) - logxinorm + logxip
!      ! Check if the density is in the limits
!      logxir(i) = max( logxir(i) , 0.d0  )
!      logxir(i) = min( logxir(i) , 4.7d0 )
! !Density profile 
!      logner(i) = logner(i) - lognenorm + lognep
     ! logner(i) = max( logner(i) , 15.d0  )
     ! logner(i) = min( logner(i) , 21.d0 )
!      write(*,*) 'logxir, logner', i, rad(i), logxir(i), logner(i)
!   enddo

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
