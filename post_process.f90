
include 'subroutines/header.h'

! Compiling and running:

! make -f postmakefile main
! ./post_process

program post_process
  !Code to extract Distance parameter out of reltransDCp fits
use env_variables
  implicit none
  logical chainmode,anynull
  real param(15)
  integer xe,verbose
  double precision dist,distance
  integer unit,status,readwrite,blocksize
  character (len=500) chainfile,comment,newchainfile
  integer steps,columns,col,k,newunit
  real chisquared,parray(15)
  integer iparam(15),i,colnum,ncols,nhdu,hdutype
  character (len=16) ttype(1),tform(1),tunit(1)

! Input parameters
                      !C     !B
  param(1)  = 11.1126    !9.6   !10.1       !h
  param(2)  = 0.998   !0.998 !0.998      !a
  param(3)  = 37.2524    !43.6  !34.1       !inc
  param(4)  = 5.43081347    !25.1  !8.35       !rin
  param(5)  = 2e4     !1e3   !1e3        !rout
  param(6)  = 0.0     !0.0   !0.0        !zcos
  param(7)  = 1.69994   !1.7   !1.714      !Gamma
  param(8)  = 3.35450   !3.07  !3.452      !logxi
  param(9)  = 2.48467     !0.88  !2.1        !Afe
  param(10) = 18.6000    !19.55 !15.6       !lognep
  param(11) = 181.115   !386.0 !122.0      !kTe
  param(12) = 0.0     !0.0   !0.0        !Nh
  param(13) = 0.292266    !0.34  !0.38       !boost
  param(14) = 25.4543    !13.9  !25.5       !Mass
  param(15) = 9.70325e-2   !0.106 !0.095      !Anorm

  chainfile = '/Users/nai47/Dropbox/reltrans/paper_dcp_sys100_mc_newxspec.out'
  newchainfile = '/Users/nai47/Dropbox/reltrans/new.out'

  !Column number corresponding to each parameter in the chain
  !(not used if chainmode=false)
  !0 means that the parameter was fixed to the value in param(:)
  iparam(1)  = 7           !h
  iparam(2)  = 0           !a
  iparam(3)  = 8           !inc
  iparam(4)  = 9           !rin
  iparam(5)  = 0           !rout
  iparam(6)  = 0           !zcos
  iparam(7)  = 10          !Gamma
  iparam(8)  = 11          !logxi
  iparam(9)  = 12          !Afe
  iparam(10) = 13          !lognep
  iparam(11) = 14          !kTe
  iparam(12) = 0           !Nh
  iparam(13) = 15          !boost
  iparam(14) = 16          !Mass
  iparam(15) = 17          !Anorm

! log(ne) Model B fit would need to have to get Dkpc = 2.2
! is 18.91853751380274
  
! Settings
  chainmode = .true.  !Reading in a chain (true) or just entering one parameter set (false)
  xe       = 20       !Number of radial zones
  adensity = 1        !1 = zone A ne; 0 = const ne
  verbose  = 0

  dist = distance(param,xe,verbose)
  write(*,*)"reltransDCp distance (kpc) = ",dist

! Read in chain
  if( chainmode )then
     
     !Open chain file
     status = 0
     call ftgiou(unit,status)
     readwrite = 0
     call ftopen(unit,chainfile,readwrite,blocksize,status)
     
     !Shift to  extension "CHAIN"
     status = 0
     call ftmnhd(unit,2,'CHAIN',0,status)
     if( status .ne. 0 ) stop 'cannot shift to extension CHAIN'

     !Read number of rows and columns
     call ftgkyj(unit,'NAXIS2',steps,comment,status)
     if(status .ne. 0) stop 'Cannot determine No of rows'
     call ftgkyj(unit,'TFIELDS',columns,comment,status)
     if(status .ne. 0) stop 'Cannot determine No of columns'

     !Copy chain file to new chain file
     status = 0
     call copyhdu(unit,newchainfile,newunit)
     if( status .ne. 0 ) write(*,*)"Could not copy chain file"
     write(*,*)"unit=",unit
     write(*,*)"newunit=",newunit
     
     !Add a new column to the new CHAIN binary table
     !First move to the CHAIN extension
     status = 0
     call ftmnhd(newunit,2,'CHAIN',0,status)
     if( status .ne. 0 ) write(*,*)"Couldn't move to CHAIN extension of newchainfile"
     !Insert the name and type of the new column.
     status = 0
     colnum = columns + 1
     ncols  = 1
     ttype(1) = 'Dkpc'
     tform(1) = '1D'
     tunit(1) = 'kpc'
     call FTICLS(newunit,colnum,ncols,ttype,tform,status)
     if( status .ne. 0 ) write(*,*)"FTICLS fucked up"
  
     !Go through each step and calculate distance for each one
     write(*,*)"Calculating distance and appending to chain..."
     do k = 1,steps
        status  = 0
        !Read in parameters
        call ftgcve(unit,columns,k,1,1,-1.0,chisquared,anynull,status)
        do i = 1,15
           if( iparam(i) .gt. 0 )then
              call ftgcve(unit,iparam(i),k,1,1,-1.0,parray(i),anynull,status)
           else
              parray(i) = param(i)
           end if
        end do
        !Calculate distance
        dist = distance(parray,xe,verbose)        
        !Append to level 1pt 5 event list
        call ftpcld(newunit,columns+1,k,1,1,dist,status)
     end do

  end if
  
! Close chain files
  call ftclos(unit, status)
  call ftfiou(unit, status)
  call ftclos(newunit, status)
  call ftfiou(newunit, status)
  
  ! write(*,*)"fort.96: logxi(r) vs r"
  
end program post_process

  
!-----------------------------------------------------------------------
function distance(param,xe,verbose)
  implicit none
  !Output
  double precision distance
  !Inputs
  real param(15)
  integer xe,verbose
  !Internal
  integer          :: refvar,Cp,dset,nex,nlp,i
  parameter (nex=2**12,nlp=1)
  real             :: par(32),dloge,earx(0:nex)
  real, parameter  :: Emin = 1e-2, Emax = 3e3
  double precision, parameter :: pi = acos(-1.d0), rnmax = 300.d0
  double precision :: a, inc, rin, rout, zcos, Gamma, h(nlp)
  real             :: logxi, Dkpc, Afe, lognep, Ecut_obs
  double precision :: eta_0, eta, honr
  real             :: beta_p, Nh, boost, Mass, Anorm
  double precision :: qboost,b1,b2,muobs
  real             :: floHz, fhiHz, DelA, DelAB(nlp), g(nlp), Ecut_s
  integer          :: ReIm, resp_matr, ndelta, npts(nlp), m
  parameter (ndelta=1000)
  double precision :: rlp(ndelta,nlp),dcosdr(ndelta,nlp),tlp(ndelta,nlp)
  double precision :: cosd(ndelta,nlp),cosdout(nlp),mudisk,rh,rmin,disco
  double precision :: gso(nlp),lens(nlp),tauso(nlp),logxip
  double precision :: cosdelta_obs(nlp),fcons,contx_int(nlp)
  integer          :: Cp_cont
  real             :: contx(nex,nlp),E,dE
  double precision :: logxir(xe),gsdr(xe),logner(xe),logxieff(xe)
  double precision :: pnorm,pnormer,re,dlogxi
  
! Hardwired setting (makes no difference)
  refvar   = 1
    
! Create general parameter array
  call dummy_reltransDCp(param,par,Cp,dset)

! Set up logarithmic energy grid
  dloge = log10( Emax / Emin ) / real(nex)
  do i = 0, nex
     earx(i) = Emin * (Emax/Emin)**(real(i)/real(nex))
  end do
  
! Set the parameters
  call set_param(dset,par,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut_obs,&
           eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
           g,Anorm,resp_matr,refvar,verbose)        
  muobs  = cos( inc * pi / 180.d0 )
  mudisk = honr / sqrt( honr**2 + 1.d0  )
  logxip = logxi
  
! Set minimum r (ISCO) and convert rin and h to rg
  if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
  rmin   = disco( a )
  if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
  rh     = 1.d0+sqrt(1.d0-a**2)  
  if( rin .lt. rmin )then
     write(*,*)"Warning! rin<ISCO! Set to ISCO"
     rin = rmin
  end if
  do m=1,nlp 
     if( h(m) .lt. 0.d0 ) h(m) = abs(h(m)) * rh
     if( h(m) .lt. 1.5d0*rh )then
        write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
        h(m) = 1.5d0 * rh
     end if 
  end do
  
! Calculate functions required for emissivity
! Calculate dcos/dr and time lags vs r for the lamppost model
  call getdcos(a,h,mudisk,ndelta,nlp,rout,npts,rlp,dcosdr,tlp,cosd,cosdout) 
  
! Calculate integral of continuum spectrum plus relative quantities (cutoff energies, lensing/gfactors, luminosity, etc)
  call init_cont(nlp,a,h,zcos,Ecut_s,Ecut_obs,gso,muobs,lens,tauso,cosdelta_obs,Cp_cont,Cp,fcons,Gamma,&
         Dkpc,Mass,earx,Emin,Emax,contx,dlogE,verbose,dset,Anorm,contx_int,eta)
  
! Call reltransDCp logxi(r) subroutine
  call alt_radfunctions_dens(verbose,xe,rin,rnmax,eta_0,logxip,dble(lognep),a,h,Gamma,honr,rlp,dcosdr&
       &,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner)
  logxir = logxir + log10(boost)

! Setup the rtdist version by integrating the spectrum
  Dkpc   = 1.0
  dset   = 1
  call init_cont(nlp,a,h,zcos,Ecut_s,Ecut_obs,gso,muobs,lens,tauso,cosdelta_obs,Cp_cont,Cp,fcons,Gamma,&
       Dkpc,Mass,earx,Emin,Emax,contx,dlogE,verbose,dset,Anorm,contx_int,eta)

! Call rtdist logxi(r) subroutine
  pnorm = pnormer(b1,b2,qboost)
  call alt_radfuncs_dist(xe, rin, rnmax, b1, b2, qboost, fcons,&
     & dble(lognep), a, h(1), honr, rlp, dcosdr, cosd, ndelta, rmin, npts(1),&
     & logxieff, gsdr, logner, pnorm)

! Write both logxi(r) functions out
  ! write(96,*)"skip on"
  ! do i = 1,xe
  !    re     = (rnmax/rin)**(real(i-1) / real(xe))
  !    re     = re + (rnmax/rin)**(real(i) / real(xe))
  !    re     = re * rin * 0.5
  !    write(96,*)re,logxir(i),logxieff(i),logxieff(i)-logxir(i)
  ! end do
  ! write(96,*)"log x"
  ! write(96,*)"la y log\gc"
  ! write(96,*)"la x r (r\dg\u)"
  ! write(96,*)"lw 5"
  ! write(96,*)"cs 1.5"
  
  dlogxi = logxir(2) - logxieff(2)
  ! write(*,*)"dlogxi = ",dlogxi
  
  Distance = 10**( 0.5*dlogxi ) / boost

  return
end function distance
!-----------------------------------------------------------------------





!-----------------------------------------------------------------------
subroutine copyhdu(inunit,outfile,outunit)
! Inputs
! inunit:  unit number of (alreay opened) input file.
! outfile: name of output file to be created
! Outputs
! outunit: unit number of output file
  implicit none
  integer status,inunit,outunit,blocksize,morekeys,hdutype
  character (len=500) outfile
  
! The STATUS parameter must always be initialized.
  status=0

! Delete the file if it already exists, so we can then recreate it
! The deletefile subroutine is listed at the end of this file.
  call deletefile(outfile,status)

! Get  unused Logical Unit Numbers to use to open the FITS file.
  call ftgiou(outunit,status)

! Create the new empty FITS file (value of blocksize is ignored)
  blocksize=1
  call ftinit(outunit,outfile,blocksize,status)

! Skip to the 2nd extension in the input file
  call ftmahd(inunit,2,hdutype,status)
  
! FTCOPY copies the current HDU from the input FITS file to the output
! file.  The MOREKEY parameter allows one to reserve space for additional
! header keywords when the HDU is created.   FITSIO will automatically
! insert more header space if required, so programmers do not have to
! reserve space ahead of time, although it is more efficient to do so if
! it is known that more keywords will be appended to the header.
  morekeys=0
  call ftcopy(inunit,outunit,morekeys,status)

 end subroutine copyhdu
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
subroutine alt_radfuncs_dist(xe, rin, rnmax, b1, b2, qboost, fcons,&
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
  
! Now loop through xe radial bins
  do i = 1,xe
     !Radius
     re     = (rnmax/rin)**(real(i-1) / real(xe))
     re     = re + (rnmax/rin)**(real(i) / real(xe))
     re     = re * rin * 0.5
     re1(i) = re
     !Density
     logner(i) = adensity * zA_logne(re,rin,lognep) + (1-adensity) * lognep
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

  ! write(*,*)"max(logxieff)=",maxval(logxieff)
  
!check max and min for both ionisation and density
  ! logxieff = max( logxieff , 0.d0  )
  ! logxieff = min( logxieff , 4.7d0 )
  ! logner   = max( logner , 15.d0  )
  ! logner   = min( logner , 22.d0 )
  !...no need to enforce limits on logne since this is done in myreflect()
  !This is needed because reflionx has a different maximum to xillverDCp

  ! verbose = myenv("REV_VERB",0)
  ! if( verbose .gt. 2 )then
  !    !Write out logxir for plots
  !    lximax = -huge(lximax)
  !    do i = 1,xe
  !       write(188,*)re1(i),Fx(i),logxir(i)
  !       lximax = max( lximax , logxieff(i) )
  !    end do
  !    write(188,*)"no no"
  !    write(*,*)"MAX LOGXIeff = ",lximax
  ! end if
  
  return
end subroutine alt_radfuncs_dist  
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine alt_radfunctions_dens(verbose,xe,rin,rnmax,eta_0,logxip,lognep,spin,h,Gamma,honr,rlp,dcosdr&
     &,cosd,contx_int,ndelta,nlp,rmin,npts,logxir,gsdr,logner)
    ! In  : xe,rin,rnmax,eta_0,logxip,spin,h,honr,rlp,dcosdr,cosd,ndelta,rmin,npts
    ! Out : logxir(1:xe), gsdr(1:xe), logner(1:xe)
    use env_variables
    implicit none
    integer         , intent(IN)   :: xe, ndelta, nlp, npts(nlp)
    double precision, intent(IN)   :: rin,rmin,rnmax,eta_0,logxip,lognep,spin,h(nlp),honr,Gamma
    real                           :: gso(nlp)
    double precision, intent(IN)   :: rlp(ndelta,nlp), dcosdr(ndelta,nlp), cosd(ndelta,nlp), contx_int(nlp)
    double precision :: rlp_column(ndelta),dcosdr_column(ndelta),cosd_column(ndelta), dgsofac
    double precision, intent(INOUT):: logxir(xe), gsdr(xe), logner(xe)
    integer          :: i, kk, get_index, myenv, l, m, verbose
    double precision :: rp, logxinorm, lognenorm,  mus, interper, newtex, mui, dinang, gsd(nlp), dglpfacthick
    double precision :: xi_lp(xe,nlp), logxi_lp(xe,nlp), logxip_lp(nlp), xitot, xiraw, mylogne, mudisk, gsd_temp
    double precision, allocatable :: rad(:)

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
            if (m .eq. 2) xi_lp(i,m) = eta_0*xi_lp(i,m)
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
    
    ! do m=1,nlp 
    !     do i=1,xe
    !         logxi_lp(i,m) = log10(xi_lp(i,m)) - logner(i) - lognenorm - logxinorm + lognep + logxip        
    !     end do
    !     logxip_lp(m) = max(maxval(logxi_lp(:,m)),0.)
    ! end do
    
    ! !Write radii, ionisation (for both and each LP), gamma factors, and log(xi(r))+log(ne(r)) (which is nearly the same as
    ! !epsilon(r) for identical coronal spectrra and gamma=2) to file. 
    ! !note 1) we need to do this before the ionisation array is set to have a minimum of 0, in order
    ! !to recover the correct scaling of the emissivity at large radii
    ! !2) in order to correctly compare the dfer_arr array with the single LP case, it has to be renormalized by (1+eta_0)
    ! if( verbose .gt. 1 ) then
    !     print*, "Peak ionisations from each LP: first " , logxip_lp(1), " second ", logxip_lp(2)
    !     open (unit = 27, file = 'Output/RadialScalings.dat', status='replace', action = 'write')
    !         do i = 1, xe
    !             write(27,*) rad(i), logxir(i), gsdr(i), logxir(i)+logner(i), logxi_lp(i,1), logxi_lp(i,2), dfer_arr(i) 
    !         end do 
    !     close(27)    
    ! end if
    
    ! !check max and min for ionisation 
    ! logxir = max( logxir , 0.d0  )
    ! logxir = min( logxir , 4.7d0 )
    
    deallocate(rad)

    return
end subroutine alt_radfunctions_dens
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
subroutine dummy_reltransDCp(param,par,Cp,dset)
  implicit none
  integer :: Cp, dset
  real    :: param(15),par(32)
! Settings
  Cp   = 2   !|Cp|=2 means nthcomp, Cp>1 means there is a density parameter     
  dset = 0   !dset=0 means distance is not set, logxi set instead
! Transfer to general parameter array
  par(1)  = param(1)         !h1
  par(2)  = 0.               !h2
  par(3)  = param(2)         !a
  par(4)  = param(3)         !inc
  par(5)  = param(4)         !rin
  par(6)  = param(5)         !rout
  par(7)  = param(6)         !zcos
  par(8)  = param(7)         !Gamma
  par(9)  = param(8)         !logxi
  par(10) = param(9)         !Afe
  par(11) = param(10)        !lognep
  par(12) = param(11)        !kTe
  par(13) = 0.               !eta_0
  par(14) = 0.               !eta
  par(15) = 0.               !beta_p
  par(16) = param(12)        !Nh
  par(17) = param(13)        !boost
  par(18) = 1.0              !qboost
  par(19) = param(14)        !Mass
  par(20) = 0.0              !honr
  par(21) = 0.0              !b1
  par(22) = 0.0              !b2
  par(23) = 0.0              !floHz
  par(24) = 0.0              !fhiHz
  par(25) = 1.0              !ReIm
  par(26) = 0.0              !DelA
  par(27) = 0.0              !DelAB
  par(28) = 0.0              !g
  par(29) = 0.               !DelAB2
  par(30) = 0.               !g2
  par(31) = param(15)        !Anorm
  par(32) = 1.0              !telescope response  
  return
end subroutine dummy_reltransDCp
!-----------------------------------------------------------------------






!-----------------------------------------------------------------------
      subroutine deletefile(filename,status)

!C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

!C  Simply return if status is greater than zero
      if (status .gt. 0)return

!C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
!C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
!C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!C  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
!-----------------------------------------------------------------------
