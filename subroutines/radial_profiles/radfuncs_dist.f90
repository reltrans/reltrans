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

  use env_variables
  implicit none
  integer         , intent(IN)   :: xe, ndelta, npts 
  double precision, intent(IN)   :: rin, rmin, rnmax, b1, b2, qboost
  double precision, intent(IN)   :: fcons, lognep, spin, h, honr
  double precision, intent(IN)   :: rlp(ndelta), dcosdr(ndelta), cosd(ndelta)
  double precision, intent(INOUT):: logxieff(xe), gsdr(xe), logner(xe)
  integer          :: i, kk, get_index, get_env_int, verbose
  double precision :: pnorm,re,re1(xe),zA_logne,cosfac,mus,interper,newtex,mudisk
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: ptf,pfunc_raw,gsd,dglpfacthick,eps_bol,Fx(xe),logxir(xe),mui,dinang
  double precision :: dareafac,lximax

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

  verbose = get_env_int("REV_VERB",0)
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
