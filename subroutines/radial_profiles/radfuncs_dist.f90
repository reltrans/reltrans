!-----------------------------------------------------------------------
subroutine radfuncs_dist(config, model_args, fcons)
! Calculates ionization parameter, g-factor and disc density as a
! function of disc radius. This subroutine is only called for rtdist
!
! In  : config, model_args, fcons
! Out :
! radial_grids arrays, indexed by 1:ion_zones:
! logxir: effective ionisation; gsdr: source-to-disc blueshift
! logner: log10 of electron density
!note: this does not work with multiple lamposts for now

  use rtconstants, only: wp
  use common_types
  use dyn_gr, only: ndelta, rlp, dcosdr, cosd, npts
  use radial_grids, only: logxir, gsdr, logner, pnorm, logxir_natural
  use env_variables, only : adensity
  implicit none
  type(t_config),          intent(in) :: config
  type(t_model_arguments), intent(in) :: model_args
  real(wp),        intent(in) :: fcons
  integer          :: i, kk, get_index, get_env_int, verbose
  real(wp) :: re,re1(config%ion_zones),zA_logne,cosfac,mus,interper,           &
      newtex,mudisk
  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp) :: ptf,pfunc_raw,gsd,dglpfacthick,eps_bol,Fx(config%ion_zones),    &
      logxir_raw(config%ion_zones),mui,dinang
  real(wp) :: dareafac,lximax
  
! Set disk opening angle
  mudisk   = model_args%honr / sqrt( model_args%honr**2 + 1.0_wp  )
! Now loop through ion_zones radial bins
  do i = 1,config%ion_zones
     !Radius
     re     = (config%rnmax/model_args%rin)**                               &
         (real(i-1, wp) / real(config%ion_zones, wp))
     re     = re + (config%rnmax/model_args%rin)**                          &
         (real(i, wp) / real(config%ion_zones, wp))
     re     = re * model_args%rin * 0.5_wp
     re1(i) = re
     !Density
     !IF adensity = 0 logner(i)=dble(model_args%lognep) for every i
     !IF adensity = 1 logner(i)=zA_logne(re,model_args%rin,dble(model_args%lognep) where re depends on i
     logner(i) = abs(model_args%lognep * (adensity - 1))
     logner(i) = logner(i) + (adensity * zA_logne(re,model_args%rin,model_args%lognep))
     !Interpolate functions from rpl grid
     kk     = get_index(rlp,ndelta,re,config%rmin,npts(1))
     cosfac = interper(rlp,dcosdr,ndelta,re,kk)
     mus    = interper(rlp,cosd  ,ndelta,re,kk)
     if( kk .eq. npts(1) )then !Extrapolation onto Newtonian grid
        cosfac = newtex(rlp,dcosdr,ndelta,re,model_args%h(1),model_args%honr,kk)
        mus    = newtex(rlp,cosd  ,ndelta,re,model_args%h(1),model_args%honr,kk)
     end if
     !Calculate 13.6 eV - 13.6 keV illuminating flux
     ptf     = pnorm * pfunc_raw(-mus,model_args%b1,model_args%b2,model_args%qboost)
     gsd     = dglpfacthick(re,model_args%a,model_args%h(1),mudisk)
     !gsd     = dglpfac(re,model_args%a,model_args%h(1))
     gsdr(i) = gsd
     eps_bol = gsd**2 * 2.0_wp * pi * ptf
     eps_bol = eps_bol * cosfac / dareafac(re,model_args%a)
     Fx(i)   = fcons * eps_bol
     !Calculate logxi(r)
     logxir_raw(i) = log10( 4.0_wp * pi * Fx(i) ) - logner(i)
     !Now adjust to effective ionization parameter
     mui       = dinang(model_args%a, re, model_args%h(1), mus)
     logxir(i) = logxir_raw(i) - 0.1505_wp - log10(mui)
  end do

  ! logxir_natural: unmodified (natural) ionization parameter before clipping
  ! logxir: effective ionization parameter after clipping
  logxir_natural = logxir  ! Need for the post processing code
  
  !check max and min for both ionisation and density
  logxir = max( logxir , 0.0_wp  )
  logxir = min( logxir , 4.7_wp )

  verbose = get_env_int("REV_VERB",0)
  if( verbose .gt. 2 )then
     !Write out logxir for plots
     lximax = -huge(lximax)
     do i = 1,config%ion_zones
        lximax = max( lximax , logxir_natural(i) )
     end do
     write(*,*)"MAX LOGXIeff = ",lximax
  end if

  return
end subroutine radfuncs_dist
!-----------------------------------------------------------------------
