!-----------------------------------------------------------------------
subroutine radfuncs_dist(config, model_args, fcons)
! Calculates ionization parameter, g-factor and disc density as a
! function of disc radius. This subroutine is only called for rtdist
!
! In  : config, model_args, fcons
! Out :
! logxir(1:xe) -- Effective ionisation as a function of r (via radial_grids)
! gsdr(1:xe)   -- Source-to-disc blueshift as a function of r (via radial_grids)
! logner(1:xe) -- Log10 of electron density as a function of r (via radial_grids)
!note: this does not work with multiple lamposts for now

  use common_types
  use dyn_gr, only: ndelta, rlp, dcosdr, cosd, npts
  use radial_grids, only: logxir, gsdr, logner, pnorm
  use env_variables
  implicit none
  type(t_config),          intent(in) :: config
  type(t_model_arguments), intent(in) :: model_args
  double precision,        intent(in) :: fcons
  integer          :: i, kk, get_index, get_env_int, verbose
  double precision :: re,re1(config%xe),zA_logne,cosfac,mus,interper,newtex,mudisk
  double precision, parameter :: pi = acos(-1.d0)
  double precision :: ptf,pfunc_raw,gsd,dglpfacthick,eps_bol,Fx(config%xe),logxir_raw(config%xe),mui,dinang
  double precision :: dareafac,lximax
  
! Set disk opening angle
  mudisk   = model_args%honr / sqrt( model_args%honr**2 + 1.d0  )
! Now loop through xe radial bins
  do i = 1,config%xe
     !Radius
     re     = (config%rnmax/model_args%rin)**(real(i-1) / real(config%xe))
     re     = re + (config%rnmax/model_args%rin)**(real(i) / real(config%xe))
     re     = re * model_args%rin * 0.5
     write(*,*)"re=",re
     re1(i) = re
     !Density
     logner(i) = zA_logne(re,model_args%rin,dble(model_args%lognep))
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
     eps_bol = gsd**2 * 2.0 * pi * ptf
     eps_bol = eps_bol * cosfac / dareafac(re,model_args%a)
     Fx(i)   = fcons * eps_bol
     !Calculate logxi(r)
     logxir_raw(i) = log10( 4.0 * pi * Fx(i) ) - logner(i)
     !Now adjust to effective ionization parameter
     mui       = dinang(model_args%a, re, model_args%h(1), mus)
     logxir(i) = logxir_raw(i) - 0.1505 - log10(mui)
  end do

!check max and min for both ionisation and density
  logxir = max( logxir , 0.d0  )
  logxir = min( logxir , 4.7d0 )
  ! logner   = max( logner , 15.d0  )
  ! logner   = min( logner , 22.d0 )
  !...no need to enforce limits on logne since this is done in myreflect()
  !This is needed because reflionx has a different maximum to xillverDCp

  verbose = get_env_int("REV_VERB",0)
  if( verbose .gt. 2 )then
     !Write out logxir for plots
     lximax = -huge(lximax)
     do i = 1,config%xe
        write(188,*)re1(i),Fx(i),logxir_raw(i)
        lximax = max( lximax , logxir(i) )
     end do
     write(188,*)"no no"
     write(*,*)"MAX LOGXIeff = ",lximax
  end if

  return
end subroutine radfuncs_dist
!-----------------------------------------------------------------------
