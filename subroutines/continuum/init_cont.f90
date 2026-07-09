subroutine init_cont(config, model_args, arrays, Cp_cont, fcons, dset)
    !!!sets up the continuum arrays/quantities depending on model parameters/flavour
    use common_types
    use dyn_gr
    use conv_mod
    use gr_continuum
    implicit none
    type(t_config)         , intent(in)      :: config
    type(t_model_arguments), intent(inout)   :: model_args
    type(t_arrays)         , intent(inout)   :: arrays
    integer                , intent(in)      :: dset
    integer                , intent(out)     :: Cp_cont
    double precision       , intent(out)     :: fcons
    
    integer :: m
    real :: Cutoff_s, Cutoff_obs, Eintegrate
    double precision :: lacc, ell13pt6, get_lacc, get_fcons
    
    Cutoff_s   = model_args%Cutoff_s
    Cutoff_obs = model_args%Cutoff_obs
    
    if (model_args%nlp .eq. 1) then
       arrays%contx_int(1) = 1. !note: for a single LP we don't need to account for this factor in the ionisation profile, so it's defaulted to 1

       if( model_args%Cp .ge. 0 ) then
          ! write(*,*) 'nthcomp illumination for nthcomp and reflionx'
          Cp_cont = 2 !This is needed since we can't use getcont(Cp,...) because in reflionx Cp = 0
          Cutoff_obs = Cutoff_s * gso(1) / real(1.d0 + model_args%zcos)

          call getcont(Cp_cont, arrays%earx, nex, model_args%Gamma,            &
              Cutoff_s, Cutoff_obs, model_args%logxi, model_args%lognep,       &
              model_args%zcos, arrays%contx(:,1))
          arrays%contx = lens(1) / real(1.d0 + model_args%zcos)**2             &
              * gso(1) / real(1.d0 + model_args%zcos) * arrays%contx
       else if (model_args%Cp .eq. -1) then
          ! write(*,*) 'powerlaw illumination'
          Cutoff_s = real(1.d0 + model_args%zcos) * Cutoff_obs / gso(1)
          call getcont(model_args%Cp, arrays%earx, nex, model_args%Gamma,      &
              Cutoff_s, Cutoff_obs, model_args%logxi, model_args%lognep,       &
              model_args%zcos, arrays%contx(:,1))
          arrays%contx = lens(1) / real(1.d0 + model_args%zcos)**2             &
              * (gso(1) / real(1.d0 + model_args%zcos))**model_args%Gamma      &
              * arrays%contx
       endif

       model_args%Cutoff_s   = Cutoff_s
       model_args%Cutoff_obs = Cutoff_obs
       
       if( dset .eq. 1 ) then
          fcons = get_fcons(model_args%h(1), model_args%a, model_args%zcos,    &
              model_args%Gamma, model_args%Dkpc, model_args%Mass,              &
              model_args%Anorm, nex, arrays%earx, arrays%contx, config%dloge)
       else
          fcons = 0.0
       end if
         
       if( config%verbose .gt. 0 )then
          if( dset .eq. 1 )then    
             lacc = get_lacc(model_args%h(1), model_args%a, model_args%zcos,   &
                 model_args%Gamma, model_args%Dkpc, model_args%Mass,           &
                 model_args%Anorm, nex, arrays%earx, arrays%contx,             &
                 config%dloge)
             write(*,*)"Lacc/Ledd=",lacc 
             ell13pt6 = fcons * model_args%Mass * 1.73152e-28
             write(*,*)"13.6eV-13.6keV luminosity of single source=",ell13pt6
          else
             call sourcelum(nex, arrays%earx, arrays%contx,                    &
                 real(model_args%Mass), real(gso(1)), real(model_args%Gamma))
          end if
          if( abs(model_args%Cp) .eq. 1 )then
             write(*,*)"Ecut in source restframe (keV)=",Cutoff_s
             write(*,*)"Ecut in observer restframe (keV)=",Cutoff_obs
          else
             write(*,*)"kTe in source restframe (keV)=", Cutoff_s
             write(*,*)"kTe in observer restframe (keV)=", Cutoff_obs
          end if
          write(*,*) 'gso factor ', gso(1)
          write(*,*) 'lensing factor ', lens(1)

       end if
       
    else
       do m=1,model_args%nlp
          !here the observed cutoffs are set from the temperature in the source frame
          Cutoff_obs = Cutoff_s * gso(m) / real(1.d0 + model_args%zcos)
          call getcont(model_args%Cp, arrays%earx, nex, model_args%Gamma,      &
              Cutoff_s, Cutoff_obs, model_args%logxi, model_args%lognep,       &
              model_args%zcos, arrays%contx(:,m))
          if (m .gt. 1) arrays%contx(:,m) = model_args%eta * arrays%contx(:,m)
          !TODO fix this section, calculate luminosities better
          if( config%verbose .gt. 0 )then
             call sourcelum(nex, arrays%earx, arrays%contx(:,m),               &
                 real(model_args%Mass), real(gso(m)), real(model_args%Gamma))
             if( abs(model_args%Cp) .eq. 1 )then
                write(*,*)"Ecut observed from source #", m, "is (keV)=" ,Cutoff_obs
             else
                write(*,*)"kTe observed from source #", m, "is (keV)=" ,Cutoff_obs
             end if
          end if
          arrays%contx_int(m) = Eintegrate(config%Emin, config%Emax, nex,      &
              arrays%earx, arrays%contx(:,m), config%dloge)
          arrays%contx(:,m) = lens(m) / real(1.d0 + model_args%zcos)**2        &
              * gso(m) / real(1.d0 + model_args%zcos) * arrays%contx(:,m)
       end do
    end if

end subroutine init_cont
