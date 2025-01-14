subroutine init_cont(nlp, a, h, zcos, Ecut_s, Ecut_obs, logxi, logne, &
     muobs, Cp_cont, Cp, fcons, Gamma, Dkpc, Mass,&
    earx, Emin, Emax, contx, dlogE, verbose, dset, Anorm, contx_int, eta)
    !!!sets up the continuum arrays/quantities depending on model parameters/flavour 
    use dyn_gr
    use conv_mod
    use gr_continuum
    implicit none
    integer         , intent(in)    :: nlp,Cp,dset,verbose
    real            , intent(in)    :: Dkpc,Anorm,Mass,dlogE,Emin,Emax, logxi, logne
    double precision, intent(in)    :: h(nlp)
    double precision, intent(in)    :: a,zcos,Gamma,muobs,eta
    integer         , intent(out)   :: Cp_cont
    real            , intent(out)   :: earx(0:nex),contx(nex,nlp)
    double precision, intent(out)   :: fcons,contx_int(nlp)
    
    integer                         :: m, i
    real                            :: Ecut_s,Ecut_obs,Eintegrate
    double precision                :: lacc,ell13pt6,get_lacc,get_fcons,dgsofac

    
    if (nlp .eq. 1) then 

       ! gso(1) = real( dgsofac(a,h(1)) ) 
       ! call getlens(a,h(1),muobs,lens(1),tauso(1),cosdelta_obs(1))
       ! if( tauso(1) .ne. tauso(1) ) stop "tauso is NaN"
       
       Ecut_s = real(1.d0+zcos) * Ecut_obs / gso(1)
       Cp_cont = Cp
       if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx
       call getcont(Cp, earx, nex, Gamma, Ecut_obs, logxi, logne, contx(:,1))
       
       if( dset .eq. 1 ) then
          fcons = get_fcons(h(1),a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)     
       else
          fcons = 0.0
       end if
         
       if( verbose .gt. 0 )then
          if( dset .eq. 1 )then    
             lacc = get_lacc(h(1),a,zcos,Gamma,Dkpc,Mass,Anorm,nex,earx,contx,dlogE)
             write(*,*)"Lacc/Ledd=",lacc 
             ell13pt6 = fcons * Mass * 1.73152e-28
             write(*,*)"13.6eV-13.6keV luminosity of single source=",ell13pt6
          else
             call sourcelum(nex,earx,contx,real(Mass),real(gso(1)),real(Gamma))
          end if
          if( abs(Cp) .eq. 1 )then
             write(*,*)"Ecut in source restframe (keV)=",Ecut_s
          else
             write(*,*)"kTe in source restframe (keV)=", Ecut_s
          end if
          write(*,*) 'gso factor ', gso(1)
          write(*,*) 'lensing factor ', lens(1)

       end if
       
       ! do i = 1, nex
       !    write(61,*) (earx(i-1)+earx(i))*0.5 , contx(i,1)
       ! enddo
       contx_int(1) = 1. !note: for a single LP we don't need to account for this factor in the ionisation profile, so it's defaulted to 1       

       if (Cp .eq. 2) then
          ! write(*,*) 'nthcomp illumination'
          ! contx = lens(1) * (gso(1)/(real(1.d0+zcos)))**Gamma * contx
          contx = lens(1) * (gso(1)/(real(1.d0+zcos))) * contx
       else
          ! write(*,*) 'powerlaw illumination'
          contx = lens(1) * (gso(1)/(real(1.d0+zcos)))**Gamma * contx
       endif
    else
       do m=1,nlp   
          !here the observed cutoffs are set from the temperature in the source frame   
          ! gso(m) = real( dgsofac(a,h(m)) )
          ! call getlens(a,h(m),muobs,lens(m),tauso(m),cosdelta_obs(m))
          ! if( tauso(m) .ne. tauso(m) ) stop "tauso is NaN"

          Ecut_obs = Ecut_s * gso(m) / real(1.d0+zcos)
          Cp_cont = Cp 
          if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx        
          call getcont(Cp, earx, nex, Gamma, Ecut_obs, logxi, logne, contx(:,m))
          if (m .gt. 1) contx(:,m) = eta*contx(:,m)  
          !TODO fix this section, calculate luminosities better
          if( verbose .gt. 0 )then
             call sourcelum(nex,earx,contx(:,m),real(mass),real(gso(m)),real(Gamma))
             if( abs(Cp) .eq. 1 )then
                write(*,*)"Ecut observed from source #", m, "is (keV)=" ,Ecut_obs
             else
                write(*,*)"kTe observed from source #", m, "is (keV)=" ,Ecut_obs
             end if
          end if
          contx_int(m) = Eintegrate(Emin,Emax,nex,earx,contx(:,m),dlogE)    

          if (Cp .eq. 2) then
             contx = lens(1) * (gso(1)/(real(1.d0+zcos))) * contx
          else
             contx(:,m) = lens(m) * (gso(m)/(real(1.d0+zcos)))**Gamma * contx(:,m)
          endif

          ! do i = 1, nex
          !    write(10,*) (earx(i-1)+earx(i))*0.5, contx(i,m)
          ! enddo
          ! write(10,*) 'no no'
       end do
    end if  
    !TBD ADD PROPAGATION LAG HERE

end subroutine init_cont
