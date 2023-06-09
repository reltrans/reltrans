subroutine init_cont(nlp,a,h,zcos,Ecut_s,Ecut_obs,gso,muobs,lens,tauso,cosdelta_obs,Cp_cont,Cp,fcons,Gamma,Dkpc,Mass,&
    earxi,Emin,Emax,dlogE,verbose,dset,Anorm,contx_int,eta)
    !!!sets up the continuum arrays/quantities depending on model parameters/flavour 
    use dyn_gr
    use dyn_en 
    use conv_mod
    implicit none
    integer, intent(in)             :: nlp,Cp,dset,verbose
    real, intent(in)                :: Dkpc,Anorm,Mass,dlogE,Emin,Emax
    double precision, intent(in)    :: h(nlp)
    double precision, intent(in)    :: a,zcos,Gamma,muobs,eta
    
    integer                         :: m
    real                            :: Ecut_s,Ecut_obs,Eintegrate,earxi(0:nexi),contxi(nexi,nlp)
    double precision                :: lacc,ell13pt6,get_lacc,get_fcons,dgsofac

    integer, intent(out)            :: Cp_cont
    !real, intent(out)               :: !earx(0:nex)!,contx(nex,nlp)
    double precision, intent(out)   :: gso(nlp),lens(nlp),tauso(nlp),cosdelta_obs(nlp),fcons,contx_int(nlp)
    
    if (nlp .eq. 1) then 
        gso(1) = real( dgsofac(a,h(1)) ) 
        Ecut_s = real(1.d0+zcos) * Ecut_obs / gso(1)
        call getlens(a,h(1),muobs,lens(1),tauso(1),cosdelta_obs(1))
        if( tauso(1) .ne. tauso(1) ) stop "tauso is NaN"
        Cp_cont = Cp
        if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx
        call getcont(nexi,earxi,Gamma,Ecut_obs,Cp_cont,contxi)   
        call rebinE(earxi,contxi,nexi,earx,contx,nex)       
        if( dset .eq. 1 )then
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
        end if         
        contx_int(1) = 1. !note: for a single LP we don't need to account for this factor in the ionisation profile, so it's defaulted to 1
        contx = lens(1) * (gso(1)/(real(1.d0+zcos)))**Gamma * contx        
    else 
        do m=1,nlp   
            !here the observed cutoffs are set from the temperature in the source frame   
            gso(m) = real( dgsofac(a,h(m)) )
            call getlens(a,h(m),muobs,lens(m),tauso(m),cosdelta_obs(m))
            if( tauso(m) .ne. tauso(m) ) stop "tauso is NaN"
            Ecut_obs = Ecut_s * gso(m) / real(1.d0+zcos)
            Cp_cont = Cp 
            if( Cp .eq. 0 ) Cp_cont = 2 !For reflection given by reflionx        
            call getcont(nexi,earxi,Gamma,Ecut_obs,Cp_cont,contxi(:,m))
            call rebinE(earxi,contxi(:,m),nexi,earx,contx(:,m),nex)    
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
            contx(:,m) = lens(m) * (gso(m)/(real(1.d0+zcos)))**Gamma * contx(:,m)   
        end do  
    end if  

end subroutine init_cont

!-----------------------------------------------------------------------
subroutine getcont(nex,earx,Gamma,Ecut_obs,Cp,contx)
    !!! Calculates continuum spectrum and sets xillver parameters !!!
    !!!  Arg:
    !  earx: energy grid
    !  nex: number of grid points
    !  Gamma: continuum spectrum inclination
    !  Afe: iron abundance (not important for the continuum)
    !  Ecut_obs: high energy cut-off or electron temperature
    !  dens: log(density) of the disc
    !  logxi: ionisation parameter (not important for the continuum)
    !  Cp: sets which xillver spectrum
    ! (output) contx: continuum spectrum
    !!! Last change: Adam 2021 March; based on Gullo's getcont code from 2020 Jul
    !
    implicit none
    integer, intent(in)  :: nex, Cp
    real   , intent(in)  :: earx(0:nex), Ecut_obs
    real   , intent(out) :: contx(nex)
    real                 :: xillpar(7), xillparDCp(8)
    double precision , intent(in) :: Gamma
    integer :: ifl   
    !Set the parameters for the continuum and not only 

    !First determine which xillver and so which array of parameters 7 or 8
    ! Remember that only xillverDCp has 8 par all the others have 7 but with different interpretation of parameter 3 (xillver(3)) 
    !to clean up the code and not pass a million arguments we just take default values for ne/csi/theta/Afe, it does not impact
    !the continuum anyway 
    if (Cp .eq. 2) then
        xillparDCp(1) = real( Gamma )
        xillparDCp(2) = 1.
        xillparDCp(3) = Ecut_obs
        xillparDCp(4) = 15.
        xillparDCp(5) = 0.0
        xillparDCp(6) = 0.0   
        xillparDCp(7) = 30.0 
        xillparDCp(8) = 0.0       !reflection fraction of 0
    else if (Cp .lt. 0 ) then
        xillpar(3) = Ecut_obs
    else if (Cp .eq. 1) then
        xillpar(3) = 15.
    endif

    xillpar(1) = real( Gamma )
    xillpar(2) = 1.
    xillpar(4) = 0.
    xillpar(5) = 0.0   !cosmological redshift is accounted for by the transfer function
    xillpar(6) = 30.0       !inclination angle (doesn't matter for continuum)
    xillpar(7) = 0.0       !reflection fraction of 0
    ifl        = 1

    call myxill(earx,nex,xillpar,xillparDCp,ifl,Cp,contx)
    return
end subroutine getcont
!-----------------------------------------------------------------------
