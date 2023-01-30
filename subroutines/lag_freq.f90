subroutine lag_freq(nex,earx,nf,fix,flo,fhi,Emin,Emax,nlp,contx,absorbx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                    h,z,Gamma,eta,beta_p,boost,g,DelAB,ionvar,ReGraw,ImGraw)

! Calculates the energy-averaged cross spectrum in the energy bins of interest. REDO ALL THIS COMMENT IT'S WRONG
!
! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
! Output: Graw(E,\nu) after multiplying by the absorption model
! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
    use constants
    implicit none
    integer, intent(in) :: nex,nf,ionvar,nlp
    integer             :: Ea1,Ea2,Eb1,Eb2       
    real   , intent(in) :: g(nlp),DelAB(nlp),boost,z,Gamma,Emin,Emax,beta_p,eta
    real   , intent(in) :: gso(nlp),tauso(nlp),h(nlp)
    real                :: gslope,ABslope
    real   , intent(in) :: earx(0:nex),contx(nex,nlp),absorbx(nex),fix(0:nf)
    real                :: ReW0(nlp,nex,nf),ImW0(nlp,nex,nf),ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),&
                           ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)                       
    real,    intent(out):: ReGraw(nf),ImGraw(nf)
    real                :: ReGrawEa,ImGrawEa,ReGrawEb,ImGrawEb
    real                :: E,fac,TempReG,TempImG 
    real                :: f,DelAB_nu,g_nu
    real                :: tau_d,phase_d,tau_p,phase_p,beta,flo,fhi
    complex             :: W0,W1,W2,W3,Sraw,cexp_p,cexp_d,cexp_phi,Stemp
    integer             :: i,j,m

    call energy_bounds(nex,Emin,Emax,Ea1,Ea2,Eb1,Eb2)

    gslope = 1.
    ABslope = 1.
    
    if(nlp .gt. 1) then
        ReW0(2,:,:) = eta*ReW0(2,:,:)
        ImW0(2,:,:) = eta*ImW0(2,:,:)
        ReW1(2,:,:) = eta*ReW1(2,:,:)
        ImW1(2,:,:) = eta*ImW1(2,:,:)
        ReW2(2,:,:) = eta*ReW2(2,:,:)
        ImW2(2,:,:) = eta*ImW2(2,:,:)
        ReW3(2,:,:) = eta*ReW3(2,:,:)
        ImW3(2,:,:) = eta*ImW3(2,:,:)
    endif 
    
    !Now calculate the cross-spectrum (/complex covariance), including absorption
    do j = 1, nf
        f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )        
        ReGrawEa = 0.0
        ImGrawEa = 0.0
        ReGrawEb = 0.0
        ImGrawEb = 0.0             
        do i = Ea1, Ea2
            phase_d = 0.
            phase_p = 0.
            tau_d = 0.
            tau_p = 0.
            Sraw = 0.
            E = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                DelAB_nu = DelAB(m) * ((fix(1) + fix(0))*0.5/f)**ABslope
                g_nu = g(m) * ((fix(1) + fix(0))*0.5/f)**gslope
                fac = log(gso(m)/((1.0+z)*E))
                !set up phase factors
                if (m .gt. 1) then  
                    tau_d = (tauso(m)-tauso(1))
                    tau_p = (h(m) - h(1))/(beta_p)             
                    phase_d = 2.*pi*tau_d*f
                    phase_p = 2.*pi*tau_p*f
                endif  
                cexp_d = cmplx(cos(phase_d),sin(phase_d))
                cexp_p = cmplx(cos(phase_p),sin(phase_p)) 
                cexp_phi = cmplx(cos(DelAB_nu),sin(DelAB_nu))  
                !print*,m,cexp_d,cexp_p,cexp_phi,f*fconv           
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = g_nu*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m))
                Stemp = Stemp + W0 + W3 + cexp_d*contx(i,m)
                Stemp = cexp_p*Stemp
                Sraw = Sraw + Stemp
            end do
            !separate into real/imaginary parts for compatibility with the rest of the code
            ReGrawEa = ReGrawEa + real(Sraw)*absorbx(i)
            ImGrawEa = ImGrawEa + aimag(Sraw)*absorbx(i)
        end do

        do i = Eb1, Eb2
            phase_d = 0.
            phase_p = 0.
            tau_d = 0.
            tau_p = 0.
            Sraw = 0.      
            E = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                DelAB_nu = DelAB(m) * ((fix(1) + fix(0))*0.5/f)**ABslope
                g_nu = g(m) * ((fix(1) + fix(0))*0.5/f)**gslope
                fac = log(gso(m)/((1.0+z)*E))
                !set up phase factors
                if (m .gt. 1) then
                    tau_d = (tauso(m)-tauso(1))
                    tau_p = (h(m) - h(1))/(beta_p)
                    phase_d = 2.*pi*tau_d*f
                    phase_p = 2.*pi*tau_p*f
                endif    
                cexp_d = cmplx(cos(phase_d),sin(phase_d))
                cexp_p = cmplx(cos(phase_p),sin(phase_p)) 
                cexp_phi = cmplx(cos(DelAB_nu),sin(DelAB_nu))             
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = g_nu*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m))
                Stemp = Stemp + W0 + W3 + cexp_d*contx(i,m)
                Stemp = cexp_p*Stemp
                Sraw = Sraw + Stemp
            end do
            !separate into real/imaginary parts for compatibility with the rest of the code
            ReGrawEb = ReGrawEb + real(Sraw)*absorbx(i)
            ImGrawEb = ImGrawEb + aimag(Sraw)*absorbx(i)            
        end do
        !Now cross-spectrum between the two energy bands
        !note: here the conjugate is b 
        ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
        ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
    end do

    return
end subroutine lag_freq

subroutine lag_freq_nocoh(nex,earx,nf,fix,flo,fhi,Emin,Emax,nlp,contx,absorbx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,&
                          ReW3,ImW3,h,z,Gamma,eta,boost,g,DelAB,ionvar,ReGraw,ImGraw)
                                
    use constants
    implicit none
    integer, intent(in) :: nex,nf,ionvar,nlp
    integer             :: Ea1,Ea2,Eb1,Eb2       
    real   , intent(in) :: g(nlp),DelAB(nlp),boost,z,Gamma,Emin,Emax,eta
    real   , intent(in) :: gso(nlp),tauso(nlp),h(nlp)
    real                :: gslope,ABslope
    real   , intent(in) :: earx(0:nex),contx(nex,nlp),absorbx(nex),fix(0:nf)
    real                :: ReW0(nlp,nex,nf),ImW0(nlp,nex,nf),ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),&
                           ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)                       
    real,    intent(out):: ReGraw(nf),ImGraw(nf)
    real                :: ReGrawEa,ImGrawEa,ReGrawEb,ImGrawEb
    real                :: E,fac,TempReG,TempImG 
    real                :: f,DelAB_nu,g_nu
    real                :: tau_d,phase_d,flo,fhi
    complex             :: W0,W1,W2,W3,Sraw,cexp_d,cexp_phi,Stemp
    integer             :: i,j,m

    call energy_bounds(nex,Emin,Emax,Ea1,Ea2,Eb1,Eb2)

    gslope = 1.
    Abslope = 1.
    
    do m=1,nlp 
        do j=1,nf 
            f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )        
            ReGrawEa = 0.0
            ImGrawEa = 0.0
            ReGrawEb = 0.0
            ImGrawEb = 0.0             
            do i=Ea1,Ea2 
                phase_d = 0.
                tau_d = 0.
                E = 0.5 * ( earx(i) + earx(i-1) )
                DelAB_nu = DelAB(m) * ((fix(1) + fix(0))*0.5/f)**ABslope
                g_nu = g(m) * ((fix(1) + fix(0))*0.5/f)**gslope
                fac = log(gso(m)/((1.0+z)*E))
                !set up phase factors
                if (m .gt. 1) then  
                    tau_d = (tauso(m)-tauso(1))          
                    phase_d = 2.*pi*tau_d*f
                endif  
                cexp_d = cmplx(cos(phase_d),sin(phase_d))
                cexp_phi = cmplx(cos(DelAB_nu),sin(DelAB_nu))  
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = g_nu*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m)) + W0 + W3 + cexp_d*contx(i,m)
                !separate into real/imaginary parts for compatibility with the rest of the code
                ReGrawEa = ReGrawEa + real(Stemp)*absorbx(i)
                ImGrawEa = ImGrawEa + aimag(Stemp)*absorbx(i)
            end do
            
            do i=Eb1,Eb2 
                phase_d = 0.
                tau_d = 0.
                E = 0.5 * ( earx(i) + earx(i-1) )
                DelAB_nu = DelAB(m) * ((fix(1) + fix(0))*0.5/f)**ABslope
                g_nu = g(m) * ((fix(1) + fix(0))*0.5/f)**gslope
                fac = log(gso(m)/((1.0+z)*E))
                !set up phase factors
                if (m .gt. 1) then  
                    tau_d = (tauso(m)-tauso(1))          
                    phase_d = 2.*pi*tau_d*f
                endif  
                cexp_d = cmplx(cos(phase_d),sin(phase_d))
                cexp_phi = cmplx(cos(DelAB_nu),sin(DelAB_nu))  
                !print*,m,cexp_d,cexp_p,cexp_phi,f*fconv           
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = g_nu*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m)) + W0 + W3 + cexp_d*contx(i,m)
                !separate into real/imaginary parts for compatibility with the rest of the code
                ReGrawEb = ReGrawEb + real(Stemp)*absorbx(i)
                ImGrawEb = ImGrawEb + aimag(Stemp)*absorbx(i)
            end do
            if(m .eq. 1) then
                ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
                ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
            else
                ReGraw(j) = ReGraw(j) + eta**2. * ((ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb))
                ImGraw(j) = ImGraw(j) + eta**2. * ((ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb))
            end if                        
        end do
    end do

    return
end subroutine lag_freq_nocoh

subroutine energy_bounds(nex,Emin,Emax,Ea1,Ea2,Eb1,Eb2)
    use telematrix 
    implicit none
    integer, intent(in) :: nex
    integer, intent(out):: Ea1,Ea2,Eb1,Eb2 
    real, intent(in)    :: Emin,Emax
    real                :: band1_Elo,band1_Ehi,band2_Elo,band2_Ehi
     
    if( needchans ) then
        write(*,*)"Enter lower energy in the first band"
        read(*,*) band1_Elo
        write(*,*)"Enter upper energy in the first band"
        read(*,*) band1_Ehi
        write(*,*)"Enter lower energy in the second band"
        read(*,*) band2_Elo
        write(*,*)"Enter upper energy in the second band"
        read(*,*) band2_Ehi
        Ea1 = ceiling( real(nex) * log10(band1_Elo / Emin) / log10(Emax / Emin))
        Ea2 = ceiling( real(nex) * log10(band1_Ehi / Emin) / log10(Emax / Emin))
        Eb1 = ceiling( real(nex) * log10(band2_Elo / Emin) / log10(Emax / Emin))
        Eb2 = ceiling( real(nex) * log10(band2_Ehi / Emin) / log10(Emax / Emin))
        needchans = .false.
    end if

    return
end subroutine energy_bounds
