subroutine rawS(nex,earx,nf,flo,fhi,nlp,contx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,h,z,Gamma,eta,&
                beta_p,boost,g,DelAB,ionvar,DC,ReSraw,ImSraw)
    ! Calculates the FT of the spectrum before multiplying by the absorption model
    !
    ! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
    ! Output: Sraw(E,\nu) before multiplying by the absorption model
    ! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
    !
    ! Inputs:
    ! nex,earx(0:nex)       Internal energy grid
    ! nf                    Number of frequency bins
    ! floHz/fhiHz           Frequency range - needed for double LP model
    ! contx(1:nex)          Continuum photar model
    ! tauso                 Source to observer lag - needed for double LP model
    ! ReW0(1:nex,1:nf)      Real part of transfer function number 0
    ! ImW0(1:nex,1:nf)      Imaginary part of transfer function number 0
    ! ReW1(1:nex,1:nf)      Real part of transfer function number 1
    ! ImW1(1:nex,1:nf)      Imaginary part of transfer function number 1
    ! ReW2(1:nex,1:nf)      Real part of transfer function number 2
    ! ImW2(1:nex,1:nf)      Imaginary part of transfer function number 2
    ! ReW3(1:nex,1:nf)      Real part of transfer function number 3
    ! ImW3(1:nex,1:nf)      Imaginary part of transfer function number 3
    ! g/g2                  The model parameter called gamma/gamma2
    ! DelAB/DelAB2/         The model parameter called phiAB/phiAB2
    ! boost                 The boost model parameter 
    ! z                     Cosmological redshift
    ! gso                   Blueshift travelling from source to observer
    ! Gamma                 Photon index
    ! eta                   Weighing factor for the two lamp posts for nlp > 1
    ! ionvar                Integer: ionvar=1 means include ionization variations, ionvar=0 means don't.
    ! DC                    Integer: DC=1 means this is the DC component, DC=0 means opposite.

    ! Outputs
    ! ReGraw(1:nex,1:nf)    Real part of Sraw(E,nu)      - in specific photon flux (photar/dE)
    ! ImGraw(1:nex,1:nf)    Imaginary part of Sraw(E,nu) - in specific photon flux (photar/dE)
    use constants
    implicit none
    integer nex,nf,ionvar,DC,nlp
    complex W0,W1,W2,W3,Sraw(nex,nf),cexp_p,cexp_d,cexp_phi,Stemp   
    real earx(0:nex),contx(nex,nlp),tauso(nlp),ReW0(nlp,nex,nf),ImW0(nlp,nex,nf)
    real ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)
    real DelAB(nlp),g(nlp),boost,z,gso(nlp),Gamma,eta,ReSraw(nex,nf),ImSraw(nex,nf),h(nlp),beta_p 
    real E,fac,tau_d,phase_d,tau_p,phase_p,f,flo,fhi
    integer i,j,m

    ReSraw = 0.
    ImSraw = 0.
    Sraw = 0.

    phase_d = 0.
    phase_p = 0.
    tau_d = 0.
    tau_p = 0.

    do m=1,nlp 
       if (boost .lt. 0 .and. DC .eq. 1) then
            do j = 1,nf
               do i = 1,nex
                  if (m .gt. 1) ReW0(m,i,j) = eta*ReW0(m,i,j)
                  ReSraw(i,j) = ReSraw(i,j) + (-boost) * ReW0(m,i,j)
                enddo
            enddo  
        else
            if( m .gt. 1 ) then
                !set up extra terms if second lamp post present
                ReW0(m,:,:) = eta*ReW0(m,:,:)
                ImW0(m,:,:) = eta*ImW0(m,:,:)
                ReW1(m,:,:) = eta*ReW1(m,:,:)
                ImW1(m,:,:) = eta*ImW1(m,:,:)
                ReW2(m,:,:) = eta*ReW2(m,:,:)
                ImW2(m,:,:) = eta*ImW2(m,:,:)
                ReW3(m,:,:) = eta*ReW3(m,:,:)
                ImW3(m,:,:) = eta*ImW3(m,:,:)
                tau_d = tauso(m)-tauso(1)
                tau_p = (h(m) - h(1))/(beta_p)
            end if
            do j = 1,nf
                if (DC .eq. 1) then
                    f = 0.
                else 
                    f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )
                endif
                do i = 1,nex
                    E   = 0.5 * ( earx(i) + earx(i-1) )
                    fac = log(gso(m)/((1.0+z)*E))
                    !set up phase factors
                    if (m .gt. 1) then
                        phase_d = 2.*pi*tau_d*f
                        phase_p = 2.*pi*tau_p*f
                    endif    
                    cexp_d = cmplx(cos(phase_d),sin(phase_d))
                    cexp_p = cmplx(cos(phase_p),sin(phase_p)) 
                    cexp_phi = cmplx(cos(DelAB(m)),sin(DelAB(m)))
                    !set up transfer functions 
                    W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                    W1 = (1-DC) * boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                    W2 = (1-DC) * boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                    W3 = ionvar * (1-DC) * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                    !calculate complex covariance
                    !note: the reason we use complex here is to ease the calculations 
                    !when we add all the extra phases from the double lamp post 
                    Stemp = g(m)*cexp_phi*(W1 + W2 + fac*cexp_d*contx(i,m))
                    Stemp = Stemp + W0 + W3 + cexp_d*contx(i,m)
                    Stemp = cexp_p*Stemp
                    Sraw(i,j) = Sraw(i,j) + Stemp
                    !separate into real/imaginary parts for compatibility with the rest of the code
                    ReSraw(i,j) = real(Sraw(i,j))
                    ImSraw(i,j) = aimag(Sraw(i,j))
                    ! write(74,*) E, fac
                    ! write(75,*) E, gso(m)
                    ! write(76,*) E, Stemp
                    ! write(77,*) E, Sraw(i,j)
                 enddo
            enddo 
        endif    
     end do

    return
end subroutine rawS
