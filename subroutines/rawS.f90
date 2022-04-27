subroutine rawS(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,&
    boost,z,gso,Gamma,ionvar,DC,ReGraw,ImGraw)
! Calculates the FT of the spectrum before multiplying by the absorption model
!
! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
! Output: Sraw(E,\nu) before multiplying by the absorption model
! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
!
! Inputs:
! nex,earx(0:nex)       Internal energy grid
! nf                    Number of frequency bins
! contx(1:nex)          Continuum photar model
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
! ionvar                Integer: ionvar=1 means include ionization variations, ionvar=0 means don't.
! DC                    Integer: DC=1 means this is the DC component, DC=0 means opposite.
  
! Outputs
! ReGraw(1:nex,1:nf)    Real part of Sraw(E,nu)      - in specific photon flux (photar/dE)
! ImGraw(1:nex,1:nf)    Imaginary part of Sraw(E,nu) - in specific photon flux (photar/dE)
    implicit none
    integer nex,nf,ionvar,DC,nlp
    real earx(0:nex),corr,contx(nex,nlp),ReW0(nex,nf),ImW0(nex,nf)
    real ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nex,nf),ImW3(nex,nf)
    real DelAB(nlp),g(nlp),boost,z,gso(nlp),Gamma,ReGraw(nex,nf),ImGraw(nex,nf)
    real sinD,cosD,E,fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s,gsoz
    integer i,j,m

    ReGraw = 0.
    ImGraw = 0.

    !Now calculate the cross-spectrum (/complex covariance)
    if (boost .lt. 0 .and. DC .eq. 1) then 
        do j = 1,nf
            do i = 1,nex
                ReGraw(i,j) = (-boost) * ReW0(i,j)
            enddo
        enddo     
    else
        do j = 1,nf
            do i = 1,nex
                !Multiply by boost parameter and group like terms
                ReW0s = boost * ReW0(i,j)
                ImW0s = boost * ImW0(i,j)                
                ReW3s = (1-DC) * ionvar * boost * ReW3(i,j)
                ImW3s = (1-DC) * ionvar * boost * ImW3(i,j)  
                E   = 0.5 * ( earx(i) + earx(i-1) )
                do m=1,nlp
                    fac = log(gso(m)/((1.0+z)*E))
                    sinD = sin(DelAB(m))
                    cosD = cos(DelAB(m)) 
                    ReWbs = (1-DC) * boost * (ReW1(m,i,j) + ReW2(m,i,j))
                    ImWbs = (1-DC) * boost * (ImW1(m,i,j) + ImW2(m,i,j))
                    !Real part contribution for each LP
                    ReGraw(i,j) = ReGraw(i,j) + g(m) * (cosD*(fac*contx(i,m) + ReWbs) - sinD*ImWbs) 
                    ReGraw(i,j) = ReGraw(i,j) + contx(i,m)         
                    !Imaginary part contribution for each LP 
                    ImGraw(i,j) = ImGraw(i,j) + g(m) * (sinD*(fac*contx(i,m) + ReWbs) + cosD*ImWbs) 
                end do 
                !Real part shared by both
                ReGraw(i,j) = ReGraw(i,j) + ReW0s + ReW3s
                !Imaginary part shared by both
                ImGraw(i,j) = ImGraw(i,j) + ImW0s + ImW3s
            end do
        end do
    endif
    

    return
end subroutine rawS
