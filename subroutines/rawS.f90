subroutine rawS(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,&
     Gamma,ionvar,DC,ReGraw,ImGraw)
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
! g                     The model parameter called \gamma
! DelAB                 The model parameter to replace \phi_{B}. This one is better. In radians.
! boost                 The boost model parameter (named afac elsewhere in the code)
! z                     Cosmological redshift
! gso                   Blueshift travelling from source to observer
! lens                  Lensing factor
! Gamma                 Photon index
! ionvar                Integer: ionvar=1 means include ionization variations, ionvar=0 means don't.
! DC                    Integer: DC=1 means this is the DC component, DC=0 means opposite.
  
! Outputs
! ReGraw(1:nex,1:nf)    Real part of Sraw(E,nu)      - in specific photon flux (photar/dE)
! ImGraw(1:nex,1:nf)    Imaginary part of Sraw(E,nu) - in specific photon flux (photar/dE)
  implicit none
  integer nex,nf,ionvar,DC
  real earx(0:nex),corr,contx(nex),ReW0(nex,nf),ImW0(nex,nf)
  real ReW1(nex,nf),ImW1(nex,nf),ReW2(nex,nf),ImW2(nex,nf),ReW3(nex,nf),ImW3(nex,nf)
  real g,DelAB,boost,z,gso,lens,Gamma,ReGraw(nex,nf),ImGraw(nex,nf)
  real sinD,cosD,E,fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s,gsoz
  integer i,j
  sinD = sin(DelAB)
  cosD = cos(DelAB)
  gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
  corr = lens * gsoz**Gamma  !Correction factor for direct component
  !Now calculate the cross-spectrum (/complex covariance)
  do j = 1,nf
     do i = 1,nex
        E   = 0.5 * ( earx(i) + earx(i-1) )
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = (1-DC) * boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = (1-DC) * boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = (1-DC) * ionvar * boost * ReW3(i,j)
        ImW3s = (1-DC) * ionvar * boost * ImW3(i,j)
        !Real part
        ReGraw(i,j) = cosD * ( fac * corr * contx(i) + ReWbs )
        ReGraw(i,j) = ReGraw(i,j) - sinD * ImWbs
        ReGraw(i,j) = ReGraw(i,j) * g
        ReGraw(i,j) = ReGraw(i,j) + corr * contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGraw(i,j) = sinD * ( fac*corr*contx(i) + ReWbs )
        ImGraw(i,j) = ImGraw(i,j) + cosD * ImWbs
        ImGraw(i,j) = ImGraw(i,j) * g
        ImGraw(i,j) = ImGraw(i,j) + ImW0s + ImW3s
     end do
  end do
  return
end subroutine rawS

