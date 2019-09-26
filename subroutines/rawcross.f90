subroutine rawcross(nex,earx,nf,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,Gamma,ionvar,ReGraw,ImGraw)
! Calculates the pre-normalisation cross-spectrum (where the norm is alpha*e^{i\phi_A}).
!
! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
! Output: un-normalized cross-spectrum Graw
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
! 
! Outputs
! ReGraw(1:nex,1:nf)    Real part of G(E,nu)      - in specific photon flux (photar/dE)
! ImGraw(1:nex,1:nf)    Imaginary part of G(E,nu) - in specific photon flux (photar/dE)
  implicit none
  integer nex,nf,ionvar
  real earx(0:nex),corr,contx(nex),ReW0(nex,nf),ImW0(nex,nf)
  real ReW1(nex,nf),ImW1(nex,nf),ReW2(nex,nf),ImW2(nex,nf),ReW3(nex,nf),ImW3(nex,nf)
  real g,DelAB,boost,z,gso,lens,Gamma,ReGraw(nex,nf),ImGraw(nex,nf)
  real sinD,cosD,E,fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s,gsoz
  real fac2,lnS,num,den
  integer i,j
  sinD = sin(DelAB)
  cosD = cos(DelAB)
  gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
  corr = lens * gsoz**Gamma  !Correction factor for direct component
  !Calculate ln(gsoz/S), the pre-factor for W3(E,nu)
  num = 0.0
  den = 0.0
  do i = 1,nex
     E   = 0.5 * ( earx(i) + earx(i-1) )
     num = num + E * log(E) * contx(i)
     den = den + E * contx(i)
  end do
  lnS = num / den
  fac2 = log(gsoz) - lnS
  !Now calculate the cross-spectrum (/complex covariance)
  do j = 1,nf
     do i = 1,nex
        E   = 0.5 * ( earx(i) + earx(i-1) )
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        ReGraw(i,j) = cosD * ( fac*corr*contx(i) + ReWbs + fac2*ReW3s )
        ReGraw(i,j) = ReGraw(i,j) - sinD * ( ImWbs + fac2*ImW3s )
        ReGraw(i,j) = ReGraw(i,j) * g
        ReGraw(i,j) = ReGraw(i,j) + corr*contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGraw(i,j) = sinD * ( fac*corr*contx(i) + ReWbs + fac2*ReW3s  )
        ImGraw(i,j) = ImGraw(i,j) + cosD * ( ImWbs + fac2*ImW3s )
        ImGraw(i,j) = ImGraw(i,j) * g
        ImGraw(i,j) = ImGraw(i,j) + ImW0s + ImW3s
     end do
  end do
  return
end subroutine rawcross

