subroutine lag_freq(nex, earx, Ea1, Ea2, Eb1, Eb2, contx, &
     ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3,ImW3,&
     absorbx, g, gslope,  DelAB, ABslope, nfx, fix, boost, z, gso, lens,&
     Gamma, ionvar, ReGraw, ImGraw)

  implicit none
  integer, intent(in) :: nex, nfx, ionvar
  integer, intent(in) :: Ea1, Ea2, Eb1, Eb2
  real   , intent(in) :: g, gslope, DelAB, ABslope, boost, z, gso, lens, Gamma
  real   , intent(in) :: earx(0:nex), contx(nex), absorbx(nex)
  real   , intent(in) :: fix(0:nfx)
  real   , intent(in) :: ReW0(nex,nfx), ImW0(nex,nfx), &
       ReW1(nex,nfx), ImW1(nex,nfx), ReW2(nex,nfx), ImW2(nex,nfx), &
       ReW3(nex,nfx), ImW3(nex,nfx)
       
                         
  real,    intent(out):: ReGraw(nfx), ImGraw(nfx)
  real                :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
  real                :: corr, sinD, cosD, E, fac, ReW0s, ImW0s,&
       ReWbs, ImWbs, ReW3s, ImW3s, gsoz
  real                :: f, DelAB_nu, g_nu
  integer             :: i, j



  gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
  corr = lens * gsoz**Gamma  !Correction factor for direct component
  !Now calculate the cross-spectrum (/complex covariance)
  do j = 1, nfx
     f = (fix(j) + fix(j - 1) ) * 0.5
     ! f = floHz * (fhiHz/floHz)**(  (real(j)-0.5) / real(nf) )
     ! DelAB_nu = DelAB * (floHz / f)
     DelAB_nu = DelAB * ( (fix(1) + fix(0)) * 0.5  / f)**ABslope
     g_nu     = g * ( (fix(1) + fix(0)) * 0.5  / f)**gslope
     sinD = sin(DelAB_nu)
     cosD = cos(DelAB_nu)
     ReGrawEa = 0.0
     ImGrawEa = 0.0
     ReGrawEb = 0.0
     ImGrawEb = 0.0     
     do i = Ea1, Ea2
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
        ReGrawEa = ReGrawEa + cosD * ( fac * corr * contx(i) + ReWbs )
        ReGrawEa = ReGrawEa - sinD * ImWbs
        ReGrawEa = ReGrawEa * g_nu
        ReGrawEa = ReGrawEa + corr * contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGrawEa = ImGrawEa + sinD * ( fac * corr * contx(i) + ReWbs )
        ImGrawEa = ImGrawEa + cosD * ImWbs
        ImGrawEa = ImGrawEa * g_nu
        ImGrawEa = ImGrawEa + ImW0s + ImW3s
        !Account for absorption
        ReGrawEa = ReGrawEa * absorbx(i)
        ImGrawEa = ImGrawEa * absorbx(i)
     end do
     
     do i = Eb1, Eb2
        fac = log(gsoz/E)
        !Multiply by boost parameter and group like terms
        ReW0s = boost * ReW0(i,j)
        ImW0s = boost * ImW0(i,j)
        ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
        ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
        ReW3s = ionvar * boost * ReW3(i,j)
        ImW3s = ionvar * boost * ImW3(i,j)
        !Real part
        ReGrawEb = ReGrawEb + cosD * ( fac * corr * contx(i) + ReWbs )
        ReGrawEb = ReGrawEb - sinD * ImWbs
        ReGrawEb = ReGrawEb * g_nu
        ReGrawEb = ReGrawEb + corr * contx(i) + ReW0s + ReW3s
        !Imaginary part
        ImGrawEb = ImGrawEb + sinD * ( fac * corr * contx(i) + ReWbs )
        ImGrawEb = ImGrawEb + cosD * ImWbs
        ImGrawEb = ImGrawEb * g_nu
        ImGrawEb = ImGrawEb + ImW0s + ImW3s
        !Account for absorption
        ReGrawEb = ReGrawEb * absorbx(i)
        ImGrawEb = ImGrawEb * absorbx(i)
     end do

     ! write(21, *) f, ReGrawEa
     ! write(22, *) f, ReGrawEb
     ! write(23, *) f, ImGrawEa
     ! write(24, *) f, ImGrawEb

!Now cross-spectrum between the two energy bands
     ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
     ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
     
  end do

end subroutine lag_freq
