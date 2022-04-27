subroutine lag_freq(nlp,nex,earx,contx,absorbx,Emin,Emax,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                    g,DelAB,nf,fix,boost,z,gso,lens,Gamma,ionvar,ReGraw,ImGraw)

! Calculates the energy-averaged cross spectrum in the energy bins of interest. REDO ALL THIS COMMENT IT'S WRONG
!
! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
! Output: Sraw(E,\nu) before multiplying by the absorption model
! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
!
! Inputs:
! nlp                   Number of lamp-posts 
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
! lens                  Lensing factor -- note: not necessary if we do the correction directly when calling the continuum
! Gamma                 Photon index
! ionvar                Integer: ionvar=1 means include ionization variations, ionvar=0 means don't.
! DC                    Integer: DC=1 means this is the DC component, DC=0 means opposite.
  
! Outputs
! ReGraw(1:nex,1:nf)    Real part of Sraw(E,nu)      - in specific photon flux (photar/dE)
! ImGraw(1:nex,1:nf)    Imaginary part of Sraw(E,nu) - in specific photon flux (photar/dE)

    implicit none
    integer, intent(in) :: nex,nf,ionvar,nlp
    integer             :: Ea1,Ea2,Eb1,Eb2       
    real   , intent(in) :: g(nlp),DelAB(nlp),boost,z,Gamma,Emin,Emax
    real   , intent(in) :: gso(nlp),lens(nlp)
    real                :: gslope,ABslope
    real   , intent(in) :: earx(0:nex),contx(nex,nlp),absorbx(nex),fix(0:nf)
    real   , intent(in) :: ReW0(nex,nf),ImW0(nex,nf),ReW3(nex,nf),ImW3(nex,nf),&
                           ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf)                       
    real,    intent(out):: ReGraw(nf),ImGraw(nf)
    real                :: ReGrawEa,ImGrawEa,ReGrawEb,ImGrawEb
    real                :: sinD,cosD,E,fac,TempReG,TempImG
    real                :: ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s 
    real                :: f,DelAB_nu,g_nu
    integer             :: i,j,m

    call energy_bounds(nex,Emin,Emax,Ea1,Ea2,Eb1,Eb2)

    gslope = 1.
    ABslope = 1.
  
    !Now calculate the cross-spectrum (/complex covariance), including absorption
    do j = 1, nf
        f = (fix(j) + fix(j-1)) * 0.5
        DelAB_nu = DelAB(nlp) * ((fix(1) + fix(0))*0.5/f)**ABslope
        g_nu     = g(nlp) * ((fix(1) + fix(0))*0.5/f)**gslope
        sinD = sin(DelAB_nu)
        cosD = cos(DelAB_nu)
        ReGrawEa = 0.0
        ImGrawEa = 0.0
        ReGrawEb = 0.0
        ImGrawEb = 0.0     
        do i = Ea1, Ea2
            !Multiply by boost parameter and group like terms
            ReW0s = boost * ReW0(i,j)
            ImW0s = boost * ImW0(i,j)                
            ReW3s = ionvar * boost * ReW3(i,j)
            ImW3s = ionvar * boost * ImW3(i,j)  
            E   = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                fac = log(gso(m)/((1.0+z)*E))
                sinD = sin(DelAB_nu)
                cosD = cos(DelAB_nu) 
                ReWbs = boost * (ReW1(m,i,j) + ReW2(m,i,j))
                ImWbs = boost * (ImW1(m,i,j) + ImW2(m,i,j))
                !Real part contribution for each LP
                TempReG = g(m) * (cosD*(fac*contx(i,m) + ReWbs) - sinD*ImWbs) 
                TempReG = TempReG + contx(i,m)         
                !Imaginary part contribution for each LP 
                TempImG = g(m) * (sinD*(fac*contx(i,m) + ReWbs) + cosD*ImWbs) 
            end do 
            !Real part shared by both
            TempReG = TempReG + ReW0s + ReW3s
            ReGrawEa = ReGrawEa + TempReG
            !Imaginary part shared by both
            TempImG = TempImG + ImW0s + ImW3s
            ImGrawEa = ImGrawEa + TempImG 
            ReGrawEa = ReGrawEa * absorbx(i)
            ImGrawEa = ImGrawEa * absorbx(i)
        end do

        do i = Eb1, Eb2
            !Multiply by boost parameter and group like terms
            ReW0s = boost * ReW0(i,j)
            ImW0s = boost * ImW0(i,j)                
            ReW3s = ionvar * boost * ReW3(i,j)
            ImW3s = ionvar * boost * ImW3(i,j)  
            E   = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                fac = log(gso(m)/((1.0+z)*E))
                sinD = sin(DelAB_nu)
                cosD = cos(DelAB_nu) 
                ReWbs = boost * (ReW1(m,i,j) + ReW2(m,i,j))
                ImWbs = boost * (ImW1(m,i,j) + ImW2(m,i,j))
                !Real part contribution for each LP
                TempReG = g(m) * (cosD*(fac*contx(i,m) + ReWbs) - sinD*ImWbs) 
                TempReG = TempReG + contx(i,m)         
                !Imaginary part contribution for each LP 
                TempImG = g(m) * (sinD*(fac*contx(i,m) + ReWbs) + cosD*ImWbs) 
            end do 
            !Real part shared by both
            TempReG = TempReG + ReW0s + ReW3s
            ReGrawEb = ReGrawEb + TempReG
            !Imaginary part shared by both
            TempImG = TempImG + ImW0s + ImW3s
            ImGrawEb = ImGrawEb + TempImG 
            ReGrawEb = ReGrawEb * absorbx(i)
            ImGrawEb = ImGrawEb * absorbx(i)
        end do
        !Now cross-spectrum between the two energy bands
        !note: here the conjugate is b 
        ReGraw(j) = (ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb)
        ImGraw(j) = (ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb)
    end do

end subroutine lag_freq

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

end subroutine energy_bounds
