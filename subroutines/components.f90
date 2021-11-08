subroutine write_components(ne,ear,nex,earx,nf,nlp,contx,absorbx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,floHz,fhiHz,&
                            ReIm,DelA,DelAB,g,boost,z,gso,lens,Gamma,ionvar,resp_matr)                      
    !this subroutine separates the components from the model, calculates each cross spectrum including the effects of absorption,
    !folds the response matrix if desired, calls the phase correction, averages over frequnecy, and prints each different components
    !to a new file. This code repeats a lot and it's a bit of a monstrosity, mostly because it's annoying to separate the transfer
    !functions W0, W1 etc into model components easily. Apologies.
    
    implicit none
    integer :: ne,nex,nf,nlp,ionvar,ReIm,resp_matr
    real :: ear(0:ne),earx(0:nex),corr,contx(nex,nlp),absorbx(nex),ReW0(nex,nf),ImW0(nex,nf)
    real :: ReW1(nex,nf),ImW1(nex,nf),ReW2(nex,nf),ImW2(nex,nf),ReW3(nex,nf),ImW3(nex,nf)
    real :: g,DelA,DelAB,boost,z,gso(nlp),lens(nlp),Gamma
    real :: fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s
    real :: tempRe,tempIm,dE
    real :: f,floHz,fhiHz
    double precision :: fc
    double precision, parameter :: pi = acos(-1.d0)
    integer :: i,j
    !individual components transfer functions (S) and cross spectrum (G) dynamic allocation
    real   , dimension(:,:)    , allocatable :: ReSPL,ImSPL,ReSLT,ImSLT,ReSPR,ImSPR
    real   , dimension(:,:)    , allocatable :: ReGPLcorr,ImGPLcorr,ReGLTcorr,ImGLTcorr,ReGPRcorr,ImGPRcorr
    real   , dimension(:,:)    , allocatable :: ReGPL,ImGPL,ReGLT,ImGLT,ReGPR,ImGPR
    real :: ReGPLbar(nex),ImGPLbar(nex),ReGLTbar(nex),ImGLTbar(nex),ReGPRbar(nex),ImGPRbar(nex)
    !Arrays for each component that make up the final output to file
    real :: ener(ne),ReSPLprint(ne),ImSPLprint(ne),ReSLTprint(ne),ImSLTprint(ne),ReSPRprint(ne),ImSPRprint(ne)

    !Allocate model component matrixes
    allocate( ReSPL(nex,nf) )
    allocate( ImSPL(nex,nf) )
    allocate( ReSLT(nex,nf) )
    allocate( ImSLT(nex,nf) )
    allocate( ReSPR(nex,nf) )
    allocate( ImSPR(nex,nf) )    
    !Allocate cross spectra
    allocate( ReGPL(nex,nf) )
    allocate( ImGPL(nex,nf) )
    allocate( ReGLT(nex,nf) )
    allocate( ImGLT(nex,nf) )
    allocate( ReGPR(nex,nf) )
    allocate( ImGPR(nex,nf) )
    allocate( ReGPLcorr(nex,nf) )
    allocate( ImGPLcorr(nex,nf) )
    allocate( ReGLTcorr(nex,nf) )
    allocate( ImGLTcorr(nex,nf) )
    allocate( ReGPRcorr(nex,nf) )
    allocate( ImGPRcorr(nex,nf) )
    
    !This stores each component contribution in the Re/Im matrices: PL includeds the continuum contributions, LT the light travel
    !time, PR the pivoting reflection, IV the ionization variations    
    call components(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,&
               Gamma,ionvar,ReSPL,ImSPL,ReSLT,ImSLT,ReSPR,ImSPR)    
    !Include the effects of absorption in each model component matrix              
    do j = 1, nf
        do i = 1, nex
            ReSPL(i,j) = ReSPL(i,j) * absorbx(i)
            ImSPL(i,j) = ImSPL(i,j) * absorbx(i)
            ReSLT(i,j) = ReSLT(i,j) * absorbx(i)
            ImSLT(i,j) = ImSLT(i,j) * absorbx(i)
            ReSPR(i,j) = ReSPR(i,j) * absorbx(i)
            ImSPR(i,j) = ImSPR(i,j) * absorbx(i)
        end do
    end do    
    !Calculate raw cross-spectrum from S(E,\nu) and the reference band parameters, for each component separately
    if (ReIm .gt. 0.0) then
        call propercross(nex,nf,earx,ReSPL,ImSPL,ReGPL,ImGPL,resp_matr)
        call propercross(nex,nf,earx,ReSLT,ImSLT,ReGLT,ImGLT,resp_matr)
        call propercross(nex,nf,earx,ReSPR,ImSPR,ReGPR,ImGPR,resp_matr)
    else
        call propercross_NOmatrix(nex,nf,earx,ReSPL,ImSPL,ReGPL,ImGPL)
        call propercross_NOmatrix(nex,nf,earx,ReSLT,ImSLT,ReGLT,ImGLT)
        call propercross_NOmatrix(nex,nf,earx,ReSPR,ImSPR,ReGPR,ImGPR)
    endif              
    !Apply phase correction parameter to the cross-spectral model (for bad calibration)
    do j = 1,nf
        do i = 1,nex
            ReGPLcorr(i,j) = cos(DelA) * ReGPL(i,j) - sin(DelA) * ImGPL(i,j)
            ImGPLcorr(i,j) = cos(DelA) * ImGPL(i,j) + sin(DelA) * ReGPL(i,j)
            ReGLTcorr(i,j) = cos(DelA) * ReGLT(i,j) - sin(DelA) * ImGLT(i,j)
            ImGLTcorr(i,j) = cos(DelA) * ImGLT(i,j) + sin(DelA) * ReGLT(i,j)
            ReGPRcorr(i,j) = cos(DelA) * ReGPR(i,j) - sin(DelA) * ImGPR(i,j)
            ImGPRcorr(i,j) = cos(DelA) * ImGPR(i,j) + sin(DelA) * ReGPR(i,j)
        end do
    end do
    !Calculate frequency-averaged spectra 
    ReGPLbar = 0.0
    ImGPLbar = 0.0
    ReGLTbar = 0.0
    ImGLTbar = 0.0
    ReGPRbar = 0.0
    ImGPRbar = 0.0
    fc = 0.5d0 * ( floHz + fhiHz )   
    fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
    do j = 1,nf
        f = floHz * (fhiHz/floHz)**( (real(j)-0.5) / real(nf) )
        do i = 1,nex
            ReGPLbar(i) = ReGPLbar(i) + ReGPLcorr(i,j) / f
            ImGPLbar(i) = ImGPLbar(i) + ImGPLcorr(i,j) / f
            ReGLTbar(i) = ReGLTbar(i) + ReGLTcorr(i,j) / f
            ImGLTbar(i) = ImGLTbar(i) + ImGLTcorr(i,j) / f
            ReGPRbar(i) = ReGPRbar(i) + ReGPRcorr(i,j) / f
            ImGPRbar(i) = ImGPRbar(i) + ImGPRcorr(i,j) / f
        end do
    end do
    ReGPLbar = ReGPLbar * fac
    ImGPLbar = ImGPLbar * fac
    ReGLTbar = ReGLTbar * fac
    ImGLTbar = ImGLTbar * fac
    ReGPRbar = ReGPRbar * fac
    ImGPRbar = ImGPRbar * fac       
    !calculate energy array for output; for the energy we take average energy of the bin
    do i=1,ne 
        ener(i) = (ear(i)+ear(i-1))/2.   
    end do    
    !Rebin and write output depending on ReIm parameter
    !note that xspec gets output in e.g. lags*dE, and we want just the lags, so a factor dE needs to be included
    open (unit = 10, file = 'Output/PivotingPL.dat', status='replace', action = 'write')
    open (unit = 11, file = 'Output/LightTravelTime.dat', status='replace', action = 'write')
    open (unit = 12, file = 'Output/PivotingReflection.dat', status='replace', action = 'write')  
    if (abs(ReIm) .le. 4) then
        call crebin(nex,earx,ReGPLbar,ImGPLbar,ne,ear,ReSPLprint,ImSPLprint)        
        call crebin(nex,earx,ReGLTbar,ImGLTbar,ne,ear,ReSLTprint,ImSLTprint)
        call crebin(nex,earx,ReGPRbar,ImGPRbar,ne,ear,ReSPRprint,ImSPRprint)
        if (abs(ReIm) .eq. 1 ) then         !Real part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), ReSPLprint(i)/dE
                write (11,*) ener(i), ReSLTprint(i)/dE
                write (12,*) ener(i), ReSPRprint(i)/dE
            end do    
        else if (abs(ReIm) .eq. 2) then     !Imaginary part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), ImSPLprint(i)/dE
                write (11,*) ener(i), ImSLTprint(i)/dE
                write (12,*) ener(i), ImSPRprint(i)/dE
            end do
        else if (abs(ReIm) .eq. 3) then     !Modulus
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), sqrt(ReSPLprint(i)**2 + ImSPLprint(i)**2)/dE
                write (11,*) ener(i), sqrt(ReSLTprint(i)**2 + ImSLTprint(i)**2)/dE
                write (12,*) ener(i), sqrt(ReSPRprint(i)**2 + ImSPRprint(i)**2)/dE
            end do
        else if (abs(ReIm) .eq. 4) then     !Time lag (s)
            do i = 1,ne
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), atan2(ImSPLprint(i),ReSPLprint(i)) / ( 2.0*pi*fc )
                write (11,*) ener(i), atan2(ImSLTprint(i),ReSLTprint(i)) / ( 2.0*pi*fc ) 
                write (12,*) ener(i), atan2(ImSPRprint(i),ReSPRprint(i)) / ( 2.0*pi*fc )
            end do
        end if
    else
        call cfoldandbin(nex,earx,ReGPLbar,ImGPLbar,ne,ear,ReSPLprint,ImSPLprint,resp_matr)
        call cfoldandbin(nex,earx,ReGLTbar,ImGLTbar,ne,ear,ReSLTprint,ImSLTprint,resp_matr)
        call cfoldandbin(nex,earx,ReGPRbar,ImGPRbar,ne,ear,ReSPRprint,ImSPRprint,resp_matr)
        if (abs(ReIm) .eq. 5) then          !Modulus
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), sqrt(ReSPLprint(i)**2 + ImSPLprint(i)**2)/dE
                write (11,*) ener(i), sqrt(ReSLTprint(i)**2 + ImSLTprint(i)**2)/dE
                write (12,*) ener(i), sqrt(ReSPRprint(i)**2 + ImSPRprint(i)**2)/dE
            end do
        else if (abs(ReIm) .eq. 6) then     !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                write (10,*) ener(i), atan2(ImSPLprint(i),ReSPLprint(i)) / ( 2.0*pi*fc ) 
                write (11,*) ener(i), atan2(ImSLTprint(i),ReSLTprint(i)) / ( 2.0*pi*fc )
                write (12,*) ener(i), atan2(ImSPRprint(i),ReSPRprint(i)) / ( 2.0*pi*fc ) 
            end do
        end if
    end if
    close(10)
    close(11)
    close(12)
end subroutine write_ComponentS

subroutine components(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,&
     Gamma,ionvar,ReSPL,ImSPL,ReSLT,ImSLT,ReSPR,ImSPR)
    ! Calculates the FT of the spectrum components before multiplying by the absorption model
    ! This is essentially the same as S, but it returns the FT for each component that contributes to the lags:
    ! FTcont for the pivoting continuum - stored in ReSPL and ImPLt
    ! W0 for the light travel lags - stored in ReSLT and ImSLT
    ! W1 and W2 for the pivoting reflection - stored in ReSPR and ImSPR
    ! W3 for the ionization variations - stored in ReSIV and ImSIV
    ! 
    ! Input: continuum model (contx is f(E) in the papers) and the reflection transfer functions
    ! Output: S(E,\nu) before multiplying by the absorption model for each separate component
    ! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
    ! 
    ! Inputs:
    ! nex,earx(0:nex)       Internal energy grid
    ! nf                    Number of frequency bins
    ! contx(1:nex)          Time-independent continuum model
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
      
    ! Outputs
    ! ReSPL(1:nex,1:nf)    Real part of S(E,nu) for the continuum only, in specific photon flux (photar/dE)
    ! ImSPL(1:nex,1:nf)    Imaginary part of S(E,nu) for the continuum only,  in specific photon flux (photar/dE)
    ! ReSLT(1:nex,1:nf)    Real part of S(E,nu) for the light travel time only, in specific photon flux (photar/dE)
    ! ImSLT(1:nex,1:nf)    Imaginary part of S(E,nu) for the light travel time only,  in specific photon flux (photar/dE)
    ! ReSPR(1:nex,1:nf)    Real part of S(E,nu) for the pivoting reflection only, in specific photon flux (photar/dE)
    ! ImSPR(1:nex,1:nf)    Imaginary part of S(E,nu) for the pivoting reflection only,  in specific photon flux (photar/dE)

    implicit none
    integer nex,nf,ionvar,nlp
    real earx(0:nex),corr,contx(nex,nlp),cont_1(nex),cont_2(nex),ReW0(nex,nf),ImW0(nex,nf)
    real ReW1(nex,nf),ImW1(nex,nf),ReW2(nex,nf),ImW2(nex,nf),ReW3(nex,nf),ImW3(nex,nf)
    real g,DelAB,boost,z,gso(nlp),lens(nlp),Gamma
    real ReSPL(nex,nf),ImSPL(nex,nf),ReSLT(nex,nf),ImSLT(nex,nf)
    real ReSPR(nex,nf),ImSPR(nex,nf)
    real sinD,cosD,E,fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s,gsoz
    integer i,j,m
    sinD = sin(DelAB)
    cosD = cos(DelAB)
    cont_1 = 0.
    cont_2 = 0.
    !gsoz = gso / (1.0+z)       !blueshift corrected for expansion of the Universe
    !corr = 1.!lens * gsoz**Gamma  !Correction factor for direct component
    
    !Now calculate the cross-spectrum (/complex covariance)
    !Redo with Gullo's convention: separate pivoting PL, pivoting reflection+light crossing time, light crossing time only.
    !this part is less fucked than it used to be
    do j = 1,nf
        do i = 1,nex
            E   = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                fac = log(gso(m)/((1.0+z)*E)) !TBD ask about this factor!!!!
                cont_1(i) = cont_1(i) + contx(i,m)
                cont_2(i) = cont_2(i) + fac*contx(i,m)   
            end do 
            !Multiply by boost parameter and group like terms
            ReW0s = boost * ReW0(i,j)
            ImW0s = boost * ImW0(i,j)
            ReWbs = boost * ( ReW1(i,j) + ReW2(i,j) )
            ImWbs = boost * ( ImW1(i,j) + ImW2(i,j) )
            ReW3s = ionvar * boost * ReW3(i,j)
            ImW3s = ionvar * boost * ImW3(i,j)
            !Real part of all the components
            ReSPL(i,j) = g * cosD * cont_2(i) + cont_1(i)
            ReSLT(i,j) = ReW0s + cont_1(i) 
            ReSPR(i,j) = g * (cosD * ReWbs - sinD * ImWbs) + ReW3s  + ReW0s + cont_1(i)             
            !Imaginary part of all the components
            ImSPL(i,j) = g * sinD * cont_2(i)
            ImSLT(i,j) = ImW0s 
            ImSPR(i,j) = g * (sinD * ReWbs + cosD * ImWbs) + ImW3s + ImW0s
        end do
    end do
    return
end subroutine










