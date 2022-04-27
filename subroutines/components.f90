subroutine write_components(ne,ear,nex,earx,nf,nlp,contx,absorbx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,floHz,fhiHz,&
                            ReIm,DelA,DelAB,g,boost,z,gso,lens,Gamma,ionvar,resp_matr)                      
    !this subroutine separates the components from the model, calculates each cross spectrum including the effects of absorRTion,
    !folds the response matrix if desired, calls the phase correction, averages over frequnecy, and prints each different components
    !to a new file. This code repeats a lot and it's a bit of a monstrosity, mostly because it's annoying to separate the transfer
    !functions W0, W1 etc into model components easily. Apologies.
    !Nomenclature for each contribution: 
    !PL - pivoting continuum 
    !LT - light travel time only
    !PR - pivoting reflection of each source
    !RT - total reflection lag due to light travel time, pivoting of each reflection signal, and ionization variations  
    
    implicit none
    integer :: ne,nex,nf,nlp,ionvar,ReIm,resp_matr
    real :: ear(0:ne),earx(0:nex),corr,contx(nex,nlp),absorbx(nex),ReW0(nex,nf),ImW0(nex,nf)
    real :: ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nex,nf),ImW3(nex,nf)
    real :: g(nlp),DelA,DelAB(nlp),boost,z,gso(nlp),lens(nlp),Gamma
    real :: fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s
    real :: tempRe,tempIm,dE
    real :: f,floHz,fhiHz
    double precision :: fc
    double precision, parameter :: pi = acos(-1.d0)
    integer :: i,j,m
    !indiRTdual components transfer functions (S) and cross spectrum (G) dynamic allocation
    real, dimension(:,:,:), allocatable :: ReSPL,ImSPL,ReSPR,ImSPR
    real, dimension(:,:)  , allocatable :: ReSLT,ImSLT,ReSRT,ImSRT    
    real, dimension(:,:,:), allocatable :: ReGPLcorr,ImGPLcorr,ReGPRcorr,ImGPRcorr
    real, dimension(:,:)  , allocatable :: ReGLTcorr,ImGLTcorr,ReGRTcorr,ImGRTcorr
    real, dimension(:,:)  , allocatable :: ReGLT,ImGLT,ReGRT,ImGRT 
    real, dimension(:,:,:), allocatable :: ReGPL,ImGPL,ReGPR,ImGPR

    real :: ReGPLbar(nlp,nex),ImGPLbar(nlp,nex),ReGPRbar(nlp,nex),ImGPRbar(nlp,nex)
    real :: ReGLTbar(nex),ImGLTbar(nex),ReGRTbar(nex),ImGRTbar(nex)
    !Arrays for each component that make up the final output to file    
    real :: ReSPLprint(nlp,ne),ImSPLprint(nlp,ne),ReSPRprint(nlp,ne),ImSPRprint(nlp,ne)
    real :: ener(ne),ReSLTprint(ne),ImSLTprint(ne),ReSRTprint(ne),ImSRTprint(ne)
    !strings to open the output files
    character (len=17) path_name
    character (len=2)  path_lp
    character (len=4)  path_ext 
    character (len=24) path_full
    
    !Allocate model component matrixes
    allocate( ReSPL(nlp,nex,nf) )
    allocate( ImSPL(nlp,nex,nf) )
    allocate( ReSPR(nlp,nex,nf) )
    allocate( ImSPR(nlp,nex,nf) )  
    allocate( ReSLT(nex,nf) )
    allocate( ImSLT(nex,nf) ) 
    allocate( ReSRT(nex,nf) )
    allocate( ImSRT(nex,nf) ) 
    !Allocate cross spectra
    allocate( ReGPL(nlp,nex,nf) )
    allocate( ImGPL(nlp,nex,nf) )
    allocate( ReGPR(nlp,nex,nf) )
    allocate( ImGPR(nlp,nex,nf) )
    allocate( ReGLT(nex,nf) )
    allocate( ImGLT(nex,nf) )
    allocate( ReGRT(nex,nf) )
    allocate( ImGRT(nex,nf) )

    allocate( ReGPLcorr(nlp,nex,nf) )
    allocate( ImGPLcorr(nlp,nex,nf) )
    allocate( ReGPRcorr(nlp,nex,nf) )
    allocate( ImGPRcorr(nlp,nex,nf) )
    allocate( ReGLTcorr(nex,nf) )
    allocate( ImGLTcorr(nex,nf) )
    allocate( ReGRTcorr(nex,nf) )
    allocate( ImGRTcorr(nex,nf) )
    
    !This stores each component contribution in the Re/Im matrices: PL includeds the continuum contributions, LT the light travel
    !time, PR the pivoting reflection, IV the ionization variations    
    call components(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,&
                    Gamma,ionvar,ReSPL,ImSPL,ReSLT,ImSLT,ReSPR,ImSPR,ReSRT,ImSRT)    
    !Include the effects of absorRTion in each model component matrix              
    do j = 1, nf
        do i = 1, nex
            do m=1,nlp 
                ReSPL(m,i,j) = ReSPL(m,i,j) * absorbx(i)
                ImSPL(m,i,j) = ImSPL(m,i,j) * absorbx(i)
                ReSPR(m,i,j) = ReSPR(m,i,j) * absorbx(i)
                ImSPR(m,i,j) = ImSPR(m,i,j) * absorbx(i)
            end do            
            ReSLT(i,j) = ReSLT(i,j) * absorbx(i)
            ImSLT(i,j) = ImSLT(i,j) * absorbx(i)
            ReSRT(i,j) = ReSRT(i,j) * absorbx(i)
            ImSRT(i,j) = ImSRT(i,j) * absorbx(i)
        end do
    end do    
    !Calculate raw cross-spectrum from S(E,\nu) and the reference band parameters, for each component separately
    if (ReIm .gt. 0.0) then
        do m=1,nlp 
            call propercross(nex,nf,earx,ReSPL(m,:,:),ImSPL(m,:,:),ReGPL(m,:,:),ImGPL(m,:,:),resp_matr)
            call propercross(nex,nf,earx,ReSPR(m,:,:),ImSPR(m,:,:),ReGPR(m,:,:),ImGPR(m,:,:),resp_matr)
        end do
        call propercross(nex,nf,earx,ReSLT,ImSLT,ReGLT,ImGLT,resp_matr)        
        call propercross(nex,nf,earx,ReSRT,ImSRT,ReGRT,ImGRT,resp_matr)
    else
        do m=1,nlp 
            call propercross_NOmatrix(nex,nf,earx,ReSPL(m,:,:),ImSPL(m,:,:),ReGPL(m,:,:),ImGPL(m,:,:))
            call propercross_NOmatrix(nex,nf,earx,ReSPR(m,:,:),ImSPR(m,:,:),ReGPR(m,:,:),ImGPR(m,:,:))
        end do
        call propercross_NOmatrix(nex,nf,earx,ReSLT,ImSLT,ReGLT,ImGLT)
        call propercross_NOmatrix(nex,nf,earx,ReSRT,ImSRT,ReGRT,ImGRT)
    endif              
    !Apply phase correction parameter to the cross-spectral model (for bad calibration)
    do j = 1,nf
        do i = 1,nex
            do m=1,nlp 
                ReGPLcorr(m,i,j) = cos(DelA) * ReGPL(m,i,j) - sin(DelA) * ImGPL(m,i,j)
                ImGPLcorr(m,i,j) = cos(DelA) * ImGPL(m,i,j) + sin(DelA) * ReGPL(m,i,j)
                ReGPRcorr(m,i,j) = cos(DelA) * ReGPR(m,i,j) - sin(DelA) * ImGPR(m,i,j)
                ImGPRcorr(m,i,j) = cos(DelA) * ImGPR(m,i,j) + sin(DelA) * ReGPR(m,i,j)
            end do
            ReGLTcorr(i,j) = cos(DelA) * ReGLT(i,j) - sin(DelA) * ImGLT(i,j)
            ImGLTcorr(i,j) = cos(DelA) * ImGLT(i,j) + sin(DelA) * ReGLT(i,j)
            ReGRTcorr(i,j) = cos(DelA) * ReGRT(i,j) - sin(DelA) * ImGRT(i,j)
            ImGRTcorr(i,j) = cos(DelA) * ImGRT(i,j) + sin(DelA) * ReGRT(i,j)
        end do
    end do
    !Calculate frequency-averaged spectra 
    ReGPLbar = 0.0
    ImGPLbar = 0.0
    ReGPRbar = 0.0
    ImGPRbar = 0.0
    ReGLTbar = 0.0
    ImGLTbar = 0.0
    ReGRTbar = 0.0
    ImGRTbar = 0.0
    fc = 0.5d0 * ( floHz + fhiHz )   
    fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
    do j = 1,nf
        f = floHz * (fhiHz/floHz)**( (real(j)-0.5) / real(nf) )
        do i = 1,nex
            do m=1,nlp 
                ReGPLbar(m,i) = ReGPLbar(m,i) + ReGPLcorr(m,i,j) / f
                ImGPLbar(m,i) = ImGPLbar(m,i) + ImGPLcorr(m,i,j) / f
                ReGPRbar(m,i) = ReGPRbar(m,i) + ReGPRcorr(m,i,j) / f
                ImGPRbar(m,i) = ImGPRbar(m,i) + ImGPRcorr(m,i,j) / f
            end do
            ReGLTbar(i) = ReGLTbar(i) + ReGLTcorr(i,j) / f
            ImGLTbar(i) = ImGLTbar(i) + ImGLTcorr(i,j) / f
            ReGRTbar(i) = ReGRTbar(i) + ReGRTcorr(i,j) / f
            ImGRTbar(i) = ImGRTbar(i) + ImGRTcorr(i,j) / f
        end do
    end do
    ReGPLbar = ReGPLbar * fac
    ImGPLbar = ImGPLbar * fac
    ReGPRbar = ReGPRbar * fac
    ImGPRbar = ImGPRbar * fac 
    ReGLTbar = ReGLTbar * fac
    ImGLTbar = ImGLTbar * fac
    ReGRTbar = ReGRTbar * fac
    ImGRTbar = ImGRTbar * fac      
    !calculate energy array for output; for the energy we take average energy of the bin
    do i=1,ne 
        ener(i) = (ear(i)+ear(i-1))/2.   
    end do    
       
    path_ext = '.dat'
    path_name = 'Output/PivotingPL'
    do m=1,nlp 
        path_lp = '_'//char(48+m)
        path_full = path_name//path_lp//path_ext
        open (unit = 10+m, file = path_full, status='replace', action = 'write')
    end do
    path_name = 'Output/PivotingRE'
    do m=1,nlp 
        path_lp = '_'//char(48+m)
        path_full = path_name//path_lp//path_ext
        open (unit = 10+nlp+m, file = path_full, status='replace', action = 'write')
    end do
    path_name = 'Output/TravelTime'
    path_full = path_name//path_ext
    open (unit = 11+2*nlp, file = path_full, status='replace', action = 'write') 
    path_name = 'Output/TotReflect'
    path_full = path_name//path_ext
    open (unit = 12+2*nlp, file = path_full, status='replace', action = 'write') 

    if (abs(ReIm) .le. 4) then
        do m=1,nlp 
            call crebin(nex,earx,ReGPLbar(m,:),ImGPLbar(m,:),ne,ear,ReSPLprint(m,:),ImSPLprint(m,:))        
            call crebin(nex,earx,ReGPRbar(m,:),ImGPRbar(m,:),ne,ear,ReSPRprint(m,:),ImSPRprint(m,:))
        end do        
        call crebin(nex,earx,ReGLTbar,ImGLTbar,ne,ear,ReSLTprint,ImSLTprint)
        call crebin(nex,earx,ReGRTbar,ImGRTbar,ne,ear,ReSRTprint,ImSRTprint)
        if (abs(ReIm) .eq. 1 ) then         !Real part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                do m=1, nlp
                    write (10+m,*)      ener(i), ReSPLprint(m,i)/dE
                    write (10+nlp+m,*)  ener(i), ReSPRprint(m,i)/dE
                end do
                write (11+2*nlp,*)      ener(i), ReSLTprint(i)/dE
                write (12+2*nlp,*)      ener(i), ReSRTprint(i)/dE
            end do    
        else if (abs(ReIm) .eq. 2) then     !Imaginary part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                do m=1,nlp 
                    write (10+m,*)      ener(i), ImSPLprint(m,i)/dE
                    write (10+nlp+m,*)  ener(i), ImSPRprint(m,i)/dE
                end do
                write (11+2*nlp,*)      ener(i), ImSLTprint(i)/dE
                write (12+2*nlp,*)      ener(i), ImSRTprint(i)/dE
            end do
        else if (abs(ReIm) .eq. 3) then     !Modulus
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                do m=1,nlp
                    write (10+m,*)      ener(i), sqrt(ReSPLprint(m,i)**2 + ImSPLprint(m,i)**2)/dE
                    write (10+nlp+m,*)  ener(i), sqrt(ReSPRprint(m,i)**2 + ImSPRprint(m,i)**2)/dE
                end do
                write (11+2*nlp,*)      ener(i), sqrt(ReSLTprint(i)**2 + ImSLTprint(i)**2)/dE
                write (12+2*nlp,*)      ener(i), sqrt(ReSRTprint(i)**2 + ImSRTprint(i)**2)/dE
            end do
        else if (abs(ReIm) .eq. 4) then     !Time lag (s)
            do i = 1,ne
                dE = ear(i) - ear(i-1)
                do m=1,nlp 
                    write (10+m,*)      ener(i), atan2(ImSPLprint(m,i),ReSPLprint(m,i)) / ( 2.0*pi*fc )
                    write (10+nlp+m,*)  ener(i), atan2(ImSPRprint(m,i),ReSPRprint(m,i)) / ( 2.0*pi*fc )
                end do
                write (11+2*nlp,*)      ener(i), atan2(ImSLTprint(i),ReSLTprint(i)) / ( 2.0*pi*fc ) 
                write (12+2*nlp,*)      ener(i), atan2(ImSRTprint(i),ReSRTprint(i)) / ( 2.0*pi*fc ) 
            end do
        end if
    else
        do m=1,nlp 
            call cfoldandbin(nex,earx,ReGPLbar(m,:),ImGPLbar(m,:),ne,ear,ReSPLprint(m,:),ImSPLprint(m,:),resp_matr)
            call cfoldandbin(nex,earx,ReGPRbar(m,:),ImGPRbar(m,:),ne,ear,ReSPRprint(m,:),ImSPRprint(m,:),resp_matr) 
        end do
        call cfoldandbin(nex,earx,ReGLTbar,ImGLTbar,ne,ear,ReSLTprint,ImSLTprint,resp_matr)
        call cfoldandbin(nex,earx,ReGRTbar,ImGRTbar,ne,ear,ReSRTprint,ImSRTprint,resp_matr)
        if (abs(ReIm) .eq. 5) then          !Modulus
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                do m=1,nlp 
                    write (10+m,*)      ener(i), sqrt(ReSPLprint(m,i)**2 + ImSPLprint(m,i)**2)/dE
                    write (10+nlp+m,*)  ener(i), sqrt(ReSPRprint(m,i)**2 + ImSPRprint(m,i)**2)/dE
                end do
                write (11+2*nlp,*)      ener(i), sqrt(ReSLTprint(i)**2 + ImSLTprint(i)**2)/dE
                write (12+2*nlp,*)      ener(i), sqrt(ReSRTprint(i)**2 + ImSRTprint(i)**2)/dE
            end do
        else if (abs(ReIm) .eq. 6) then     !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                do m=1,nlp 
                    write (10+m,*)      ener(i), atan2(ImSPLprint(m,i),ReSPLprint(m,i)) / ( 2.0*pi*fc )
                    write (10+nlp+m,*)  ener(i), atan2(ImSPRprint(m,i),ReSPRprint(m,i)) / ( 2.0*pi*fc )
                end do
                write (11+2*nlp,*)      ener(i), atan2(ImSLTprint(i),ReSLTprint(i)) / ( 2.0*pi*fc ) 
                write (12+2*nlp,*)      ener(i), atan2(ImSRTprint(i),ReSRTprint(i)) / ( 2.0*pi*fc )
            end do
        end if
    end if
    do m=1,nlp 
        close(10+m)
        close(10+nlp+m)
    end do
    close(11+2*nlp)
    close(12+2*nlp)    
end subroutine write_ComponentS

subroutine components(nex,earx,nf,nlp,contx,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,g,DelAB,boost,z,gso,lens,&
     Gamma,ionvar,ReSPL,ImSPL,ReSLT,ImSLT,ReSPR,ImSPR,ReSRT,ImSRT)
    ! Calculates the FT of the spectrum components before multiplying by the absorRTion model
    ! This is essentially the same as S, but it returns the FT for each component that contributes to the lags:
    ! FTcont for the pivoting continuum - stored in ReSPL and ImPLt
    ! W0 for the light travel lags - stored in ReSLT and ImSLT
    ! W1 and W2 for the pivoting reflection - stored in ReSPR and ImSPR
    ! W3 for the ionization variations - stored in ReSRT and ImSRT
    !     
    ! Outputs
    ! ReSPL(nlp,1:nex,1:nf)     Real part of S(E,nu) for each continuum only
    ! ImSPL(nlp,1:nex,1:nf)     Imaginary part of S(E,nu) for each continuum only
    ! ReSLT(1:nex,1:nf)         Real part of S(E,nu) for the light travel time only
    ! ImSLT(1:nex,1:nf)         Imaginary part of S(E,nu) for the light travel time only
    ! ReSPR(nlp,1:nex,1:nf)     Real part of S(E,nu) for the pivoting reflection only for each continuum
    ! ImSPR(nlp,1:nex,1:nf)     Imaginary part of S(E,nu) for the pivoting reflection only for each continuum
    ! ReSRT(1:nex,1:nf)         Real part of S(E,nu) for the varying ionization only
    ! ImSRT(1:nex,1:nf)         Imaginary part of S(E,nu) for the varying ionization only

    implicit none
    integer nex,nf,ionvar,nlp
    real earx(0:nex),corr,contx(nex,nlp),ReW0(nex,nf),ImW0(nex,nf),ReW3(nex,nf),ImW3(nex,nf)
    real ReW1(nlp,nex,nf),ImW1(nlp,nex,nf),ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),contx_sum(nex)
    real g(nlp),DelAB(nlp),boost,z,gso(nlp),lens(nlp),Gamma
    real ReSPL(nlp,nex,nf),ImSPL(nlp,nex,nf),ReSPR(nlp,nex,nf),ImSPR(nlp,nex,nf)
    real ReSLT(nex,nf),ImSLT(nex,nf),ReSRT(nex,nf),ImSRT(nex,nf)
    real sinD,cosD,E,fac,ReW0s,ImW0s,ReWbs,ImWbs,ReW3s,ImW3s,gsoz
    integer i,j,m

    ReSPL = 0.
    ReSLT = 0. 
    ReSPR = 0.
    ReSRT = 0.
    ImSPL = 0.
    ImSPR = 0.
    ImSLT = 0.
    ImSRT = 0.
    
    do j = 1,nf
        contx_sum = 0.
        do i = 1,nex
        !Multiply by boost parameter and group like terms
            ReW0s = boost * ReW0(i,j)
            ImW0s = boost * ImW0(i,j)                
            ReW3s = ionvar * boost * ReW3(i,j)
            ImW3s = ionvar * boost * ImW3(i,j)  
            E   = 0.5 * ( earx(i) + earx(i-1) )
            do m=1,nlp
                fac = log(gso(m)/((1.0+z)*E))
                sinD = sin(DelAB(m))
                cosD = cos(DelAB(m)) 
                ReWbs = boost * (ReW1(m,i,j) + ReW2(m,i,j))
                ImWbs = boost * (ImW1(m,i,j) + ImW2(m,i,j))
                !pivoting lags
                ReSPL(m,i,j) = g(m)*cosD*fac*contx(i,m) + contx(i,m) !TBD FIX THIS INDEXING HERE AND EVERYWHERE ELSE IN THE CODE
                ImSPL(m,i,j) = g(m)*sinD*fac*contx(i,m)
                !individual pivoting reflection lag
                ReSPR(m,i,j) = g(m)*(cosD*ReWbs - sinD*ImWbs) + contx(i,m) 
                ImSPR(m,i,j) = g(m)*(sinD*ReWbs + cosD*ImWbs)
                !total pivoting reflection lag 
                ReSRT(i,j) = ReSRT(i,j) + ReSPR(m,i,j)
                ImSRT(i,j) = ImSRT(i,j) + ImSPR(m,i,j)
                contx_sum(i) = contx_sum(i) + contx(i,m)
            end do 
            !light travel time lags
            ReSLT(i,j) = ReW0s + contx_sum(i)
            ImSLT(i,j) = ImW0s 
            !total pivoting reflection+light travel time lags
            ReSRT(i,j) = ReSRT(i,j) + ReW0s + ReW3s
            ImSRT(i,j) = ImSRT(i,j) + ImW0s + ImW3s
        end do
    end do

    return
end subroutine
