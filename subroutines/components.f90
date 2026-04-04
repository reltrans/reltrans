module m_components
    implicit none

contains

subroutine write_components(config, model_args, arrays, ear, earx, nf, floHz, fhiHz, absorbx, tauso, gso, &
                            ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3, ImW3)             
    !this subroutine separates the components from the model, calculates each cross spectrum including the effects of absorRTion,
    !folds the response matrix if desired, calls the phase correction, averages over frequnecy, and prints each different components
    !to a new file. This code repeats a lot and it's a bit of a monstrosity, mostly because it's annoying to separate the transfer
    !functions W0, W1 etc into model components easily. Apologies.
    !Nomenclature for each contribution: 
    !PL - pivoting continuum 
    !LT - light travel time only
    !PR - pivoting reflection of each source
    !RT - total reflection lag due to light travel time, pivoting of each reflection signal, and ionization variations  

    use common_types             
    implicit none
    ! Pass structured types instead of individual variables
    type(t_config), intent(IN)     :: config
    type(t_model_arguments), intent(IN) :: model_args
    type(t_arrays), intent(IN)      :: arrays

    ! Dimensions and energy grids
    integer, intent(IN) :: nf
    real   , intent(IN) :: ear(0:), earx(0:)
    real   , intent(IN) :: absorbx(:), tauso(:), gso(:)
    real   , intent(IN) :: floHz, fhiHz

    ! Transfer function arrays (W0-W3)
    real   , intent(INOUT) :: ReW0(:,:,:), ImW0(:,:,:)
    real   , intent(INOUT) :: ReW1(:,:,:), ImW1(:,:,:)
    real   , intent(INOUT) :: ReW2(:,:,:), ImW2(:,:,:)
    real   , intent(INOUT) :: ReW3(:,:,:), ImW3(:,:,:)

    ! Internal calculation variables
    real :: fac, tempRe, tempIm, dE, corr
    real :: f, flo, fhi
    double precision :: fc
    double precision, parameter :: pi = acos(-1.d0)
    integer :: i, j, m, ne

    ! Component transfer functions and cross spectrum allocation
    real, dimension(:,:), allocatable :: ReScont, ImScont, ReSrev, ImSrev
    real, dimension(:,:), allocatable :: ReSpiv, ImSpiv, ImSion, ReSion
    real, dimension(:,:), allocatable :: ReGcont, ImGcont, ReGrev, ImGrev 
    real, dimension(:,:), allocatable :: ReGpiv, ImGpiv, ReGion, ImGion

    real, dimension(:), allocatable :: ReGcont_bar, ImGcont_bar, ReGpiv_bar, ImGpiv_bar
    real, dimension(:), allocatable :: ReGrev_bar, ImGrev_bar, ReGion_bar, ImGion_bar

    ! Final output arrays
    real, dimension(:), allocatable :: ReScont_print, ImScont_print, ReSpiv_print, ImSpiv_print
    real, dimension(:), allocatable :: ener, ReSrev_print, ImSrev_print, ReSion_print, ImSion_print
    
    character (len=30) path

    ne = size(ear) - 1
    
    ! Allocate model component matrixes
    allocate( ReScont(config%nex,nf) )
    allocate( ImScont(config%nex,nf) )
    allocate( ReSrev (config%nex,nf) )
    allocate( ImSrev (config%nex,nf) )
    allocate( ReSpiv (config%nex,nf) )
    allocate( ImSpiv (config%nex,nf) )
    allocate( ReSion (config%nex,nf) )
    allocate( ImSion (config%nex,nf) )

    ! Allocate cross spectra
    allocate( ReGcont(config%nex,nf) )
    allocate( ImGcont(config%nex,nf) )
    allocate( ReGrev (config%nex,nf) )
    allocate( ImGrev (config%nex,nf) )
    allocate( ReGpiv (config%nex,nf) )
    allocate( ImGpiv (config%nex,nf) )
    allocate( ReGion (config%nex,nf) )
    allocate( ImGion (config%nex,nf) )

    allocate(ReGcont_bar(config%nex), ImGcont_bar(config%nex), ReGpiv_bar(config%nex), ImGpiv_bar(config%nex))
    allocate(ReGrev_bar(config%nex),  ImGrev_bar(config%nex),  ReGion_bar(config%nex), ImGion_bar(config%nex))
    allocate(ReScont_print(ne), ImScont_print(ne), ReSpiv_print(ne), ImSpiv_print(ne))
    allocate(ener(ne), ReSrev_print(ne), ImSrev_print(ne), ReSion_print(ne), ImSion_print(ne))
    
    ! open (unit = 20, file = 'fort.20', status='replace', action = 'write')
    ! do m = 1, nlp
    !    do j = 1, nf
    !       do i = 1, nex
    !          write(*,*) m, j, i
    !          write(20,*) m, j, i, ReW0(m,i,j),ImW0(m,i,j),ReW1(m,i,j),ImW1(m,i,j), ReW2(m,i,j),ImW2(m,i,j),ReW3(m,i,j),ImW3(m,i,j)
    !       end do
    !    end do
    ! enddo
    ! close(20)
    
 
    write(*,*) 'inside components'
    !This stores each component contribution in the Re/Im matrices 
    ! Logic updated to check components via model_args object
    if (model_args%nlp .gt. 1 .and. model_args%beta_p .eq. 0.) then  
        call components_nocoh(config%nex, earx, nf, real(config%flo), real(config%fhi), model_args%nlp, arrays%contx, absorbx, &
                              tauso, gso, ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3, ImW3, &
                              real(model_args%h), real(model_args%zcos), real(model_args%Gamma), &
                              real(model_args%eta), model_args%boost, model_args%g, &
                              model_args%DelAB, config%ionvar, model_args%ReIm, model_args%resp_matr, &
                              ReGcont, ImGcont, ReGrev, ImGrev, ReGpiv, ImGpiv, ReGion, ImGion)
        ! open (unit = 21, file = 'fort.21', status='replace', action = 'write')
        ! do j = 1, nf
        !    write(21,*) '-----------------------------', j
        !    do i = 1, nex
        !       write(21,*) i, ReGcont(i,j),ImGcont(i,j),ReGrev(i,j),ImGrev(i,j),&
        !                       ReGpiv(i,j),ImGpiv(i,j),ReGion(i,j),ImGion(i,j)
        !    end do
        ! end do
        ! close(21)
     else 
        call components(config%nex, earx, nf, real(config%flo), real(config%fhi), model_args%nlp, arrays%contx, &
                        tauso, gso, ReW0, ImW0, ReW1, ImW1, ReW2, ImW2, ReW3, ImW3, &
                        real(model_args%h), real(model_args%zcos), real(model_args%Gamma), &
                        real(model_args%eta), model_args%beta_p, &
                        model_args%boost, model_args%g, model_args%DelAB, config%ionvar, &
                        ReScont, ImScont, ReSrev, ImSrev, ReSpiv, ImSpiv, ReSion, ImSion) 
        ! do i = 1, nex
        !    write(20,*) ReSion(i),ImSion(i)
        ! enddo

        !Include the effects of absorRTion in each model component matrix           
       do j = 1, nf
            do i = 1, config%nex
                ReScont(i,j) = ReScont(i,j) * absorbx(i)
                ImScont(i,j) = ImScont(i,j) * absorbx(i)
                ReSrev(i,j)  = ReSrev(i,j)  * absorbx(i)
                ImSrev(i,j)  = ImSrev(i,j)  * absorbx(i)
                ReSpiv(i,j)  = ReSpiv(i,j)  * absorbx(i)
                ImSpiv(i,j)  = ImSpiv(i,j)  * absorbx(i)
                ReSion(i,j)  = ReSion(i,j)  * absorbx(i)
                ImSion(i,j)  = ImSion(i,j)  * absorbx(i)
            end do
        end do 
        !Calculate raw cross-spectrum from S(E,\nu) and the reference band parameters, for each component separately
        if (model_args%ReIm .gt. 0.0) then
            call propercross(config%nex, nf, earx, ReScont, ImScont, ReGcont, ImGcont, model_args%resp_matr)
            call propercross(config%nex, nf, earx, ReSrev,  ImSrev,  ReGrev,  ImGrev,  model_args%resp_matr)     
            call propercross(config%nex, nf, earx, ReSpiv,  ImSpiv,  ReGpiv,  ImGpiv,  model_args%resp_matr)   
            call propercross(config%nex, nf, earx, ReSion,  ImSion,  ReGion,  ImGion,  model_args%resp_matr)
        else
            call propercross_NOmatrix(config%nex, nf, earx, ReScont, ImScont, ReGcont, ImGcont)
            call propercross_NOmatrix(config%nex, nf, earx, ReSrev,  ImSrev,  ReGrev,  ImGrev)
            call propercross_NOmatrix(config%nex, nf, earx, ReSpiv,  ImSpiv,  ReGpiv,  ImGpiv)
            call propercross_NOmatrix(config%nex, nf, earx, ReSion,  ImSion,  ReGion,  ImGion)
        endif   
    endif           
    !Apply phase correction parameter to the cross-spectral model (for bad calibration)
    do j = 1,nf
        do i = 1, config%nex
            tempRe = ReGcont(i,j)
            tempIm = ImGcont(i,j)
            ReGcont(i,j) = cos(model_args%DelA) * tempRe - sin(model_args%DelA) * tempIm
            ImGcont(i,j) = cos(model_args%DelA) * tempIm + sin(model_args%DelA) * tempRe

            tempRe = ReGrev(i,j)
            tempIm = ImGrev(i,j)
            ReGrev(i,j) = cos(model_args%DelA) * tempRe - sin(model_args%DelA) * tempIm
            ImGrev(i,j) = cos(model_args%DelA) * tempIm + sin(model_args%DelA) * tempRe

            tempRe = ReGpiv(i,j)
            tempIm = ImGpiv(i,j)
            ReGpiv(i,j) = cos(model_args%DelA) * tempRe - sin(model_args%DelA) * tempIm
            ImGpiv(i,j) = cos(model_args%DelA) * tempIm + sin(model_args%DelA) * tempRe

            tempRe = ReGion(i,j)
            tempIm = ImGion(i,j)
            ReGion(i,j) = cos(model_args%DelA) * tempRe - sin(model_args%DelA) * tempIm
            ImGion(i,j) = cos(model_args%DelA) * tempIm + sin(model_args%DelA) * tempRe
        end do
    end do
    !Calculate frequency-averaged spectra 
    ReGcont_bar = 0.0
    ImGcont_bar = 0.0
    ReGrev_bar = 0.0
    ImGrev_bar = 0.0
    ReGpiv_bar = 0.0
    ImGpiv_bar = 0.0
    ReGion_bar = 0.0
    ImGion_bar = 0.0
    fc = 0.5d0 * ( floHz + fhiHz )   
    fac = 2.302585* fc**2 * log10(fhiHz/floHz) / ((fhiHz-floHz) * real(nf))
    do j = 1,nf
        f = floHz * (fhiHz/floHz)**( (real(j)-0.5) / real(nf) )
        do i = 1,config%nex 
            ReGcont_bar(i) = ReGcont_bar(i) + ReGcont(i,j) / f
            ImGcont_bar(i) = ImGcont_bar(i) + ImGcont(i,j) / f
            ReGrev_bar(i) = ReGrev_bar(i) + ReGrev(i,j) / f
            ImGrev_bar(i) = ImGrev_bar(i) + ImGrev(i,j) / f
            ReGpiv_bar(i) = ReGpiv_bar(i) + ReGpiv(i,j) / f
            ImGpiv_bar(i) = ImGpiv_bar(i) + ImGpiv(i,j) / f
            ReGion_bar(i) = ReGion_bar(i) + ReGion(i,j) / f
            ImGion_bar(i) = ImGion_bar(i) + ImGion(i,j) / f
        end do
    end do
    ReGcont_bar = ReGcont_bar * fac
    ImGcont_bar = ImGcont_bar * fac
    ReGpiv_bar = ReGpiv_bar * fac
    ImGpiv_bar = ImGpiv_bar * fac 
    ReGrev_bar = ReGrev_bar * fac
    ImGrev_bar = ImGrev_bar * fac
    ReGion_bar = ReGion_bar * fac
    ImGion_bar = ImGion_bar * fac      
    !calculate energy array for output; for the energy we take average energy of the bin
    do i=1,ne 
        ener(i) = (ear(i)+ear(i-1))/2.   
    end do    
       
    path = 'Output/PivPL.dat'
    open (unit = 11, file = path, status='replace', action = 'write')
    
    path = 'Output/Reverb.dat'
    open (unit = 12, file = path, status='replace', action = 'write') 

    path = 'Output/PivRef.dat'
    open (unit = 13, file = path, status='replace', action = 'write')
    
    path = 'Output/IonVar.dat'
    open (unit = 14, file = path, status='replace', action = 'write') 

    if (abs(model_args%ReIm) .le. 4) then
        call crebin(config%nex,earx,ReGcont_bar,ImGcont_bar,ne,ear,ReScont_print,ImScont_print)
        call crebin(config%nex,earx,ReGrev_bar,ImGrev_bar,ne,ear,ReSrev_print,ImSrev_print)
        call crebin(config%nex,earx,ReGpiv_bar,ImGpiv_bar,ne,ear,ReSpiv_print,ImSpiv_print)
        call crebin(config%nex,earx,ReGion_bar,ImGion_bar,ne,ear,ReSion_print,ImSion_print)
        if (abs(model_args%ReIm) .eq. 1 ) then         !Real part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), ReScont_print(i)/dE
                write (12,*) ener(i), ReSrev_print(i)/dE
                write (13,*) ener(i), ReSpiv_print(i)/dE
                write (14,*) ener(i), ReSion_print(i)/dE
            end do    
        else if (abs(model_args%ReIm) .eq. 2) then     !Imaginary part
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), ImScont_print(i)/dE
                write (12,*) ener(i), ImSrev_print(i)/dE
                write (13,*) ener(i), ImSpiv_print(i)/dE
                write (14,*) ener(i), ImSion_print(i)/dE
            end do
        else if (abs(model_args%ReIm) .eq. 3) then     !Modulus
            do i = 1,ne 
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), sqrt(ReScont_print(i)**2 + ImScont_print(i)**2)/dE
                write (12,*) ener(i), sqrt(ReSrev_print(i)**2 + ImSrev_print(i)**2)/dE
                write (13,*) ener(i), sqrt(ReSpiv_print(i)**2 + ImSpiv_print(i)**2)/dE
                write (14,*) ener(i), sqrt(ReSion_print(i)**2 + ImSion_print(i)**2)/dE
            end do
        else if (abs(model_args%ReIm) .eq. 4) then     !Time lag (s)
            do i = 1,ne
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), atan2(ImScont_print(i),ReScont_print(i)) / ( 2.0*pi*fc )
                write (12,*) ener(i), atan2(ImSrev_print(i),ReSrev_print(i)) / ( 2.0*pi*fc )
                write (13,*) ener(i), atan2(ImSpiv_print(i),ReSpiv_print(i)) / ( 2.0*pi*fc ) 
                write (14,*) ener(i), atan2(ImSion_print(i),ReSion_print(i)) / ( 2.0*pi*fc )
            end do
        end if
    else
        call cfoldandbin(config%nex,earx,ReGcont_bar,ImGcont_bar,ne,ear,ReScont_print,ImScont_print,model_args%resp_matr)
        call cfoldandbin(config%nex,earx,ReGpiv_bar,ImGpiv_bar,ne,ear,ReSpiv_print,ImSpiv_print,model_args%resp_matr)
        call cfoldandbin(config%nex,earx,ReGrev_bar,ImGrev_bar,ne,ear,ReSrev_print,ImSrev_print,model_args%resp_matr)
        call cfoldandbin(config%nex,earx,ReGion_bar,ImGion_bar,ne,ear,ReSion_print,ImSion_print,model_args%resp_matr)
        if (abs(model_args%ReIm) .eq. 5) then          !Modulus
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), sqrt(ReScont_print(i)**2 + ImScont_print(i)**2)/dE
                write (13,*) ener(i), sqrt(ReSpiv_print(i)**2 + ImSpiv_print(i)**2)/dE
                write (12,*) ener(i), sqrt(ReSrev_print(i)**2 + ImSrev_print(i)**2)/dE
                write (14,*) ener(i), sqrt(ReSion_print(i)**2 + ImSion_print(i)**2)/dE
            end do
        else if (abs(model_args%ReIm) .eq. 6) then     !Time lag (s)
            do i = 1, ne
                dE = ear(i) - ear(i-1)
                write (11,*) ener(i), atan2(ImScont_print(i),ReScont_print(i)) / ( 2.0*pi*fc )
                write (12,*) ener(i), atan2(ImSrev_print(i),ReSrev_print(i)) / ( 2.0*pi*fc ) 
                write (13,*) ener(i), atan2(ImSpiv_print(i),ReSpiv_print(i)) / ( 2.0*pi*fc )
                write (14,*) ener(i), atan2(ImSion_print(i),ReSion_print(i)) / ( 2.0*pi*fc )
            end do
        end if
    end if

    close(11)
    close(12)
    close(13)
    close(14)  

    if (allocated( ReScont )) deallocate( ReScont )
    if (allocated( ImScont )) deallocate( ImScont )
    if (allocated( ReSrev  )) deallocate( ReSrev  )
    if (allocated( ImSrev  )) deallocate( ImSrev  )
    if (allocated( ReSpiv  )) deallocate( ReSpiv  )
    if (allocated( ImSpiv  )) deallocate( ImSpiv  )
    if (allocated( ReSion  )) deallocate( ReSion  )
    if (allocated( ImSion  )) deallocate( ImSion  )
    if (allocated( ReGcont )) deallocate( ReGcont )
    if (allocated( ImGcont )) deallocate( ImGcont )
    if (allocated( ReGrev  )) deallocate( ReGrev  )
    if (allocated( ImGrev  )) deallocate( ImGrev  )
    if (allocated( ReGpiv  )) deallocate( ReGpiv  )
    if (allocated( ImGpiv  )) deallocate( ImGpiv  )    
    if (allocated( ReGion  )) deallocate( ReGion  )
    if (allocated( ImGion  )) deallocate( ImGion  )
    if (allocated( ReGcont_bar )) deallocate( ReGcont_bar )
    if (allocated( ImGcont_bar )) deallocate( ImGcont_bar )
    if (allocated( ReGrev_bar  )) deallocate( ReGrev_bar  )
    if (allocated( ImGrev_bar  )) deallocate( ImGrev_bar  )
    if (allocated( ReGpiv_bar  )) deallocate( ReGpiv_bar  )
    if (allocated( ImGpiv_bar  )) deallocate( ImGpiv_bar  )
    if (allocated( ReGion_bar  )) deallocate( ReGion_bar  )
    if (allocated( ImGion_bar  )) deallocate( ImGion_bar  )
    if (allocated( ReScont_print )) deallocate( ReScont_print )
    if (allocated( ImScont_print )) deallocate( ImScont_print )
    if (allocated( ReSpiv_print  )) deallocate( ReSpiv_print  )
    if (allocated( ImSpiv_print  )) deallocate( ImSpiv_print  )
    if (allocated( ener          )) deallocate( ener          )
    if (allocated( ReSrev_print  )) deallocate( ReSrev_print  )
    if (allocated( ImSrev_print  )) deallocate( ImSrev_print  )
    if (allocated( ReSion_print  )) deallocate( ReSion_print  )
    if (allocated( ImSion_print  )) deallocate( ImSion_print  )
	    
    return  
end subroutine write_ComponentS

subroutine components(nex,earx,nf,flo,fhi,nlp,contx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                      h,z,Gamma,eta,beta_p,boost,g,DelAB,ionvar,ReScont,ImScont,ReSrev,ImSrev,&
                      ReSpiv,ImSpiv,ReSion,ImSion)
    ! Calculates the FT of the spectrum components before multiplying by the absorRTion model
    ! This is essentially the same as S, but it returns the FT for each component that contributes to the lags: 
    ! OutputssinD
    ! Scont(1:nex,1:nf)     Continuum pivoting cross spectrum for all LPs
    ! Sreverb(1:nex,1:nf)   Continuum pivoting+reverberation cross spectrum for all LPs
    ! Spivot(1:nex,1:nf)    Continuum and reflection pivoting+reverberation cross spectrum for all LPs
    ! Sion(1:nex,1:nf)      Continuum pivoting+ionization fluctuation cross spectrum for all LPs

    use constants
    implicit none
    integer, intent(IN) :: nex,nf,ionvar,nlp
    real   , intent(IN) :: earx(0:nex),contx(nex,nlp)
    real   , intent(IN) :: g(nlp),DelAB(nlp),boost,z,gso(nlp),Gamma,eta,h(nlp),tauso(nlp), beta_p, flo, fhi
    real   , intent(INOUT) :: ReW0(nlp,nex,nf),ImW0(nlp,nex,nf),ReW1(nlp,nex,nf),ImW1(nlp,nex,nf)
    real   , intent(INOUT) :: ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)
    real   , intent(INOUT) :: ReScont(nex,nf),ImScont(nex,nf),ReSrev(nex,nf),ImSrev(nex,nf)
    real   , intent(INOUT) :: ReSpiv(nex,nf),ImSpiv(nex,nf),ReSion(nex,nf),ImSion(nex,nf)
    real E,fac,beta,f,phase_d,phase_p,tau_d,tau_p
    real, parameter :: pi_single = acos(-1.0)
    real corr, contx_sum(nex)
    complex, dimension(:,:), allocatable :: Scont,Sreverb,Spivot,Sion
    ! complex Stemp,Scont(nex,nf),Sreverb(nex,nf),Spivot(nex,nf),Sion(nex,nf)
    complex Stemp,cexp_d,cexp_p,cexp_phi,W0,W1,W2,W3
    integer i,j,m

    if(.not. allocated(Scont  )) allocate(Scont   (nex,nf) )
    if(.not. allocated(Sreverb)) allocate(Sreverb (nex,nf) )
    if(.not. allocated(Spivot )) allocate(Spivot  (nex,nf) )
    if(.not. allocated(Sion   )) allocate(Sion    (nex,nf) )  
    
    Scont = 0.
    Sreverb = 0.
    Spivot = 0.
    Sion = 0.
       
    !TBD initialize new lags to 0/1 as appropriate here
    phase_d = 0.
    phase_p = 0.
    tau_d = 0.
    tau_p = 0.
    
    do m=1,nlp 
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
            tau_p = (h(m) - h(1))/(beta_p) !I think this is fine, but may need an extra factor c? double check the sign
        end if
        do j = 1,nf
            f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )
            do i = 1, nex
                E   = 0.5 * ( earx(i) + earx(i-1) )
                fac = log(gso(m)/((1.0+z)*E))
                !set up phase factors
                if (m .gt. 1) then
                    phase_d = 2.0*pi_single*tau_d*f
                    phase_p = 2.0*pi_single*tau_p*f
                endif    
                cexp_d = cmplx(cos(phase_d),sin(phase_d))
                cexp_p = cmplx(cos(phase_p),sin(phase_p)) 
                cexp_phi = cmplx(cos(DelAB(m)),sin(DelAB(m)))             
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = cexp_p*(g(m)*cexp_phi*fac*cexp_d*contx(i,m) + cexp_d*contx(i,m))
                Scont(i,j) = Scont(i,j) + Stemp
                Stemp = cexp_p*(cexp_d*contx(i,m)+W0)
                Sreverb(i,j) = Sreverb(i,j) + Stemp
                Stemp =  g(m)*cexp_phi*(fac*cexp_d*contx(i,m)+W1+W2) + cexp_d*contx(i,m)
                Stemp = cexp_p*Stemp
                Spivot(i,j) = Spivot(i,j) + Stemp
                Stemp = cexp_p*(cexp_d*contx(i,m)+W3)
                Sion(i,j) = Sion(i,j) + Stemp
            enddo
        enddo 
    end do

    ReScont = real(Scont)
    ImScont = aimag(Scont)
    ReSrev = real(Sreverb)
    ImSrev = aimag(Sreverb)
    ReSpiv = real(Spivot)
    ImSpiv = aimag(Spivot)
    ReSion = real(Sion)
    ImSion = aimag(Sion)
    if(allocated(Scont  )) deallocate(Scont  )
    if(allocated(Sreverb)) deallocate(Sreverb)
    if(allocated(Spivot )) deallocate(Spivot )
    if(allocated(Sion   )) deallocate(Sion   )  
    
    return
end subroutine

subroutine components_nocoh(nex,earx,nf,flo,fhi,nlp,contx,absorbx,tauso,gso,ReW0,ImW0,ReW1,ImW1,ReW2,ImW2,ReW3,ImW3,&
                            h,z,Gamma,eta,boost,g,DelAB,ionvar,ReIm,resp_matr,ReGcont,ImGcont,ReGrev,ImGrev,&
                            ReGpiv,ImGpiv,ReGion,ImGion)
    use constants
    implicit none
    integer, intent(IN) :: nex,nf,ionvar,nlp,ReIm,resp_matr
    real, intent(IN) :: earx(0:nex),contx(nex,nlp),absorbx(nex)
    real, intent(IN) :: ReW0(nlp,nex,nf),ImW0(nlp,nex,nf),ReW1(nlp,nex,nf),ImW1(nlp,nex,nf)
    real, intent(IN) :: ReW2(nlp,nex,nf),ImW2(nlp,nex,nf),ReW3(nlp,nex,nf),ImW3(nlp,nex,nf)
    real, intent(IN) :: g(nlp),DelAB(nlp),boost,z,gso(nlp),Gamma,eta,h(nlp),tauso(nlp),flo,fhi
    real E,fac,beta,f,phase_d,tau_d
    real, parameter :: pi_single = acos(-1.0)
    real contx_sum(nex),corr
    !these are the component arrays for each lamp post separately, before the cross spectrum
    real, dimension(:,:,:), allocatable :: ReSpiv,ImSpiv,ReScont,ImScont,ReSrev,ImSrev, ReSion,ImSion
    ! real ReScont(nlp,nex,nf),ImScont(nlp,nex,nf),ReSrev(nlp,nex,nf),ImSrev(nlp,nex,nf)
    !these are the  component arrays for each lamp post separately, after the cross spectrum

    real, dimension(:,:,:), allocatable :: ReGcont_temp,ImGcont_temp,ReGrev_temp,ImGrev_temp
    real, dimension(:,:,:), allocatable :: ReGpiv_temp,ImGpiv_temp,ReGion_temp,ImGion_temp
    ! real ReGcont_temp(nlp,nex,nf),ImGcont_temp(nlp,nex,nf),ReGrev_temp(nlp,nex,nf),ImGrev_temp(nlp,nex,nf)
    ! real ReGpiv_temp(nlp,nex,nf),ImGpiv_temp(nlp,nex,nf),ReGion_temp(nlp,nex,nf),ImGion_temp(nlp,nex,nf)    
    !these are the arrays containing the sum of the two lamp posts for zero coherence 
    real, intent(INOUT) :: ReGcont(nex,nf),ImGcont(nex,nf),ReGrev(nex,nf),ImGrev(nex,nf)
    real, intent(INOUT) :: ReGpiv(nex,nf),ImGpiv(nex,nf),ReGion(nex,nf),ImGion(nex,nf)
    !this complex stuff is purely to make writing the phases easier
    complex, dimension(:,:,:), allocatable :: Scont,Sreverb,Spivot,Sion
    ! complex Stemp,Scont(nlp,nex,nf),Sreverb(nlp,nex,nf),Spivot(nlp,nex,nf),Sion(nlp,nex,nf)
    complex Stemp, cexp_d,cexp_phi,W0,W1,W2,W3
    integer i,j,m

    if(.not. allocated(ReGcont_temp)) allocate(ReGcont_temp(nlp,nex,nf) )
    if(.not. allocated(ImGcont_temp)) allocate(ImGcont_temp(nlp,nex,nf) )
    if(.not. allocated(ReGrev_temp )) allocate(ReGrev_temp (nlp,nex,nf) )
    if(.not. allocated(ImGrev_temp )) allocate(ImGrev_temp (nlp,nex,nf) )  
    if(.not. allocated(ReGpiv_temp )) allocate(ReGpiv_temp (nlp,nex,nf) )
    if(.not. allocated(ImGpiv_temp )) allocate(ImGpiv_temp (nlp,nex,nf) )
    if(.not. allocated(ReGion_temp )) allocate(ReGion_temp (nlp,nex,nf) )
    if(.not. allocated(ImGion_temp )) allocate(ImGion_temp (nlp,nex,nf) )  
    if(.not. allocated(ReScont     )) allocate(ReScont     (nlp,nex,nf) )
    if(.not. allocated(ImScont     )) allocate(ImScont     (nlp,nex,nf) )
    if(.not. allocated(ReSrev      )) allocate(ReSrev      (nlp,nex,nf) )
    if(.not. allocated(ImSrev      )) allocate(ImSrev      (nlp,nex,nf) )  
    if(.not. allocated(ReSpiv      )) allocate(ReSpiv      (nlp,nex,nf) )
    if(.not. allocated(ImSpiv      )) allocate(ImSpiv      (nlp,nex,nf) )
    if(.not. allocated(ReSion      )) allocate(ReSion      (nlp,nex,nf) )
    if(.not. allocated(ImSion      )) allocate(ImSion      (nlp,nex,nf) )  

    if(.not. allocated(Scont       )) allocate(Scont       (nlp,nex,nf) )
    if(.not. allocated(Sreverb     )) allocate(Sreverb     (nlp,nex,nf) )
    if(.not. allocated(Spivot      )) allocate(Spivot      (nlp,nex,nf) )
    if(.not. allocated(Sion        )) allocate(Sion        (nlp,nex,nf) )  

    Scont   = 0.
    Sreverb = 0.
    Spivot = 0.
    Sion    = 0.
    
    ReGcont = 0.
    ImGcont = 0.
    ReGrev = 0.
    ImGrev = 0.
    ReGpiv = 0.
    ImGpiv = 0.
    ReGion = 0.
    ImGion = 0.
    
    phase_d = 0.
    tau_d = 0.
    
    do m=1,nlp 
        do j = 1,nf 
            f = flo * (fhi/flo)**(  (real(j)-0.5) / real(nf) )
            do i = 1,nex
                E   = 0.5 * ( earx(i) + earx(i-1) )
                fac = log(gso(m)/((1.0+z)*E))
                if (m .gt. 1) then
                    tau_d = tauso(m)-tauso(1)
                    phase_d = 2.0*pi_single*tau_d*f  
                endif
                cexp_d = cmplx(cos(phase_d),sin(phase_d))     
                cexp_phi = cmplx(cos(DelAB(m)),sin(DelAB(m)))
                !set up transfer functions 
                W0 = boost * cmplx(ReW0(m,i,j),ImW0(m,i,j))
                W1 = boost * cmplx(ReW1(m,i,j),ImW1(m,i,j))
                W2 = boost * cmplx(ReW2(m,i,j),ImW2(m,i,j))                       
                W3 = ionvar * boost * cmplx(ReW3(m,i,j),ImW3(m,i,j))
                !calculate complex covariance
                !note: the reason we use complex here is to ease the calculations 
                !when we add all the extra phases from the double lamp post 
                Stemp = (g(m)*cexp_phi*fac*cexp_d*contx(i,m) + cexp_d*contx(i,m))*absorbx(i)
                Scont(m,i,j) = Scont(m,i,j) + Stemp
                Stemp = (cexp_d*contx(i,m)+W0)*absorbx(i)
                Sreverb(m,i,j) = Sreverb(m,i,j) + Stemp
                Stemp =  (g(m)*cexp_phi*(fac*cexp_d*contx(i,m)+W1+W2) + cexp_d*contx(i,m))*absorbx(i)
                Spivot(m,i,j) = Spivot(m,i,j) + Stemp
                Stemp = (cexp_d*contx(i,m)+W3)*absorbx(i)
                Sion(m,i,j) = Sion(m,i,j) + Stemp
            enddo 
        enddo    
    end do
    
    ReScont = real(Scont)
    ImScont = aimag(Scont)
    ReSrev = real(Sreverb)
    ImSrev = aimag(Sreverb)
    ReSpiv = real(Spivot)
    ImSpiv = aimag(Spivot)
    ReSion = real(Sion)
    ImSion = aimag(Sion)

    ! open (unit = 22, file = 'fort.22', status='replace', action = 'write')
    ! do m = 1, nlp
    !    do j = 1, nf
    !       do i = 1, nex
    !          write(22,*) m,j,i, Sion(m,i,j), ReSion(m,i,j),ImSion(m,i,j)
    !       end do
    !    end do
    ! enddo
    ! close(22)
    
    do m=1,nlp 
        if (ReIm .gt. 0.0) then 
            call propercross(nex,nf,earx,ReScont(m,:,:),ImScont(m,:,:),ReGcont_temp(m,:,:),ImGcont_temp(m,:,:),resp_matr)
            call propercross(nex,nf,earx,ReSrev(m,:,:),ImSrev(m,:,:),ReGrev_temp(m,:,:),ImGrev_temp(m,:,:),resp_matr)
            call propercross(nex,nf,earx,ReSpiv(m,:,:),ImSpiv(m,:,:),ReGpiv_temp(m,:,:),ImGpiv_temp(m,:,:),resp_matr)
            call propercross(nex,nf,earx,ReSion(m,:,:),ImSion(m,:,:),ReGion_temp(m,:,:),ImGion_temp(m,:,:),resp_matr)
        else 
            call propercross_NOmatrix(nex,nf,earx,ReScont(m,:,:),ImScont(m,:,:),ReGcont_temp(m,:,:),ImGcont_temp(m,:,:))
            call propercross_NOmatrix(nex,nf,earx,ReSrev(m,:,:),ImSrev(m,:,:),ReGrev_temp(m,:,:),ImGrev_temp(m,:,:))
            call propercross_NOmatrix(nex,nf,earx,ReSpiv(m,:,:),ImSpiv(m,:,:),ReGpiv_temp(m,:,:),ImGpiv_temp(m,:,:))
            call propercross_NOmatrix(nex,nf,earx,ReSion(m,:,:),ImSion(m,:,:),ReGion_temp(m,:,:),ImGion_temp(m,:,:))
        endif
        do j=1,nf 
            do i=1,nex 
                if (m .gt. 1) then
                    ReGcont_temp(m,i,j) = eta**2.*ReGcont_temp(m,i,j)
                    ImGcont_temp(m,i,j) = eta**2.*ImGcont_temp(m,i,j)
                    ReGrev_temp(m,i,j) = eta**2.*ReGrev_temp(m,i,j)
                    ImGrev_temp(m,i,j) = eta**2.*ImGrev_temp(m,i,j)
                    ReGpiv_temp(m,i,j) = eta**2.*ReGpiv_temp(m,i,j)
                    ImGpiv_temp(m,i,j) = eta**2.*ImGpiv_temp(m,i,j)
                    ReGion_temp(m,i,j) = eta**2.*ReGion_temp(m,i,j)
                    ImGion_temp(m,i,j) = eta**2.*ImGion_temp(m,i,j)
                endif
                ReGcont(i,j) = ReGcont(i,j) + ReGcont_temp(m,i,j)
                ImGcont(i,j) = ImGcont(i,j) + ImGcont_temp(m,i,j)
                ReGrev(i,j) = ReGrev(i,j) + ReGrev_temp(m,i,j)
                ImGrev(i,j) = ImGrev(i,j) + ImGrev_temp(m,i,j)
                ReGpiv(i,j) = ReGpiv(i,j) + ReGpiv_temp(m,i,j)
                ImGpiv(i,j) = ImGpiv(i,j) + ImGpiv_temp(m,i,j)
                ReGion(i,j) = ReGion(i,j) + ReGion_temp(m,i,j)
                ImGion(i,j) = ImGion(i,j) + ImGion_temp(m,i,j)
            end do
        end do
    end do          

    if(allocated(ReGcont_temp)) deallocate(ReGcont_temp)
    if(allocated(ImGcont_temp)) deallocate(ImGcont_temp)
    if(allocated(ReGrev_temp )) deallocate(ReGrev_temp )
    if(allocated(ImGrev_temp )) deallocate(ImGrev_temp )  
    if(allocated(ReGpiv_temp )) deallocate(ReGpiv_temp )
    if(allocated(ImGpiv_temp )) deallocate(ImGpiv_temp )
    if(allocated(ReGion_temp )) deallocate(ReGion_temp )
    if(allocated(ImGion_temp )) deallocate(ImGion_temp )  
    if(allocated(ReScont     )) deallocate(ReScont     )
    if(allocated(ImScont     )) deallocate(ImScont     )
    if(allocated(ReSrev      )) deallocate(ReSrev      )
    if(allocated(ImSrev      )) deallocate(ImSrev      )  
    if(allocated(ReSpiv      )) deallocate(ReSpiv      )
    if(allocated(ImSpiv      )) deallocate(ImSpiv      )
    if(allocated(ReSion      )) deallocate(ReSion      )
    if(allocated(ImSion      )) deallocate(ImSion      )  
    if(allocated(Scont       )) deallocate(Scont       )
    if(allocated(Sreverb     )) deallocate(Sreverb     )
    if(allocated(Spivot      )) deallocate(Spivot      )
    if(allocated(Sion        )) deallocate(Sion        )  
                     
    return
end subroutine
                      
end module m_components


