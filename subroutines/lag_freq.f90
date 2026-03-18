subroutine lag_freq(config, model_args, arrays, absorbx)

! Calculates the energy-averaged cross spectrum in the energy bins of
! interest.
!
! Input: continuum model (contx is f(E) in the papers) and the reflection
!        transfer functions
! Output: arrays%ReGbar(nf), arrays%ImGbar(nf) after multiplying by the
!         absorption model
! These inputs and outputs are all in terms of (dN/dE)*dE; i.e. photar
    use common_types
    use conv_mod, only: nex
    use gr_continuum, only: tauso, gso
    use constants
    implicit none
    type(t_config),          intent(in)    :: config
    type(t_model_arguments), intent(in)    :: model_args
    type(t_arrays),          intent(inout) :: arrays
    real,                    intent(in)    :: absorbx(nex)

    integer             :: Ea1, Ea2, Eb1, Eb2
    real                :: gslope, ABslope
    real                :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
    real                :: E, fac, TempReG, TempImG
    real                :: f, DelAB_nu, g_nu
    real                :: tau_d, phase_d, tau_p, phase_p, beta
    complex             :: W0, W1, W2, W3, Sraw, cexp_p, cexp_d, cexp_phi,     &
                           Stemp
    integer             :: i, j, m

    call energy_bounds(config, Ea1, Ea2, Eb1, Eb2)

    gslope = 1.
    ABslope = 1.

    if(model_args%nlp .gt. 1) then
        arrays%ReW0(2,:,:) = real(model_args%eta) * arrays%ReW0(2,:,:)
        arrays%ImW0(2,:,:) = real(model_args%eta) * arrays%ImW0(2,:,:)
        arrays%ReW1(2,:,:) = real(model_args%eta) * arrays%ReW1(2,:,:)
        arrays%ImW1(2,:,:) = real(model_args%eta) * arrays%ImW1(2,:,:)
        arrays%ReW2(2,:,:) = real(model_args%eta) * arrays%ReW2(2,:,:)
        arrays%ImW2(2,:,:) = real(model_args%eta) * arrays%ImW2(2,:,:)
        arrays%ReW3(2,:,:) = real(model_args%eta) * arrays%ReW3(2,:,:)
        arrays%ImW3(2,:,:) = real(model_args%eta) * arrays%ImW3(2,:,:)
    endif

    !Now calculate the cross-spectrum (/complex covariance), including
    !absorption
    do j = 1, config%nf
        f = real(config%flo) * (real(config%fhi) / real(config%flo))            &
            ** ((real(j) - 0.5) / real(config%nf))
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
            E = 0.5 * (arrays%earx(i) + arrays%earx(i - 1))
            do m = 1, model_args%nlp
                DelAB_nu = model_args%DelAB(m)                                 &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**ABslope
                g_nu = model_args%g(m)                                         &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**gslope
                fac = log(real(gso(m))                                         &
                    / ((1.0 + real(model_args%zcos)) * E))
                !set up phase factors
                if (m .gt. 1) then
                    tau_d = real(tauso(m) - tauso(1))
                    tau_p = (real(model_args%h(m))                             &
                        - real(model_args%h(1))) / (model_args%beta_p)
                    phase_d = 2. * pi * tau_d * f
                    phase_p = 2. * pi * tau_p * f
                endif
                cexp_d = cmplx(cos(phase_d), sin(phase_d))
                cexp_p = cmplx(cos(phase_p), sin(phase_p))
                cexp_phi = cmplx(cos(DelAB_nu), sin(DelAB_nu))
                !set up transfer functions
                W0 = model_args%boost                                          &
                    * cmplx(arrays%ReW0(m, i, j), arrays%ImW0(m, i, j))
                W1 = model_args%boost                                          &
                    * cmplx(arrays%ReW1(m, i, j), arrays%ImW1(m, i, j))
                W2 = model_args%boost                                          &
                    * cmplx(arrays%ReW2(m, i, j), arrays%ImW2(m, i, j))
                W3 = config%ionvar * model_args%boost                          &
                    * cmplx(arrays%ReW3(m, i, j), arrays%ImW3(m, i, j))
                !calculate complex covariance
                Stemp = g_nu * cexp_phi                                        &
                    * (W1 + W2 + fac * cexp_d * arrays%contx(i, m))
                Stemp = Stemp + W0 + W3 + cexp_d * arrays%contx(i, m)
                Stemp = cexp_p * Stemp
                Sraw = Sraw + Stemp
            end do
            !separate into real/imaginary parts
            ReGrawEa = ReGrawEa + real(Sraw) * absorbx(i)
            ImGrawEa = ImGrawEa + aimag(Sraw) * absorbx(i)
        end do

        do i = Eb1, Eb2
            phase_d = 0.
            phase_p = 0.
            tau_d = 0.
            tau_p = 0.
            Sraw = 0.
            E = 0.5 * (arrays%earx(i) + arrays%earx(i - 1))
            do m = 1, model_args%nlp
                DelAB_nu = model_args%DelAB(m)                                 &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**ABslope
                g_nu = model_args%g(m)                                         &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**gslope
                fac = log(real(gso(m))                                         &
                    / ((1.0 + real(model_args%zcos)) * E))
                !set up phase factors
                if (m .gt. 1) then
                    tau_d = real(tauso(m) - tauso(1))
                    tau_p = (real(model_args%h(m))                             &
                        - real(model_args%h(1))) / (model_args%beta_p)
                    phase_d = 2. * pi * tau_d * f
                    phase_p = 2. * pi * tau_p * f
                endif
                cexp_d = cmplx(cos(phase_d), sin(phase_d))
                cexp_p = cmplx(cos(phase_p), sin(phase_p))
                cexp_phi = cmplx(cos(DelAB_nu), sin(DelAB_nu))
                !set up transfer functions
                W0 = model_args%boost                                          &
                    * cmplx(arrays%ReW0(m, i, j), arrays%ImW0(m, i, j))
                W1 = model_args%boost                                          &
                    * cmplx(arrays%ReW1(m, i, j), arrays%ImW1(m, i, j))
                W2 = model_args%boost                                          &
                    * cmplx(arrays%ReW2(m, i, j), arrays%ImW2(m, i, j))
                W3 = config%ionvar * model_args%boost                          &
                    * cmplx(arrays%ReW3(m, i, j), arrays%ImW3(m, i, j))
                !calculate complex covariance
                Stemp = g_nu * cexp_phi                                        &
                    * (W1 + W2 + fac * cexp_d * arrays%contx(i, m))
                Stemp = Stemp + W0 + W3 + cexp_d * arrays%contx(i, m)
                Stemp = cexp_p * Stemp
                Sraw = Sraw + Stemp
            end do
            !separate into real/imaginary parts
            ReGrawEb = ReGrawEb + real(Sraw) * absorbx(i)
            ImGrawEb = ImGrawEb + aimag(Sraw) * absorbx(i)
        end do
        !Now cross-spectrum between the two energy bands
        !note: here the conjugate is b
        arrays%ReGbar(j) = (ReGrawEa * ReGrawEb)                              &
            + (ImGrawEa * ImGrawEb)
        arrays%ImGbar(j) = (ReGrawEb * ImGrawEa)                              &
            - (ReGrawEa * ImGrawEb)
    end do

    return
end subroutine lag_freq

subroutine lag_freq_nocoh(config, model_args, arrays, absorbx)

    use common_types
    use conv_mod, only: nex
    use gr_continuum, only: tauso, gso
    use constants
    implicit none
    type(t_config),          intent(in)    :: config
    type(t_model_arguments), intent(in)    :: model_args
    type(t_arrays),          intent(inout) :: arrays
    real,                    intent(in)    :: absorbx(nex)

    integer             :: Ea1, Ea2, Eb1, Eb2
    real                :: gslope, ABslope
    real                :: ReGrawEa, ImGrawEa, ReGrawEb, ImGrawEb
    real                :: E, fac, TempReG, TempImG
    real                :: f, DelAB_nu, g_nu
    real                :: tau_d, phase_d
    complex             :: W0, W1, W2, W3, Sraw, cexp_d, cexp_phi, Stemp
    integer             :: i, j, m

    call energy_bounds(config, Ea1, Ea2, Eb1, Eb2)

    gslope = 1.
    Abslope = 1.

    do m = 1, model_args%nlp
        do j = 1, config%nf
            f = real(config%flo) * (real(config%fhi) / real(config%flo))        &
                ** ((real(j) - 0.5) / real(config%nf))
            ReGrawEa = 0.0
            ImGrawEa = 0.0
            ReGrawEb = 0.0
            ImGrawEb = 0.0
            do i = Ea1, Ea2
                phase_d = 0.
                tau_d = 0.
                E = 0.5 * (arrays%earx(i) + arrays%earx(i - 1))
                DelAB_nu = model_args%DelAB(m)                                 &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**ABslope
                g_nu = model_args%g(m)                                         &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**gslope
                fac = log(real(gso(m))                                         &
                    / ((1.0 + real(model_args%zcos)) * E))
                !set up phase factors
                if (m .gt. 1) then
                    tau_d = real(tauso(m) - tauso(1))
                    phase_d = 2. * pi * tau_d * f
                endif
                cexp_d = cmplx(cos(phase_d), sin(phase_d))
                cexp_phi = cmplx(cos(DelAB_nu), sin(DelAB_nu))
                !set up transfer functions
                W0 = model_args%boost                                          &
                    * cmplx(arrays%ReW0(m, i, j), arrays%ImW0(m, i, j))
                W1 = model_args%boost                                          &
                    * cmplx(arrays%ReW1(m, i, j), arrays%ImW1(m, i, j))
                W2 = model_args%boost                                          &
                    * cmplx(arrays%ReW2(m, i, j), arrays%ImW2(m, i, j))
                W3 = config%ionvar * model_args%boost                          &
                    * cmplx(arrays%ReW3(m, i, j), arrays%ImW3(m, i, j))
                !calculate complex covariance
                Stemp = g_nu * cexp_phi                                        &
                    * (W1 + W2 + fac * cexp_d * arrays%contx(i, m))            &
                    + W0 + W3 + cexp_d * arrays%contx(i, m)
                !separate into real/imaginary parts
                ReGrawEa = ReGrawEa + real(Stemp) * absorbx(i)
                ImGrawEa = ImGrawEa + aimag(Stemp) * absorbx(i)
            end do

            do i = Eb1, Eb2
                phase_d = 0.
                tau_d = 0.
                E = 0.5 * (arrays%earx(i) + arrays%earx(i - 1))
                DelAB_nu = model_args%DelAB(m)                                 &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**ABslope
                g_nu = model_args%g(m)                                         &
                    * ((arrays%fix(1) + arrays%fix(0)) * 0.5 / f)**gslope
                fac = log(real(gso(m))                                         &
                    / ((1.0 + real(model_args%zcos)) * E))
                !set up phase factors
                if (m .gt. 1) then
                    tau_d = real(tauso(m) - tauso(1))
                    phase_d = 2. * pi * tau_d * f
                endif
                cexp_d = cmplx(cos(phase_d), sin(phase_d))
                cexp_phi = cmplx(cos(DelAB_nu), sin(DelAB_nu))
                !set up transfer functions
                W0 = model_args%boost                                          &
                    * cmplx(arrays%ReW0(m, i, j), arrays%ImW0(m, i, j))
                W1 = model_args%boost                                          &
                    * cmplx(arrays%ReW1(m, i, j), arrays%ImW1(m, i, j))
                W2 = model_args%boost                                          &
                    * cmplx(arrays%ReW2(m, i, j), arrays%ImW2(m, i, j))
                W3 = config%ionvar * model_args%boost                          &
                    * cmplx(arrays%ReW3(m, i, j), arrays%ImW3(m, i, j))
                !calculate complex covariance
                Stemp = g_nu * cexp_phi                                        &
                    * (W1 + W2 + fac * cexp_d * arrays%contx(i, m))            &
                    + W0 + W3 + cexp_d * arrays%contx(i, m)
                !separate into real/imaginary parts
                ReGrawEb = ReGrawEb + real(Stemp) * absorbx(i)
                ImGrawEb = ImGrawEb + aimag(Stemp) * absorbx(i)
            end do
            if(m .eq. 1) then
                arrays%ReGbar(j) = (ReGrawEa * ReGrawEb)                       &
                    + (ImGrawEa * ImGrawEb)
                arrays%ImGbar(j) = (ReGrawEb * ImGrawEa)                       &
                    - (ReGrawEa * ImGrawEb)
            else
                arrays%ReGbar(j) = arrays%ReGbar(j)                            &
                    + real(model_args%eta)**2.                                  &
                    * ((ReGrawEa * ReGrawEb) + (ImGrawEa * ImGrawEb))
                arrays%ImGbar(j) = arrays%ImGbar(j)                            &
                    + real(model_args%eta)**2.                                  &
                    * ((ReGrawEb * ImGrawEa) - (ReGrawEa * ImGrawEb))
            end if
        end do
    end do

    return
end subroutine lag_freq_nocoh

subroutine energy_bounds(config, Ea1, Ea2, Eb1, Eb2)
    use common_types
    use conv_mod, only: nex
    use telematrix
    implicit none
    type(t_config), intent(in) :: config
    integer, intent(out) :: Ea1, Ea2, Eb1, Eb2
    real                :: band1_Elo, band1_Ehi, band2_Elo, band2_Ehi
    real     :: get_env_real, dum

    if( needchans ) then
        band1_Elo = get_env_real("EMIN_REF", 0.0)
        band1_Ehi = get_env_real("EMAX_REF", 0.0)
        band2_Elo = get_env_real("EMIN_REF2", 0.0)
        band2_Ehi = get_env_real("EMAX_REF2", 0.0)
        if (band1_Elo .eq. 0.0) then
            write(*,*)"Enter lower energy in the first band"
            read(*,*) band1_Elo
        endif
        if (band1_Ehi .eq. 0.0) then
            write(*,*)"Enter upper energy in the first band"
            read(*,*) band1_Ehi
        endif
        if (band2_Elo .eq. 0.0) then
            write(*,*)"Enter lower energy in the second band"
            read(*,*) band2_Elo
        endif
        if (band2_Ehi .eq. 0.0) then
            write(*,*)"Enter upper energy in the second band"
            read(*,*) band2_Ehi
        endif
        if( band1_Elo .gt. band1_Ehi )then
           dum  = band1_Elo
           band1_Elo = band1_Ehi
           band1_Ehi = dum
           write(*,*)"Elo1>Ehi1! Switched!"
        end if
        if( band2_Elo .gt. band2_Ehi )then
           dum  = band2_Elo
           band2_Elo = band2_Ehi
           band2_Ehi = dum
           write(*,*)"Elo2>Ehi2! Switched!"
        end if
        Ea1 = ceiling(real(nex) * log10(band1_Elo / config%Emin)               &
            / log10(config%Emax / config%Emin))
        Ea2 = ceiling(real(nex) * log10(band1_Ehi / config%Emin)               &
            / log10(config%Emax / config%Emin))
        Eb1 = ceiling(real(nex) * log10(band2_Elo / config%Emin)               &
            / log10(config%Emax / config%Emin))
        Eb2 = ceiling(real(nex) * log10(band2_Ehi / config%Emin)               &
            / log10(config%Emax / config%Emin))
        needchans = .false.
    end if

    return
end subroutine energy_bounds

