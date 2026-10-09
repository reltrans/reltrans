subroutine sum_continuum_reflection_transfer_functions(config, model_args,     &
    arrays, nex, tauso, gso)
    !> Calculates the FT of the spectrum before multiplying by the absorption
    !> model. This was once called `rawS` and refers to the `S` in the Fourier
    !> transform of the total signal in Equation (8) of Mastroserio et al. 2019.
    !>
    !> The additional arguments are:
    !>
    !>     nex              The length of the fine energy grid.
    !>     tauso            The source-to-disc light travel time.
    !>     gso              The source-to-disc energyshift.
    !>
    !> Writes the output to `arrays%ReGraw` and `arrays%ImGraw`.
    use common_types, only: t_config, t_model_arguments, t_arrays
    use rtconstants, only: pi, wp
    implicit none
    type(t_arrays), intent(inout) :: arrays
    type(t_config), intent(in) :: config
    type(t_model_arguments), intent(in) :: model_args
    real(wp), intent(in) :: tauso(2), gso(2)
    integer, intent(in) :: nex
    ! Specific to the double lamp post model:
    real(wp) :: tau_d, tau_p, phase_d, phase_p
    complex(wp) :: cexp_d, cexp_p
    ! Various loop variables:
    integer :: i, j, m
    real(wp) :: f, E, fac
    complex(wp) :: W0, W1, W2, W3, Stemp, cexp_phi

    ! Zero outputs
    arrays%ReSraw = 0.0_wp
    arrays%ImSraw = 0.0_wp

    if (model_args%boost .lt. 0 .and. config%DC .eq. 1) then
        ! Reflection only:
        do j = 1,config%nf
           do i = 1,nex
              do m=1,model_args%nlp
                  if (m .gt. 1) then
                      arrays%ReW0(m,i,j) = model_args%eta *                    &
                          arrays%ReW0(m,i,j)
                  end if
                  arrays%ReSraw(i,j) = arrays%ReSraw(i,j) +                    &
                      (-model_args%boost) * arrays%ReW0(m,i,j)
              end do
            end do
        end do
        ! And we're done!
        return
    end if

    ! Calculating both the reflection and continuum components.

    ! Set up extra terms if second lamp post is presented:
    if( model_args%nlp > 1 ) then
        ! A sanity / safety check:
        if (model_args%nlp .ne. 2) then
            print *, "Assertion failed: nlp must be 2, but is", model_args%nlp
            error stop 1
        end if
        arrays%ReW0(2,:,:) = model_args%eta * arrays%ReW0(2,:,:)
        arrays%ImW0(2,:,:) = model_args%eta * arrays%ImW0(2,:,:)
        arrays%ReW1(2,:,:) = model_args%eta * arrays%ReW1(2,:,:)
        arrays%ImW1(2,:,:) = model_args%eta * arrays%ImW1(2,:,:)
        arrays%ReW2(2,:,:) = model_args%eta * arrays%ReW2(2,:,:)
        arrays%ImW2(2,:,:) = model_args%eta * arrays%ImW2(2,:,:)
        arrays%ReW3(2,:,:) = model_args%eta * arrays%ReW3(2,:,:)
        arrays%ImW3(2,:,:) = model_args%eta * arrays%ImW3(2,:,:)
        tau_d = tauso(2) - tauso(1)
        tau_p = (model_args%h(2) - model_args%h(1))/(model_args%beta_p)
    end if

    do j=1,config%nf
        if (config%DC == 1) then
            f = 0.0_wp
        else
            f = config%flo *                                                   &
                (config%fhi/config%flo)**((real(j, wp)-0.5_wp) / real(config%nf, wp))
        endif
        do i=1,nex
            do m=1,model_args%nlp
                E = 0.5_wp * ( arrays%earx(i) + arrays%earx(i-1) )
                fac = log(gso(m)/((1.0_wp+model_args%zcos)*E))

                ! Set up phase factors
                if (model_args%nlp > 1) then
                    phase_d = 2.0_wp * pi * tau_d * f
                    phase_p = 2.0_wp * pi * tau_p * f
                    cexp_d = cmplx(cos(phase_d),sin(phase_d), kind=wp)
                    cexp_p = cmplx(cos(phase_p),sin(phase_p), kind=wp)
                else
                    cexp_d = cmplx(1.0_wp, 0.0_wp, kind=wp)
                    cexp_p = cmplx(1.0_wp, 0.0_wp, kind=wp)
                end if

                ! Set up transfer functions
                W0 = model_args%boost *                                        &
                    cmplx(arrays%ReW0(m,i,j),arrays%ImW0(m,i,j), kind=wp)
                W1 = (1-config%DC) * model_args%boost *                        &
                    cmplx(arrays%ReW1(m,i,j),arrays%ImW1(m,i,j), kind=wp)
                W2 = (1-config%DC) * model_args%boost *                        &
                    cmplx(arrays%ReW2(m,i,j),arrays%ImW2(m,i,j), kind=wp)
                W3 = config%ionvar * (1-config%DC) * model_args%boost *        &
                    cmplx(arrays%ReW3(m,i,j),arrays%ImW3(m,i,j), kind=wp)

                cexp_phi =                                                     &
                    cmplx(cos(model_args%DelAB(m)),sin(model_args%DelAB(m)), kind=wp)

                ! Calculate complex covariance
                ! Note: the reason we use complex here is to ease the
                ! calculations when we add all the extra phases from the double
                ! lamp post
                Stemp = model_args%g(m) * cexp_phi *                     &
                    (W1 + W2 + fac * cexp_d * arrays%contx(i,m))
                Stemp = Stemp + W0 + W3 + cexp_d*arrays%contx(i,m)
                Stemp = cexp_p * Stemp
                ! Separate into real/imaginary parts for compatibility with the
                ! rest of the code
                arrays%ReSraw(i,j) = arrays%ReSraw(i,j) + real(Stemp, wp)
                arrays%ImSraw(i,j) = arrays%ImSraw(i,j) + aimag(Stemp)
            end do
         end do
    end do
end subroutine sum_continuum_reflection_transfer_functions
