  ! FFTW-based convolution of the reflection kernels with the restframe
  ! spectrum. Lives outside conv_mod (as a free subroutine) so that it can
  ! `use common_types`: conv_mod is a dependency of common_types, so a
  ! procedure inside conv_mod cannot itself use the derived types.
  subroutine conv_one_FFTw(config, model_args, dyn, earx, Gamma, photarx,      &
       reline, imline, ReW_conv, ImW_conv)
    use conv_mod
    use common_types
    implicit none
    type(t_config),          intent(in)    :: config
    type(t_model_arguments), intent(in)    :: model_args
    real                :: dyn
    real, intent(in)    :: photarx(nex), earx(0:nex), Gamma
    real, intent(in)    :: reline(model_args%nlp,nex), imline(model_args%nlp,nex)
    real, intent(inout) :: ReW_conv(model_args%nlp,nex), ImW_conv(model_args%nlp,nex)
    complex :: conv(nec),padFT_photarx(nec)
    complex :: padFT_reline(nec),  padFT_imline(nec)
    integer :: m, i
    real    :: photmax, depad_conv(nex), E

    do m=1,model_args%nlp
       if (config%DC .eq. 1 ) then
          ! call padding4FT(photarx,padFT_photarx)
          call padding4FT_xillver(photarx,padFT_photarx)

          call padding4FT(reline(m,:),padFT_reline)

          conv = (padFT_photarx * padFT_reline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do

       else
          call padding4FT(photarx       , padFT_photarx)
          call padding4FT(reline(m,:),padFT_reline)
          call padding4FT(imline(m,:),padFT_imline)

          conv = (padFT_photarx * padFT_reline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do

          conv = (padFT_photarx * padFT_imline) * nexm1
          call de_paddingFT(dyn, conv, depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ImW_conv(m,i) = ImW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do

       endif
    end do

  end subroutine conv_one_FFTw
