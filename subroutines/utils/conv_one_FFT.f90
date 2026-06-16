  subroutine conv_one_FFT(config, model_args, dyn, earx, Gamma, photarx,       &
       reline, imline, ReW_conv, ImW_conv)
    use conv_mod
    use common_types
    implicit none
    type(t_config),          intent(in)    :: config
    type(t_model_arguments), intent(in)    :: model_args
    real                :: dyn
    real, intent(in)    :: photarx(nex),earx(nex),Gamma
    real, intent(in)    :: reline(model_args%nlp,nex), imline(model_args%nlp,nex)
    real, intent(inout) :: ReW_conv(model_args%nlp,nex), ImW_conv(model_args%nlp,nex)
    complex :: FTphotarx(nex_conv), FTreline(nex_conv), FTimline(nex_conv)
    complex :: FTreconv(4*nex),FTimconv(4*nex)
    integer :: m, i
    real    :: photmax, depad_conv(nex), E
    ! real, parameter :: nexm1 = 1. / real(nex_conv)
    
    do m=1,model_args%nlp
       if (config%DC .eq. 1 ) then
          call pad4FFT(nex, photarx, FTphotarx)
          call pad4FFT(nex, reline(m,:), FTreline)
          FTreconv = (FTreline * FTphotarx) !* nexm1
          call pad4invFFT(dyn,nex,FTreconv,depad_conv)

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
           
       else
          call pad4FFT(nex,photarx, FTphotarx)
          call pad4FFT(nex,reline(m,:),FTreline)
          call pad4FFT(nex,imline(m,:),FTimline)
          FTreconv = (FTreline * FTphotarx) !* nexm1
          FTimconv = (FTimline * FTphotarx) !* nexm1
          call pad4invFFT(dyn,nex,FTreconv,depad_conv)
           
          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
           
          call pad4invFFT(dyn,nex,FTimconv,depad_conv) 

          do i = 1,nex
             E             = 0.5 * ( earx(i) + earx(i-1) )
             ImW_conv(m,i) = ImW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
          
        endif
    end do

  end subroutine conv_one_FFT


              ! !Convolve with line profile
              ! !First FFTs
              ! call pad4FFT(nex,photarx,FTphotarx)
              ! call pad4FFT(nex,photarx_delta,FTphotarx_delta)
              ! call pad4FFT(nex,reline,FTreline)
              ! call pad4FFT(nex,imline,FTimline)
              ! call pad4FFT(nex,reline_a,FTreline_a)
              ! call pad4FFT(nex,imline_a,FTimline_a)
              ! call pad4FFT(nex,photarx_dlogxi,FTphotarx_dlogxi)
              ! !Then the multiplications and inverse FFTs
              ! FTreconv = FTreline * FTphotarx
              ! FTimconv = FTimline * FTphotarx
              ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
              ! call pad4invFFT(dyn,nex,FTimconv,imconvmu) 
              ! do i = 1,nex
              !    ReW0(i,j) = ReW0(i,j) + reconvmu(i)
              !    ImW0(i,j) = ImW0(i,j) + imconvmu(i)
              ! end do              
              ! FTreconv = FTreline_a * FTphotarx
              ! FTimconv = FTimline_a * FTphotarx
              ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
              ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
              ! do i = 1,nex
              !    ReW1(i,j) = ReW1(i,j) + reconvmu(i)
              !    ImW1(i,j) = ImW1(i,j) + imconvmu(i)
              ! end do
              ! FTreconv = FTreline * FTphotarx_delta
              ! FTimconv = FTimline * FTphotarx_delta
              ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
              ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
              ! do i = 1,nex
              !    ReW2(i,j) = ReW2(i,j) + reconvmu(i)
              !    ImW2(i,j) = ImW2(i,j) + imconvmu(i)
              ! end do
              ! FTreconv = FTreline * FTphotarx_dlogxi
              ! FTimconv = FTimline * FTphotarx_dlogxi
              ! call pad4invFFT(dyn,nex,FTreconv,reconvmu)
              ! call pad4invFFT(dyn,nex,FTimconv,imconvmu)
              ! do i = 1,nex
              !    ReW3(i,j) = ReW3(i,j) + reconvmu(i)
              !    ImW3(i,j) = ImW3(i,j) + imconvmu(i)
              ! end do
