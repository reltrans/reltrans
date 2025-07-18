  subroutine conv_one_FFT(dyn,photarx,reline,imline,ReW_conv,ImW_conv,DC,nlp)
    use conv_mod
    implicit none
    integer, intent(in) :: DC, nlp 
    real                :: dyn
    real, intent(in)    :: photarx(nex)
    real, intent(in)    :: reline(nlp,nex), imline(nlp,nex)
    real, intent(inout) :: ReW_conv(nlp,nex), ImW_conv(nlp,nex)
    complex :: FTphotarx(nex_conv), FTreline(nex_conv), FTimline(nex_conv)
    complex :: FTreconv(4*nex),FTimconv(4*nex)
    integer :: m, i
    real    :: photmax, depad_conv(nex)
    ! real, parameter :: nexm1 = 1. / real(nex_conv)
    
    do m=1,nlp  
        if (DC .eq. 1 ) then
           call pad4FFT(nex, photarx, FTphotarx)
           call pad4FFT(nex, reline(m,:), FTreline)
           FTreconv = (FTreline * FTphotarx) !* nexm1
           call pad4invFFT(dyn,nex,FTreconv,depad_conv)
 
           ReW_conv(m,:) = ReW_conv(m,:) + depad_conv

        else
           call pad4FFT(nex,photarx, FTphotarx)
           call pad4FFT(nex,reline(m,:),FTreline)
           call pad4FFT(nex,imline(m,:),FTimline)
           FTreconv = (FTreline * FTphotarx) !* nexm1
           FTimconv = (FTimline * FTphotarx) !* nexm1
           call pad4invFFT(dyn,nex,FTreconv,depad_conv)
           ReW_conv(m,:) = ReW_conv(m,:) + depad_conv
           call pad4invFFT(dyn,nex,FTimconv,depad_conv) 
           ImW_conv(m,:) = ImW_conv(m,:) + depad_conv 

           ! call FTcnv(nex, reline(m,:), photarx, ReW_conv(m,:))
           ! call FTcnv(nex, imline(m,:), photarx, ImW_conv(m,:))

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
