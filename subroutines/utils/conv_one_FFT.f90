  subroutine conv_one_FFT(dyn,earx,Gamma,photarx,reline,imline,ReW_conv,ImW_conv,DC,nlp)
    use rtconstants, only: wp
    use conv_mod
    implicit none
    integer, intent(in) :: DC, nlp 
    real(wp)            :: dyn
    real(wp), intent(in)    :: photarx(nex),earx(0:nex),Gamma
    real(wp), intent(in)    :: reline(nlp,nex), imline(nlp,nex)
    real(wp), intent(inout) :: ReW_conv(nlp,nex), ImW_conv(nlp,nex)
    complex(wp) :: FTphotarx(nex_conv), FTreline(nex_conv), FTimline(nex_conv)
    complex(wp) :: FTreconv(4*nex),FTimconv(4*nex)
    integer :: m, i
    real(wp) :: depad_conv(nex), E
    
    do m=1,nlp  
       if (DC .eq. 1 ) then
          call pad4FFT(nex, photarx, FTphotarx)
          call pad4FFT(nex, reline(m,:), FTreline)
          FTreconv = (FTreline * FTphotarx) !* nexm1
          call pad4invFFT(dyn,nex,FTreconv,depad_conv)

          do i = 1,nex
             E             = 0.5_wp * ( earx(i) + earx(i-1) )
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
             E             = 0.5_wp * ( earx(i) + earx(i-1) )
             ReW_conv(m,i) = ReW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
           
          call pad4invFFT(dyn,nex,FTimconv,depad_conv) 

          do i = 1,nex
             E             = 0.5_wp * ( earx(i) + earx(i-1) )
             ImW_conv(m,i) = ImW_conv(m,i) + depad_conv(i) * E**(1-Gamma)
          end do
          
        endif
    end do

  end subroutine conv_one_FFT