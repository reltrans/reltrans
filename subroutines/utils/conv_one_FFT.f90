  subroutine conv_one_FFT(dyn,photarx,reline,imline,ReW_conv,ImW_conv,DC,nlp)
    use conv_mod
    implicit none
    integer, intent(in) :: DC, nlp 
    real                :: dyn
    real, intent(in)    :: photarx(nex)
    real, intent(in)    :: reline(nlp,nex), imline(nlp,nex)
    real, intent(inout) :: ReW_conv(nlp,nex), ImW_conv(nlp,nex)
    complex :: conv(nec),padFT_photarx(nec)
    complex :: padFT_reline(nec),  padFT_imline(nec)            
    integer :: m, i
    real    :: photmax, depad_conv(nex)
    
    do m=1,nlp  
        if (DC .eq. 1 ) then   
            call FTcnv(nex, reline(m,:), photarx, ReW_conv(m,:))
            write(*,*) 'FTcnv FT'
         else
            
            call FTcnv(nex, reline(m,:), photarx, ReW_conv(m,:))
            call FTcnv(nex, imline(m,:), photarx, ImW_conv(m,:))
            write(*,*) 'FTcnv FT no DC'

        endif
    end do

    ! do i = 1, nex
    !    write(10,*) i,  reline(1,i)   
    ! enddo
    ! write(10, *) 'no no'
    ! do i = 1, nex
    !    write(11,*) i,  photarx(i)   
    ! enddo
    ! write(11, *) 'no no'
    ! do i = 1, nex
    !    write(12,*) i,  ReW_conv(1,i)             
    ! enddo
    ! write(12, *) 'no no'

  end subroutine conv_one_FFT

