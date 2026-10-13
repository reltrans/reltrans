!-----------------------------------------------------------------------
subroutine xilimits(nex,earx,nlp,contx,DeltaGamma,gso,lens,z,dlogxi1,dlogxi2)
! Inputs: nex,earx,contx,DeltaGamma,gso
! Outputs: dlogxi1,dlogxi2
! logxi1 = logxi0 + dlogxi1
    use rtconstants, only: wp
    implicit none
    integer nex,i,m,nlp
    real(wp) earx(0:nex),contx(nex,nlp),contx_sum(nex),DeltaGamma,gso(nlp),lens(nlp),z,dlogxi1,dlogxi2
    real(wp) num1,num2,den,logS1,logS2,gsoz,E,gso_avg

    !before calculating the differential of the ionisation, set up the arrays properly depending on the number of lampposts
    !note: if we have multiple LPs we have to a) calculate an effective gso factor by averaging over the lensing factor (ie, how
    !luminous each LP appears given their height) and b) get a total continuum flux by summing over each LPs array.
    !for a single lamp posts there is no need to do all this stuff, but we need to read in the contx/gso factors properly to be
    !able to call the following code in the same way
    if(nlp .eq. 1 ) then
        contx_sum = contx(:,1)
        gso_avg = gso(1)
    else
        contx_sum = 0.0_wp
        gso_avg = 0.0_wp
        do m=1,nlp
            contx_sum = contx_sum + contx(:,m)
            gso_avg = gso_avg + lens(m)*gso(m)
        end do
        gso_avg = gso_avg/sum(lens)
    end if

    gsoz = gso_avg / (1.0_wp+z)
    num1 = 0.0_wp
    num2 = 0.0_wp
    den = 0.0_wp
    do i = 1,nex
       E   = 0.5_wp * ( earx(i) + earx(i-1) )
       if (E .ge. 0.1_wp .and. E .le. 1e3_wp) then
          num2 = num2 + E**(1.0_wp-0.5_wp*DeltaGamma) * contx_sum(i)
          num1 = num1 + E**(1.0_wp+0.5_wp*DeltaGamma) * contx_sum(i)
          den  = den  + E * contx_sum(i)
       endif
     end do
    logS1 = log10(num1/den)
    logS2 = log10(num2/den)
    dlogxi1 = -0.5_wp*DeltaGamma * log10(gsoz) + logS1
    dlogxi2 =  0.5_wp*DeltaGamma * log10(gsoz) + logS2
    return
end subroutine xilimits
!-----------------------------------------------------------------------  
      
