
!-----------------------------------------------------------------------
subroutine fold(nex, earx, photarx, spec)
! Initmatrix must have alreadt been called
! Input: photarx(1:nex); i.e. (dN/dE)*dE
! Output: spec(1:numchn); in count rate vs channel number
  use rtconstants, only: wp
  use telematrix
  implicit none
  integer nex,i,j,k
  real(wp) earx(0:nex),photarx(nex),spec(numchn)
  real(wp) Si(nenerg),E,dE,E2Sx(nex)
  
  !Convert to E^2*dN/dE for better accuracy
  do i = 1,nex
     E  = 0.5_wp * ( earx(i) + earx(i-1) )
     dE = earx(i) - earx(i-1)
     E2Sx(i)   = E**2 * photarx(i) / dE      
  end do
  
  !Rebin input arrays onto internal telescope energy grid
  call rebinE(earx,E2Sx,nex,En,Si,nenerg)
  
  !Convert back to (dN/dE)*dE
  do i = 1,nenerg
     E  = 0.5_wp * ( En(i) + En(i-1) )
     dE = En(i) - En(i-1)
     Si(i) = Si(i) / E**2 * dE
  end do

  !Fold around response
  spec = 0.0_wp
  do J = 1, NENERG
     do K = 1,NGRP(J)
        do I = FCHAN(J,K) + 1, LCHAN(J,K)
           spec(I) = spec(I) + Si(J) * RESP(I,J)
        end do
     end do
  end do
  
  return
end subroutine fold
!-----------------------------------------------------------------------
