!-----------------------------------------------------------------------
subroutine crebin(nex,earx,ReGx,ImGx,ne,ear,ReG,ImG)
! Re-bins ReGx and ImGx from the earx(0:nex) energy array onto the
! ear(0:ne) energy array. Output is ReG and ImG.
! Gx and G are in units of (dN/dE)*dE (i.e. photar)
  implicit none
  integer nex,ne
  real earx(0:nex),ReGx(nex),ImGx(nex),ear(0:ne),ReG(ne),ImG(ne)
  integer i
  real E,dE,E2ReGx(nex),E2ImGx(nex)

  !Convert to E^2*dN/dE for better accuracy
  do i = 1,nex
     E  = 0.5 * ( earx(i) + earx(i-1) )
     dE = earx(i) - earx(i-1)
     E2ReGx(i) = E**2 * ReGx(i) / dE
     E2ImGx(i) = E**2 * ImGx(i) / dE
  end do
  
  !Re-bin
  call rebinE(earx,E2ReGx,nex,ear,ReG,ne)
  call rebinE(earx,E2ImGx,nex,ear,ImG,ne)

  !Convert back to (dN/dE)*dE
  do i = 1,ne
     E  = 0.5 * ( ear(i) + ear(i-1) )
     dE = ear(i) - ear(i-1)
     ReG(i) = ReG(i) / E**2 * dE
     ImG(i) = ImG(i) / E**2 * dE
  end do
  
  return
end subroutine crebin
!-----------------------------------------------------------------------
