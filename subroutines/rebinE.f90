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


!-----------------------------------------------------------------------
subroutine rebinE(earx,px,nex,ear,p,ne)
!General rebinning scheme, should be nice and robust - BUT IT FUCKING ISN'T
!i,nex,earx,px = input
!j,ne,ear,p    = output
  implicit none
  integer i,nex,j,ne,ilo,ihi
  real earx(0:nex),ear(0:ne),px(nex),p(ne),Ehigh,Elow,upper,lower
  real FRAC,Ej,Ei,pi,Ei2,pi2,grad,cons,Ehi,Elo,phi,plo
  logical interp
  ilo = 1
  do j = 1,ne
     do while( earx(ilo) .le. ear(j-1) .and. ilo .lt. nex )
        ilo = ilo + 1
     end do
     ihi = ilo
     do while( earx(ihi) .le. ear(j) .and. ihi .lt. nex )
        ihi = ihi + 1
     end do
     if( ihi .gt. ilo )then
        p(j) = 0.0
        do i = ilo,ihi
           lower = MAX( earx(i-1) , ear(j-1)  )
           upper = MIN( earx(i)   , ear(j)    )
           p(j) = p(j) + px(i) * ( upper - lower )
        end do
        p(j) = p(j) / ( ear(j) - ear(j-1) )
     else
        !Interpolate (or extrapolate)
        i = ilo
        Ei = 0.5 * ( earx(i) + earx(i-1) )
        if( Ei .gt. Ej ) i = ilo - 1
        i = max( i , 2     )
        i = min( i , nex-1 )
        Ej  = 0.5 * ( ear(j) + ear(j-1) )
        Ehi = 0.5 * ( earx(i+1) + earx(i)   )
        Elo = 0.5 * ( earx(i)   + earx(i-1) )
        phi = px(i+1)
        plo = px(i)
        p(j) = plo + (phi-plo)*(Ej-Elo)/(Ehi-Elo)
     end if
     if( ilo .gt. 1 ) ilo = ilo - 1
  end do
  RETURN   
END subroutine rebinE
!-----------------------------------------------------------------------


