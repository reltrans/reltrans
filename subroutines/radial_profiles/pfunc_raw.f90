!-----------------------------------------------------------------------
function pfunc_raw(mu,b1,b2,boost)
  implicit none
  double precision pfunc_raw,mu,b1,b2,boost
  double precision calB,norm,pm,mup,p
  calB = 1.0/boost
  if( mu .le. 0.0 )then
     norm = 1.0
  else
     norm = calB
  end if
  pm = mu/sign(mu,1.d0)
  if( mu .eq. 0.d0 ) pm = 1.0
  mup = pm * ( norm**2 * (1.0/mu**2-1.0) + 1.0 )**(-0.5)
  p   = 1.0 + (b1+abs(b2))*abs(mup) + b2*mup**2
  p   = p * sqrt( 1.0 + mup**2 * (norm**2-1.0) )
  pfunc_raw = p  
  return
end function pfunc_raw  
!-----------------------------------------------------------------------

