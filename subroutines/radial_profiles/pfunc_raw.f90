!-----------------------------------------------------------------------
function pfunc_raw(mu,b1,b2,boost)
  use rtconstants, only: wp
  implicit none
  real(wp) pfunc_raw,mu,b1,b2,boost
  real(wp) calB,norm,pm,mup,p
  calB = 1.0_wp/boost
  if( mu .le. 0.0_wp )then
     norm = 1.0_wp
  else
     norm = calB
  end if
  pm = mu/sign(mu,1.0_wp)
  if( mu .eq. 0.0_wp ) pm = 1.0_wp
  mup = pm * ( norm**2 * (1.0_wp/mu**2-1.0_wp) + 1.0_wp )**(-0.5_wp)
  p   = 1.0_wp + (b1+abs(b2))*abs(mup) + b2*mup**2
  p   = p * sqrt( 1.0_wp + mup**2 * (norm**2-1.0_wp) )
  pfunc_raw = p  
  return
end function pfunc_raw  
!-----------------------------------------------------------------------

