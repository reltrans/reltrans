!-----------------------------------------------------------------------
function sysfref(rin,rlp,cosd,ndelta,cosdout)
! Calculates the relxill definition of reflection fraction        
! In: rin,rlp,ndelta,cosd,cosdout
! Out: stsfref
  implicit none
  integer ndelta
  double precision sysfref,rin,rlp(ndelta),cosd(ndelta),cosdout
  integer n
  double precision cosdin
  n = 1
  do while( rin .gt. rlp(n) )
    n = n + 1
  end do
  !Inperpolate
  if( n .eq. 1 )then
    !Need to extrapolate
    cosdin = (cosd(n+1)-cosd(n))*(rin-rlp(n))/(rlp(n+1)-rlp(n))
    cosdin = cosdin + cosd(n)
  else
    !Can actually interpolate
    cosdin = (cosd(n)-cosd(n-1))*(rin-rlp(n-1))/(rlp(n)-rlp(n-1))
    cosdin = cosdin + cosd(n-1)
  end if
  sysfref = ( cosdin - cosdout ) / ( 1.0 + cosdout )
  return
end function sysfref  
!-----------------------------------------------------------------------
