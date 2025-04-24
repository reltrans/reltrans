!-----------------------------------------------------------------------
function mylogne(r,rin)
! Calculates log10(ne). Don't let r = rin
  implicit none
  double precision mylogne,r,rin
  mylogne = 1.5 * log10(r) - 2.0 * log10( 1.0 - sqrt(rin/r) )
  return
end function mylogne 
!-----------------------------------------------------------------------
