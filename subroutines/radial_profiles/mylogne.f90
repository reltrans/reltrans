!-----------------------------------------------------------------------
function mylogne_zoneA(r,rin)
! Calculates log10(ne). Don't let r = rin
  implicit none
  double precision mylogne_zoneA,r,rin
  mylogne_zoneA = 1.5 * log10(r) - 2.0 * log10( 1.0 - sqrt(rin/r) )
  return
end function mylogne_zoneA
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function mylogne_zoneB(r,rin)
! Calculates log10(ne). Don't let r = rin
  implicit none
  double precision mylogne_zoneB,r,rin
  mylogne_zoneB =  log10(r**(-33./20.)) ! + 2./5. * log10( 1.0 - sqrt(rin/r) )
  return
end function mylogne_zoneB
!-----------------------------------------------------------------------
