!-----------------------------------------------------------------------
function zA_logne(r,rin,lognep)
! log(ne), where ne is the density.
! This function is normalised to have a maximum of lognep.
! We can therefore calculate the ionization parameter by taking
! 4 \pi * Fx / nemin * zA_one_on_ne
  implicit none
  double precision zA_logne,r,rin,lognep,rp
  rp       = 25./9. * rin
  zA_logne = lognep + 1.5*log10(r/rp)
  zA_logne = zA_logne + 2.0*log10( 1.0 - sqrt( rin / rp ) )
  zA_logne = zA_logne - 2.0*log10( 1.0 - sqrt( rin / r  ) )
  return
end function zA_logne
!-----------------------------------------------------------------------
