!-----------------------------------------------------------------------
function zA_logne(r,rin,lognep)
! log(ne), where ne is the density.
! This function is normalised to have a maximum of lognep.
! We can therefore calculate the ionization parameter by taking
! 4 \pi * Fx / nemin * zA_one_on_ne
  use rtconstants, only: wp
  implicit none
  real(wp) zA_logne,r,rin,lognep,rp
  rp       = 25.0_wp/9.0_wp * rin
  zA_logne = lognep + 1.5_wp*log10(r/rp)
  zA_logne = zA_logne + 2.0_wp*log10( 1.0_wp - sqrt( rin / rp ) )
  zA_logne = zA_logne - 2.0_wp*log10( 1.0_wp - sqrt( rin / r  ) )
  return
end function zA_logne
!-----------------------------------------------------------------------
