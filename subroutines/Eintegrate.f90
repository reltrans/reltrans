!-----------------------------------------------------------------------
function Eintegrate(Elo,Ehi,nex,earx,photarx,dlogE)
! Integrates specific flux from Elo to Ehi
! Photarx is specific photon flux, therefore integral is
! sum_i E(i) * photarx(i)
  implicit none
  integer nex
  real Eintegrate,Elo,Ehi,earx(0:nex),photarx(nex),Emin,dlogE
  integer ilo,ihi,i
  real E
! Calculate integration bounds
  Emin = earx(0)
  ilo  = ceiling( log10(Elo/Emin) / dlogE )
  ilo  = min(ilo,nex)
  ilo  = max(ilo,1)
  ihi  =   floor( log10(Ehi/Emin) / dlogE )
  ihi  = min(ihi,nex)
  ihi  = max(ihi,1)
  ihi  = max(ihi,ilo)
! Do the integration
  Eintegrate = 0.0
  do i = ilo,ihi
     E = 0.5 * ( earx(i) + earx(i-1) )
     Eintegrate = Eintegrate + E * photarx(i)
  end do
  return
end function Eintegrate
!-----------------------------------------------------------------------
