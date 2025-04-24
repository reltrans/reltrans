!-----------------------------------------------------------------------
function logxiraw(re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts,mudisk,gsd)
! In: re,spin,h,honr,rlp,dcosdr,ndelta,rmin,npts
! Out: logxiraw,gsd
  implicit none
  integer ndelta,npts,kk,get_index
  double precision re,spin,h,honr,rlp(ndelta),dcosdr(ndelta),rmin,gsd
  double precision cosfac,interper,logxiraw,dareafac,dglpfacthick,newtex
  double precision mudisk
  !Calculate source to disc blueshift at this radius
  gsd = dglpfacthick(re,spin,h,mudisk)
  !gsd = dglpfac(re,spin,h)
  !Find the rlp bin that corresponds to re
  kk = get_index(rlp,ndelta,re,rmin,npts)
  !Interpolate to get |d\cos\delta/dr| at r=re
  cosfac = interper(rlp,dcosdr,ndelta,re,kk)
  !Extrapolate to Newtonian if needs be
  if( kk .eq. npts ) cosfac = newtex(rlp,dcosdr,ndelta,re,h,honr,kk)
  !Now can do the calculation
  logxiraw = log10( gsd**2 * cosfac / dareafac(re,spin) )
  return
end function logxiraw  
!-----------------------------------------------------------------------
