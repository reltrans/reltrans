subroutine set_param(dset, param, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
     lognep, Ecut_obs, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, Anorm, resp, Ea1keV, Ea2keV, Eb1keV, Eb2keV, ABslope, gslope)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  integer         , intent(in)   :: dset
  real            , intent(in)   :: param(32)
  double precision, intent(out)  :: h, a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost
  real            , intent(out)  :: logxi, Afe, lognep, Ecut_obs
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB, g, Anorm, Dkpc
  integer         , intent(out)  :: ReIm, resp

  real                :: Ea1keV, Ea2keV, Eb1keV, Eb2keV
  real                :: gslope, ABslope

! Read in basic parameter array  
  h        = dble( param(1) )
  a        = dble( param(2) )
  inc      = dble( param(3) )
  rin      = dble( param(4) )
  rout     = dble( param(5) )
  zcos     = dble( param(6) )
  Gamma    = dble( param(7) )
  logxi    = param(8)  ! or distance
  Afe      = param(9)
  lognep   = param(10) ! will just be a dummy variable for Cp<0
  Ecut_obs = param(11) ! or kTe, but we never need to know the difference
  Nh       = param(12)
  boost    = param(13)
  qboost   = dble( param(14) )
  Mass     = dble( param(15) )
  honr     = dble( param(16) )
  b1       = dble( param(17) )
  b2       = dble( param(18) )
  floHz    = param(19)
  fhiHz    = param(20)
  ReIm     = int( param(21) )
  DelA     = param(22)
  DelAB    = param(23)
  g        = param(24)
  Anorm    = param(25)
  resp     = param(26)

  if( dset .eq. 1 )then
     Dkpc = param(8)
  else
     Dkpc = 0.0
  end if

     Ea1keV   = param(27)
     Ea2keV   = param(28)
     Eb1keV   = param(29)
     Eb2keV   = param(30)
     ABslope  = param(31)
     gslope   = param(32)
  
  return
end subroutine set_param
!-----------------------------------------------------------------------
