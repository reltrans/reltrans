subroutine set_param(dset, param, nlp, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
     lognep, Ecut, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, Anorm, resp)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  integer         , intent(in)   :: dset, nlp
  real            , intent(in)   :: param(27)
  double precision, intent(out)  :: h(nlp), a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost
  real            , intent(out)  :: logxi, Afe, lognep, Ecut
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB, g, Anorm, Dkpc
  integer         , intent(out)  :: ReIm, resp

  
! Read in basic parameter array  
  if (nlp .eq. 1) then
    h(1)   = dble( param(1) )
  else 
    h(1)   = dble( param(1) )
    h(2)   = dble( param(2) )
  endif  
  a        = dble( param(3) )
  inc      = dble( param(4) )
  rin      = dble( param(5) )
  rout     = dble( param(6) )
  zcos     = dble( param(7) )
  Gamma    = dble( param(8) )
  logxi    = param(9)  ! or distance
  Afe      = param(10)
  lognep   = param(11) ! will just be a dummy variable for Cp<0
  Ecut     = param(12) ! or kTe, but we never need to know the difference. NOTE: this is the corona frame temperature/cutoff for the double LP model, and the observed one otherwise
  Nh       = param(13)
  boost    = param(14)
  qboost   = dble( param(15) )
  Mass     = dble( param(16) )
  honr     = dble( param(17) )
  b1       = dble( param(18) )
  b2       = dble( param(19) )
  floHz    = param(20)
  fhiHz    = param(21)
  ReIm     = int( param(22) )
  DelA     = param(23)
  DelAB    = param(24)
  g        = param(25)
  Anorm    = param(26)
  resp     = param(27)

  if( dset .eq. 1 )then
     Dkpc = param(9)
  else
     Dkpc = 0.0
  end if

  return
end subroutine set_param
!-----------------------------------------------------------------------
