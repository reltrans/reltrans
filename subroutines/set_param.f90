subroutine set_param(dset, param, nlp, h, a, inc, rin, rout, zcos, Gamma, logxi, Dkpc, Afe, &
     lognep, Ecut, eta_0, eta, Nh, boost, qboost, Mass, honr, b1, b2, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, Anorm, resp)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  integer         , intent(in)   :: dset, nlp
  real            , intent(in)   :: param(31)
  double precision, intent(out)  :: h(nlp), a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost, eta_0, eta
  real            , intent(out)  :: logxi, Afe, lognep, Ecut
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB(nlp), g(nlp), Anorm, Dkpc
  integer         , intent(out)  :: ReIm, resp
  integer m

  !TBD: DelAB, g also arryas of size nlp 
! Read in basic parameter array   
  do m=1,nlp 
    h(m) = dble(param(m))
  end do 
  a        = dble( param(3) )
  inc      = dble( param(4) )
  rin      = dble( param(5) )
  rout     = dble( param(6) )
  zcos     = dble( param(7) )
  Gamma    = dble( param(8) )
  logxi    = param(9)  ! or distance
  Afe      = param(10)
  lognep   = param(11) ! will just be a dummy variable for Cp<0
  Ecut     = param(12) ! or kTe, but we never need to know the difference. NOTE: this is the corona frame temperature for the double LP model, and the observed one otherwise
  eta_0    = param(13)
  eta      = param(14)
  Nh       = param(15)
  boost    = param(16)
  qboost   = dble( param(17) )
  Mass     = dble( param(18) )
  honr     = dble( param(19) )
  b1       = dble( param(20) )
  b2       = dble( param(21) )
  floHz    = param(22)
  fhiHz    = param(23)
  ReIm     = int( param(24) )
  DelA     = param(25)
  do m=1,nlp 
    DelAB(m) = param(26+(m-1)*nlp) 
    g(m)     = param(27+(m-1)*nlp)   
  end do
  Anorm    = param(30)
  resp     = param(31)

  if( dset .eq. 1 )then
     Dkpc = param(9)
  else
     Dkpc = 0.0
  end if

  return
end subroutine set_param
!-----------------------------------------------------------------------
