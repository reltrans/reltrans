subroutine set_param(dset,aset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut,&
                     eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                     g,Anorm,resp,Oab,Feab,za,nh_xiab, logxi_xiab, fcov_xiab, z_xiab, refvar,verbose)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  integer         , intent(in)   :: dset, aset, nlp, verbose
  real            , intent(in)   :: param(39)
  double precision, intent(out)  :: h(nlp), a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost, eta_0, eta
  real            , intent(out)  :: logxi, Afe, lognep, Ecut, beta_p
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB(nlp), g(nlp), Anorm, Dkpc
  real            , intent(out)  :: Oab, Feab, za
  real            , intent(out)  :: nh_xiab, logxi_xiab, fcov_xiab, z_xiab
  integer         , intent(out)  :: ReIm, resp, refvar
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
  lognep   = param(11) 
  Ecut     = param(12) !This is the corona frame temperature for the double LP model, and the observed one otherwise
  !if(nlp .gt. 1 .and. param(13) .eq. 0) then
  !    eta_0 = 1.e-4
  !    if (verbose .gt. 0) print*,"WARNING: low eta_0 for double LP, careful with pivoting!!!"
  !else
  eta_0    = param(13)
  !end if  
  eta      = param(14)
  beta_p   = param(15)
  Nh       = param(16)
  boost    = param(17)
  qboost   = dble( param(18) )
  Mass     = dble( param(19) )
  honr     = dble( param(20) )
  b1       = dble( param(21) )
  b2       = dble( param(22) )
  floHz    = param(23)
  fhiHz    = param(24)
  ReIm     = int( param(25) )
  DelA     = param(26)
  do m=1,nlp 
    DelAB(m) = param(27+(m-1)*nlp) 
    g(m)     = param(28+(m-1)*nlp)   
  end do
  Anorm    = param(31)
  resp     = param(32)
  Oab      = param(33)
  Feab     = param(34)
  za     = param(35)
  nh_xiab = param(36)
  logxi_xiab = param(37)
  fcov_xiab = param(38)
  z_xiab   = param(39)

  if( dset .eq. 1 )then
     Dkpc = param(9)
  else
     Dkpc = 0.0
  end if
  
  !WIP optimization, need to think how best to handle this for cases when phiab/g change
  !but we don't want to re-do all the convolutions
 ! if(all(g .eq. 0 ) .or. all(DelAB .eq. 0)) then
 !    refvar = 0  
 ! end if

  return
end subroutine set_param
!-----------------------------------------------------------------------
