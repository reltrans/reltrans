subroutine set_param(dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Ecut,&
                     eta_0,eta,beta_p,Nh,boost,qboost,t_diff_sec,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                     g,Anorm,resp,refvar,verbose)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  integer         , intent(in)   :: dset, nlp, verbose
  real            , intent(in)   :: param(33)
  double precision, intent(out)  :: h(nlp), a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost, eta_0, eta, t_diff_sec
  real            , intent(out)  :: logxi, Afe, lognep, Ecut, beta_p
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB(nlp), g(nlp), Anorm, Dkpc
  integer         , intent(out)  :: ReIm, resp, refvar
  integer m
  !relativistic parameters and limit on rin and h
  double precision :: rmin, rh , disco

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
  t_diff_sec = param(19)
  Mass     = dble( param(20) )
  honr     = dble( param(21) )
  b1       = dble( param(22) )
  b2       = dble( param(23) )
  floHz    = param(24)
  fhiHz    = param(25)
  ReIm     = int( param(26) )
  DelA     = param(27)
  do m=1,nlp 
    DelAB(m) = param(28+(m-1)*nlp) 
    g(m)     = param(29+(m-1)*nlp)   
  end do
  Anorm    = param(32)
  resp     = param(33)

  if( dset .eq. 1 )then
     Dkpc = param(9)
  else
     Dkpc = 0.0
  end if

!tbd: avoid REim=-5/-6 
  
  !Set minimum r (ISCO) and convert rin and h to rg
  if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
  rmin   = disco( a )
  if( rin .lt. 0.d0 ) rin = abs(rin) * rmin
  rh     = 1.d0+sqrt(1.d0-a**2)
  if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin
  if( rin .lt. rmin )then
     write(*,*)"Warning! rin<ISCO! Set to ISCO"
        rin = rmin
  end if
  do m=1,nlp 
     if( h(m) .lt. 0.d0 ) h(m) = abs(h(m)) * rh
     if( verbose .gt. 0 ) write(*,*)"h (Rg)=",h(m)
     if( h(m) .lt. 1.5d0*rh )then
        write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
        h(m) = 1.5d0 * rh
     end if 
  end do
  
  !WIP optimization, need to think how best to handle this for cases when phiab/g change
  !but we don't want to re-do all the convolutions
 ! if(all(g .eq. 0 ) .or. all(DelAB .eq. 0)) then
 !    refvar = 0  
 ! end if

  return
end subroutine set_param
!-----------------------------------------------------------------------
