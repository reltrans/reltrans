subroutine set_param(Cp,dset,param,nlp,h,a,inc,rin,rout,zcos,Gamma,logxi,Dkpc,Afe,lognep,Cutoff_s, Cutoff_obs, &
                     eta_0,eta,beta_p,Nh,boost,qboost,Mass,honr,b1,b2,floHz,fhiHz,ReIm,DelA,DelAB,&
                     g,Anorm,resp,refvar,verbose)
!!! Sets the parameters of reltrans depending on the Cp variable
  use dyn_gr
  use isco
  implicit none
  integer         , intent(in)   :: Cp, dset, nlp, verbose
  real            , intent(in)   :: param(32)
  double precision, intent(out)  :: h(nlp), a, inc, rin, rout, zcos, Gamma
  double precision, intent(out)  :: honr, b1, b2, qboost, eta_0, eta
  real            , intent(out)  :: logxi, Afe, lognep, Cutoff_s, Cutoff_obs 
  real            , intent(out)  :: Nh, boost, Mass, floHz, fhiHz, beta_p
  real            , intent(out)  :: DelA, DelAB(nlp), g(nlp), Anorm, Dkpc
  integer         , intent(out)  :: ReIm, resp, refvar
  integer m
  double precision               :: disco

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
  if (Cp .gt. -0.5) then 
     Cutoff_s   = param(12) !In the case of the nthcomp continuum the kTe is in the corona frame temperature (reflionx falls into this)
  else if (Cp .eq. -1) then
     Cutoff_obs = param(12) !In the case of the powerlaw continuum Ecut is in the observer restframe
  endif
  eta_0    = param(13)
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

  if( dset .eq. 1 )then
     Dkpc = param(9)
  else
     Dkpc = 0.0
  end if

  !Check the inner radius and the heigh of the source 
    !Set minimum r  and convert rin and h to rg
    rh     = 1.d0+sqrt(1.d0-a**2)
    ! write(*,*) 'parameters'
    if( abs(a) .gt. 0.999 ) a = sign(a,1.d0) * 0.999
    ! rmin   = disco( a )
    risco = disco(a)
    call isco_KE_AM(a)
    if( rin .lt. 0.d0 ) rin = abs(rin) * risco
    rmin = rh + 0.0001
    if( rin .lt. rmin )then
        write(*,*)"Warning! rin < Event Horizon! Set to the Event Horizon (+0.0001Rg)"
        rin = rmin
    end if
    do m=1,nlp 
        if( h(m) .lt. 0.d0 ) h(m) = abs(h(m)) * rh

    if( verbose .gt. 0 ) then
       write(*,*) "isco  = ",risco
       write(*,*) "Event horizon = ", rh
       write(*,*) "Minimum radius alloweded = ", rmin
       write(*,*) "h (Rg) = ",h(m)
    endif

        if( h(m) .lt. 1.5d0*rh )then
            write(*,*)"Warning! h<1.5*rh! Set to 1.5*rh"
            h(m) = 1.5d0 * rh
        end if 
    end do
    if( verbose .gt. 0 ) write(*,*)"rin (Rg)=",rin

  return
end subroutine set_param
!-----------------------------------------------------------------------
