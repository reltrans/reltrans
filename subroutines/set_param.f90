subroutine set_param(h, a, inc, rin, rout, zcos, Gamma, logxi, Afe, &
     lognep, kTe, Ecut_obs, Nh, afac, Mass, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, par_reltrans19, par_reltrans20, Cp)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  double precision, intent(out)  :: h, a, inc, rin, rout, zcos, Gamma
  real            , intent(out)  :: logxi, Afe, lognep, kTe, Ecut_obs
  real            , intent(out)  :: Nh, afac, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB, g
  integer         , intent(out)  :: ReIm
  integer         , intent(in)   :: Cp
  real            , intent(inout):: par_reltrans19(19), par_reltrans20(20)


  if( Cp .eq. -1 )then       ! reltrans
     h        = dble( par_reltrans19(1) )
     a        = dble( par_reltrans19(2) )
     inc      = dble( par_reltrans19(3) )
     rin      = dble( par_reltrans19(4) )
     rout     = dble( par_reltrans19(5) )
     zcos     = dble( par_reltrans19(6) )
     Gamma    = dble( par_reltrans19(7) )
     logxi    =       par_reltrans19(8)
     Afe      =       par_reltrans19(9)
     Ecut_obs =       par_reltrans19(10)
     Nh       =       par_reltrans19(11)
     afac     =       par_reltrans19(12)
     Mass     = dble( par_reltrans19(13) )
     floHz    =       par_reltrans19(14)
     fhiHz    =       par_reltrans19(15)
     ReIm     =  int( par_reltrans19(16) )
     DelA     =       par_reltrans19(17)
     DelAB    =       par_reltrans19(18)
     g        =       par_reltrans19(19)

     kTe = Ecut_obs / 2.5
     par_reltrans20 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if ( Cp .eq. -2 )then ! xillverCp
     h        = dble( par_reltrans19(1) )
     a        = dble( par_reltrans19(2) )
     inc      = dble( par_reltrans19(3) )
     rin      = dble( par_reltrans19(4) )
     rout     = dble( par_reltrans19(5) )
     zcos     = dble( par_reltrans19(6) )
     Gamma    = dble( par_reltrans19(7) )
     logxi    =       par_reltrans19(8)
     Afe      =       par_reltrans19(9)
     kTe      =       par_reltrans19(10)
     Nh       =       par_reltrans19(11)
     afac     =       par_reltrans19(12)
     Mass     = dble( par_reltrans19(13) )
     floHz    =       par_reltrans19(14)
     fhiHz    =       par_reltrans19(15)
     ReIm     =  int( par_reltrans19(16) )
     DelA     =       par_reltrans19(17)
     DelAB    =       par_reltrans19(18)
     g        =       par_reltrans19(19)

     Ecut_obs = 2.5 * kTe 
     par_reltrans20 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if( Cp .eq. 1 )then   ! reltransD
     h        = dble( par_reltrans19(1) )
     a        = dble( par_reltrans19(2) )
     inc      = dble( par_reltrans19(3) )
     rin      = dble( par_reltrans19(4) )
     rout     = dble( par_reltrans19(5) )
     zcos     = dble( par_reltrans19(6) )
     Gamma    = dble( par_reltrans19(7) )
     logxi    =       par_reltrans19(8)
     Afe      =       par_reltrans19(9)
     lognep   =       par_reltrans19(10)
     Nh       =       par_reltrans19(11)
     afac     =       par_reltrans19(12)
     Mass     = dble( par_reltrans19(13) )
     floHz    =       par_reltrans19(14)
     fhiHz    =       par_reltrans19(15)
     ReIm     =  int( par_reltrans19(16) )
     DelA     =       par_reltrans19(17)
     DelAB    =       par_reltrans19(18)
     g        =       par_reltrans19(19) 

     Ecut_obs = 300.0    ! In this model the Ecut is fixed to 300
     kTe = Ecut_obs / 2.5
      par_reltrans20 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if ( Cp .eq. 2 )then  ! reltransDCp
     h        = dble( par_reltrans20(1) )
     a        = dble( par_reltrans20(2) )
     inc      = dble( par_reltrans20(3) )
     rin      = dble( par_reltrans20(4) )
     rout     = dble( par_reltrans20(5) )
     zcos     = dble( par_reltrans20(6) )
     Gamma    = dble( par_reltrans20(7) )
     logxi    =       par_reltrans20(8)
     Afe      =       par_reltrans20(9)
     lognep   =       par_reltrans20(10)
     kTe      =       par_reltrans20(11)
     Nh       =       par_reltrans20(12)
     afac     =       par_reltrans20(13)
     Mass     = dble( par_reltrans20(14) )
     floHz    =       par_reltrans20(15)
     fhiHz    =       par_reltrans20(16)
     ReIm     =  int( par_reltrans20(17) )
     DelA     =       par_reltrans20(18)
     DelAB    =       par_reltrans20(19)
     g        =       par_reltrans20(20)
  
     Ecut_obs = 2.5 * kTe !work out the energy cut off for the powerlaw starting form kTe because we are using xillverDCp but also cut-off powerlaw for the pivoting effects
     par_reltrans19 = 0.0 !set paramters of reltrans, reltransD, reltransCp to zero to avoid mistake in the check subroutine

  else
     write(*,*) 'No reltrans  model available for this configuration'
     stop 
  end if
  return
end subroutine set_param
!-----------------------------------------------------------------------
