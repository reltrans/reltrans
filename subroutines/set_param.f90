subroutine set_param(h, a, inc, rin, rout, zcos, Gamma, logxi, Afe, &
     lognep, kTe, Ecut_obs, Nh, afac, Mass, floHz, fhiHz, ReIm,&
     DelA, DelAB, g, par_reltrans20, par_reltrans21, Cp, resp)
!!! Sets the parameters of reltrans depending on the Cp variable
  implicit none
  double precision, intent(out)  :: h, a, inc, rin, rout, zcos, Gamma
  real            , intent(out)  :: logxi, Afe, lognep, kTe, Ecut_obs
  real            , intent(out)  :: Nh, afac, Mass, floHz, fhiHz
  real            , intent(out)  :: DelA, DelAB, g
  integer         , intent(out)  :: ReIm, resp
  integer         , intent(in)   :: Cp
  real            , intent(inout):: par_reltrans20(20), par_reltrans21(21)


  if( Cp .eq. -1 )then       ! reltrans
     h        = dble( par_reltrans20(1) )
     a        = dble( par_reltrans20(2) )
     inc      = dble( par_reltrans20(3) )
     rin      = dble( par_reltrans20(4) )
     rout     = dble( par_reltrans20(5) )
     zcos     = dble( par_reltrans20(6) )
     Gamma    = dble( par_reltrans20(7) )
     logxi    =       par_reltrans20(8)
     Afe      =       par_reltrans20(9)
     Ecut_obs =       par_reltrans20(10)
     Nh       =       par_reltrans20(11)
     afac     =       par_reltrans20(12)
     Mass     = dble( par_reltrans20(13) )
     floHz    =       par_reltrans20(14)
     fhiHz    =       par_reltrans20(15)
     ReIm     =  int( par_reltrans20(16) )
     DelA     =       par_reltrans20(17)
     DelAB    =       par_reltrans20(18)
     g        =       par_reltrans20(19)
     resp     =  int( par_reltrans20(20) )

     kTe = Ecut_obs / 2.5
     par_reltrans21 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if ( Cp .eq. -2 )then ! xillverCp
     h        = dble( par_reltrans20(1) )
     a        = dble( par_reltrans20(2) )
     inc      = dble( par_reltrans20(3) )
     rin      = dble( par_reltrans20(4) )
     rout     = dble( par_reltrans20(5) )
     zcos     = dble( par_reltrans20(6) )
     Gamma    = dble( par_reltrans20(7) )
     logxi    =       par_reltrans20(8)
     Afe      =       par_reltrans20(9)
     kTe      =       par_reltrans20(10)
     Nh       =       par_reltrans20(11)
     afac     =       par_reltrans20(12)
     Mass     = dble( par_reltrans20(13) )
     floHz    =       par_reltrans20(14)
     fhiHz    =       par_reltrans20(15)
     ReIm     =  int( par_reltrans20(16) )
     DelA     =       par_reltrans20(17)
     DelAB    =       par_reltrans20(18)
     g        =       par_reltrans20(19)
     resp     =  int( par_reltrans20(20) )

     Ecut_obs = 2.5 * kTe 
     par_reltrans21 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if( Cp .eq. 1 )then   ! reltransD
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
     Nh       =       par_reltrans20(11)
     afac     =       par_reltrans20(12)
     Mass     = dble( par_reltrans20(13) )
     floHz    =       par_reltrans20(14)
     fhiHz    =       par_reltrans20(15)
     ReIm     =  int( par_reltrans20(16) )
     DelA     =       par_reltrans20(17)
     DelAB    =       par_reltrans20(18)
     g        =       par_reltrans20(19) 
     resp     =  int( par_reltrans20(20) )

     Ecut_obs = 300.0    ! In this model the Ecut is fixed to 300
     kTe = Ecut_obs / 2.5
      par_reltrans21 = 0.0 !set paramters of reltransDCp to zero to avoid mistake in the check subroutine

  else if ( Cp .eq. 2 )then  ! reltransDCp
     h        = dble( par_reltrans21(1) )
     a        = dble( par_reltrans21(2) )
     inc      = dble( par_reltrans21(3) )
     rin      = dble( par_reltrans21(4) )
     rout     = dble( par_reltrans21(5) )
     zcos     = dble( par_reltrans21(6) )
     Gamma    = dble( par_reltrans21(7) )
     logxi    =       par_reltrans21(8)
     Afe      =       par_reltrans21(9)
     lognep   =       par_reltrans21(10)
     kTe      =       par_reltrans21(11)
     Nh       =       par_reltrans21(12)
     afac     =       par_reltrans21(13)
     Mass     = dble( par_reltrans21(14) )
     floHz    =       par_reltrans21(15)
     fhiHz    =       par_reltrans21(16)
     ReIm     =  int( par_reltrans21(17) )
     DelA     =       par_reltrans21(18)
     DelAB    =       par_reltrans21(19)
     g        =       par_reltrans21(20)
     resp     =  int( par_reltrans21(21) )
  
     Ecut_obs = 2.5 * kTe !work out the energy cut off for the powerlaw starting form kTe because we are using xillverDCp but also cut-off powerlaw for the pivoting effects
     par_reltrans20 = 0.0 !set paramters of reltrans, reltransD, reltransCp to zero to avoid mistake in the check subroutine

  else
     write(*,*) 'No reltrans  model available for this configuration'
     stop 
  end if
  return
end subroutine set_param
!-----------------------------------------------------------------------
