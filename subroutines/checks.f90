!-----------------------------------------------------------------------
subroutine need_check(Cp,Cpsave,param,paramsave,fhi,flo,fhisave,flosave,nf,nfsave,needtrans,needconv)
!
! Checks if reltrans needs to calculate the kernel
!
! Parameters that rtrans() is sensitive to:
! (1-8):   h,a,inc,rin,rout,zcos,Gamma,logxi/Dkpc
! (10):    lognep
! (14):    qboost
! (16-18): honr,b1,b2
! (20):    Anorm
! Also need to check if the frequency range changes
!
! Parameters the restframe spec is sensitive to
! (9):     Afe
! (10):    Ecut/kTe  
  
!!! Arg:
  ! INPUTS
  !   Cp:        defines which model
  !   Cpsave:    saved Cp
  !   param:     parameter array
  !   paramsave: saved array
  !   fhi:       high frequency range
  !   flo:       low frequency range
  !   fhisave:   saved frequency 
  !   flosave:   saved frequency 
  !   nf:        number of frequency bins
  !   nfsave:    saved number
  ! OUTPUTS
  !   needtrans: if true, we must do the kernel calculation
  !   neecconv:  if true, we must do the convolution
  implicit none 
  integer         , intent(in)  :: Cp, Cpsave, nf, nfsave
  real            , intent(in)  :: param(32), paramsave(32)
  real            , parameter   :: tol = 1e-7
  double precision, intent(in)  :: fhi, flo, fhisave, flosave
  double precision, parameter   :: dtol = 1e-10
  logical         , intent(out) :: needtrans,needconv
  integer :: i
  needtrans = .false.
  needconv  = .false.
! First check the parameter entries
  do i = 1,8
     if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  end do
  i = 10
  if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  do i = 14,18
     if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  end do
  i = 25
  if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
! Now check if frequency range and frequency grid have changed 
  if( nf .ne. nfsave ) then
     needtrans = .true.
  else if( abs( fhi - fhisave ) .gt. dtol) then
     needtrans = .true.
  else if( abs( flo - flosave ) .gt. dtol) then
     needtrans = .true.
  end if
!Now for needconv
  if( needtrans ) needconv = .true.
  if( Cp .ne. Cpsave ) needconv = .true.
  if( abs( param(9) - paramsave(9) ) .gt. tol ) needconv = .true.
  if( abs( param(10) - paramsave(10) ) .gt. tol ) needconv = .true.
end subroutine need_check
!-----------------------------------------------------------------------


! !-----------------------------------------------------------------------
! subroutine conv_check(Cp, Cpsave, param19, param20, param19save, param20save, needconv)
! !!! Checks if reltrans needs to re-convolve !!!
! !!! Arg:
!   !   Cp: variable to set which reltrans 
!   !   param19: parameter array for reltrans, reltransCp, reltransD
!   !   param20: parameter array for reltransDCp
!   !   param19save: saved array
!   !   param20save: saved array
!   !   (output) needtrans: logical variable to check the kernel calculation 
!   implicit none
!   integer, intent(in)  :: Cp, Cpsave
!   real   , intent(in)  :: param19(19), param20(20)
!   real   , intent(in)  :: param19save(19), param20save(20)
!   logical, intent(out) :: needconv

!   if ( abs(Cp - Cpsave) .gt. 1e-7  ) then
!      needconv = .true.
!   else if ( abs( param19(9) - param19save(9) ) .gt. 1e-7 ) then
!      needconv = .true. ! Afe iron abundance 
!   else if ( abs( param20(9) - param20save(9) ) .gt. 1e-7 ) then
!      needconv = .true. ! Afe iron abundance 
!   else if ( abs( param19(10) - param19save(10) ) .gt. 1e-7 ) then
!      needconv = .true. ! Either kTe or Ecut (in reltransCp or reltrans)
! !     NOTE in the case of reltransD param19 is lognep which has been already checked before, so no problem 
!   else if ( abs( param20(11) - param20save(11) ) .gt. 1e-7 ) then
!      needconv = .true. ! kTe for reltransDCp 
!   end if


! end subroutine conv_check
! !-----------------------------------------------------------------------
