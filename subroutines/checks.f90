!-----------------------------------------------------------------------
subroutine need_check(Cp,Cpsave,param,paramsave,fhi,flo,fhisave,flosave,nf,nfsave,needtrans,needconv)
!
! Checks if reltrans needs to calculate the kernel
!
! Parameters that rtrans() is sensitive to:
! (1-9):   h1,h2,a,inc,rin,rout,zcos,Gamma,logxi/Dkpc
! (11):    lognep
! (15):    qboost
! (17-19): honr,b1,b2
! (21):    Anorm
! Also need to check if the frequency range changes
!
! Parameters the restframe spec is sensitive to
! (10):     Afe
! (12):    Ecut/kTe  
  
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
  real            , intent(in)  :: param(26), paramsave(26)
  real            , parameter   :: tol = 1e-7
  double precision, intent(in)  :: fhi, flo, fhisave, flosave
  double precision, parameter   :: dtol = 1e-10
  logical         , intent(out) :: needtrans,needconv
  integer :: i
  needtrans = .false.
  needconv  = .false.
! First check the parameter entries
  do i = 1,9
     if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  end do
  i = 11
  if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  do i = 15,19
     if( abs( param(i) - paramsave(i) ) .gt. tol ) needtrans = .true.
  end do
  i = 26
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
  if( abs( param(10) - paramsave(10) ) .gt. tol ) needconv = .true.
  if( abs( param(12) - paramsave(12) ) .gt. tol ) needconv = .true.
end subroutine need_check
!-----------------------------------------------------------------------
