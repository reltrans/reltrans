!-----------------------------------------------------------------------
subroutine kernel_check(Cp, param20, param21, param20save, param21save, fhi, flo, fhisave, flosave, nf, nfsave, needtrans)
!!! Checks if reltrans needs to calculate the kernel!!!
!!! Arg:
  !   Cp: variable to set which reltrans 
  !   param20: parameter array for reltrans, reltransCp, reltransD
  !   param21: parameter array for reltransDCp
  !   param20save: saved array
  !   param21save: saved array
  !   fhi: high frequency range
  !   flo: low frequency range
  !   fhisave: saved frequency 
  !   flosave: saved frequency 
  !   nf: number of frequency bins
  !   nfsave: saved number
  !   (output) needtrans: logical variable to check the kernel calculation 
  implicit none 
  integer         , intent(in)  :: nf, nfsave, Cp
  real            , intent(in)  :: param20(20), param21(21)
  real            , intent(in)  :: param20save(20), param21save(21)
  double precision, intent(in)  :: fhi, flo, fhisave, flosave
  logical         , intent(out) :: needtrans
  integer :: i

  if( .not. needtrans )then
! loop to check the first 8 parameters which are the same for all reltrans
     do i = 1,8
        if( abs( param20(i) - param20save(i) ) .gt. 1e-7 ) then
           needtrans = .true.
           goto 666
        end if           
        if( abs( param21(i) - param21save(i) ) .gt. 1e-7 ) then
           needtrans = .true.
           goto 666
        end if
     end do
     if (Cp .eq. 1 .or. Cp .eq. 2) then ! if reltransD or reltransDCp
        if( abs( param20(10) - param20save(10) ) .gt. 1e-7 ) then 
           needtrans = .true. ! density (lognep) in reltransD
           goto 666
        end if
        if( abs( param21(10) - param21save(10) ) .gt. 1e-7 ) then 
           needtrans = .true.  ! density (lognep) in reltransDCp
           goto 666
        end if
     endif

666  continue 
     
!check if frequency range and frequency grid have changed 
     if     ( abs( nf  - nfsave  ) .gt. 1e-7) then
        needtrans = .true.
     else if( abs( fhi - fhisave ) .gt. 1e-7) then
        needtrans = .true.
     else if( abs( flo - flosave ) .gt. 1e-7) then
        needtrans = .true.
     end if

  end if
end subroutine kernel_check
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine conv_check(Cp, Cpsave, param20, param21, param20save, param21save, needconv)
!!! Checks if reltrans needs to re-convolve !!!
!!! Arg:
  !   Cp: variable to set which reltrans 
  !   param20: parameter array for reltrans, reltransCp, reltransD
  !   param21: parameter array for reltransDCp
  !   param20save: saved array
  !   param21save: saved array
  !   (output) needtrans: logical variable to check the kernel calculation 
  implicit none
  integer, intent(in)  :: Cp, Cpsave
  real   , intent(in)  :: param20(19), param21(20)
  real   , intent(in)  :: param20save(19), param21save(20)
  logical, intent(out) :: needconv

  if ( abs(Cp - Cpsave) .gt. 1e-7  ) then
     needconv = .true.
  else if ( abs( param20(9) - param20save(9) ) .gt. 1e-7 ) then
     needconv = .true. ! Afe iron abundance 
  else if ( abs( param21(9) - param21save(9) ) .gt. 1e-7 ) then
     needconv = .true. ! Afe iron abundance 
  else if ( abs( param20(10) - param20save(10) ) .gt. 1e-7 ) then
     needconv = .true. ! Either kTe or Ecut (in reltransCp or reltrans)
!     NOTE in the case of reltransD param20 is lognep which has been already checked before, so no problem 
  else if ( abs( param21(11) - param21save(11) ) .gt. 1e-7 ) then
     needconv = .true. ! kTe for reltransDCp 
  end if


end subroutine conv_check
!-----------------------------------------------------------------------
