!-----------------------------------------------------------------------
subroutine get_reflionx(ear, ne, param, ifl, photar)
  implicit none
  integer, intent(in)  :: ne, ifl
  real,    intent(in)  :: ear(0:ne), param(7)
  real,    intent(out) :: photar(ne)
  real                 :: photer(ne)
  character (len=500)  :: filenm,strenv
  character (len=200)  :: envnm
  logical              :: needfile
  data needfile/.true./
  save needfile
  
! Get the reflionx table  
  if( needfile )then
     envnm  = 'REFLIONX_FILE'
     filenm = strenv(envnm)
     if( trim(filenm) .eq. 'none' )then
        write(*,*)"Enter reflionx file (with full path)"
        read(*,'(a)')filenm
     end if
     needfile = .false.
  end if
  
! Interpolate a spectrum from it
  call xsatbl(ear, ne, param, filenm, ifl, photar, photer) 

  return
end subroutine get_reflionx
!-----------------------------------------------------------------------

