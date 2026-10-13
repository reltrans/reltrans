!-----------------------------------------------------------------------
subroutine get_reflionx(ear, ne, param, ifl, photar)
  use rtconstants, only: wp
  use xspec_interface, only: table_model
  use xillver_tables
  implicit none
  integer, intent(in)  :: ne, ifl
  real(wp),    intent(in)  :: ear(0:ne), param(7)
  real(wp),    intent(out) :: photar(ne)
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
     path_name_reflionx_table = trim(filenm)
     write(*,*) 'Set the reflionx table at ', path_name_reflionx_table
     needfile = .false.
  end if
! Interpolate a spectrum from it
! Pass the null terminator for C compatability, as xsatbl is an external C function.
  call table_model(ear, ne, param, [trim(path_name_reflionx_table), char(0)],  &
    ifl, photar)
  return
end subroutine get_reflionx
!-----------------------------------------------------------------------

