!-----------------------------------------------------------------------
! Error handling module
!-----------------------------------------------------------------------
module error_handler
    implicit none
    
    logical, save :: error_occurred = .false.
    character(len=256), save :: error_message = ""
    
contains

    subroutine log_error(message)
        character(len=*), intent(in) :: message
        error_occurred = .true.
        error_message = trim(message)
        write(*,'(A)') "ERROR: " // trim(message)
    end subroutine log_error
    
    subroutine log_warning(message)
        character(len=*), intent(in) :: message
        write(*,'(A)') "WARNING: " // trim(message)
    end subroutine log_warning
    
    subroutine reset_error()
        error_occurred = .false.
        error_message = ""
    end subroutine reset_error
    
    function has_error() result(err)
        logical :: err
        err = error_occurred
    end function has_error
    
    subroutine return_zero_spectrum(photar, ne)
        integer, intent(in) :: ne
        real, intent(out) :: photar(ne)
        photar = 0.0
        write(*,'(A)') "Returning zero spectrum due to error"
    end subroutine return_zero_spectrum

end module error_handler
