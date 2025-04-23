!-----------------------------------------------------------------------
    function get_env_char(name, default)
      implicit none
      character (len=*), intent(in) :: name, default

      integer :: length, stat
      character (len=200) get_env_char
      logical :: trim_name
      stat = 0
      call get_environment_variable(trim(name), get_env_char, length, stat, trim_name)
      if( stat .eq. 1 )then
         get_env_char = default
         write(*,'(A,A,A)') 'You did not set ', trim(name), ' environmanet variable'
         write(*,'(A)') ' the code assumes the tables are in this folder (./)'
      else if (stat .eq. 2) then
         get_env_char = default
         write(*,'(A)') 'The processor does not support environment variables'
         write(*,'(A)') ' the code assumes the tables are in this folder (./)'
      else if (stat .eq. -1) then 
         get_env_char = default
         write(*,'(A)') 'The path is too long for this code contact the developers'            
         write(*,'(A)') ' the code assumes the tables are in this folder (./)'
      end if
      return
    end function get_env_char
!-----------------------------------------------------------------------    
