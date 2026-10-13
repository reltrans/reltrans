!-----------------------------------------------------------------------
      function get_env_real(name, default)
        use rtconstants, only: wp
        implicit none
        integer           :: stat 
        real(wp)          :: get_env_real, default
        character (len=5) :: str
        character (len=*) :: name
        stat = 0
        CALL get_environment_variable(trim(name),str)
        call str2real(str, get_env_real, stat)
        if( stat .ne. 0 )then
           get_env_real = default
        end if
        stat = 0
        return
      end function get_env_real
!-----------------------------------------------------------------------    
