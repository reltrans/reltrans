!-----------------------------------------------------------------------
      function get_env_int(name,default)
        implicit none
        integer get_env_int, stat, default
        character (len=5) str
        character (len=*) name
        stat = 0        
        call get_environment_variable(trim(name),str)
        call str2int(str, get_env_int, stat)
        if( stat .ne. 0 )then
           get_env_int = default
        end if
        stat = 0
        return
      end function get_env_int
!-----------------------------------------------------------------------    
