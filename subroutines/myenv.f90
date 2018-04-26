!-----------------------------------------------------------------------
      function myenv(name,default)
      implicit none
      integer myenv,stat,default
      character (len=3) str
      character (len=*) name
      stat = 0        
      CALL get_environment_variable(trim(name),str)
      call mystr2int(str,myenv,stat)
      if( stat .ne. 0 )then
        myenv = default
      end if
      stat = 0
      return
      end function myenv
!-----------------------------------------------------------------------    
