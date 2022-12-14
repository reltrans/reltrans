!-----------------------------------------------------------------------
      function myenv(name,default)
      implicit none
      integer myenv,stat,default
      character (len=5) str
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

!-----------------------------------------------------------------------
      function myenv_real(name, default)
      implicit none
      integer           :: stat 
      real              :: myenv_real, default
      character (len=5) :: str
      character (len=*) :: name
      stat = 0
      CALL get_environment_variable(trim(name),str)
      call mystr2real(str, myenv_real, stat)
      if( stat .ne. 0 )then
        myenv_real = default
      end if
      stat = 0
      return
    end function myenv_real
!-----------------------------------------------------------------------    


!-----------------------------------------------------------------------
      function strenv(name)
! STATUS is -1 if VALUE is present but too short for the environment variable;
! it is 1 if the environment variable does not exist and 2 if the processor does
! not support environment variables; in all other cases STATUS is zero.        
      implicit none
      character (len=500) strenv
      character (len=200) name
      integer length,status
      status = 0
      !CALL get_environment_variable(trim(name),strenv)
      CALL GET_ENVIRONMENT_VARIABLE(trim(name),strenv,LENGTH,STATUS)
      if( status .ne. 0 .or. length .eq. 0 )then
         strenv = 'none'         
      end if
      return
      end function strenv
!-----------------------------------------------------------------------    
