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
