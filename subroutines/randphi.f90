!-----------------------------------------------------------------------
      subroutine randphi(alpha,beta,cos0,r,phi)
      implicit none
      real alpha,beta,cos0,r,phi
      real tanphi,pi
      pi = acos(-1.0)
      tanphi = -cos0*alpha/beta
      if( abs(tanphi) .ge. HUGE(tanphi) )then
        if( alpha .gt. 0.0 )then
          phi =  0.5 * pi
        else
          phi = -0.5 * pi
        end if
      else if( abs(tanphi) .lt. TINY(tanphi) )then
        if( beta .lt. 0.0 )then
          phi = 0.0
        else
          phi = pi
        end if
      else if( alpha .gt. 0.0 .and. beta .lt. 0.0 )then
        phi = atan( tanphi )         
        do while( phi .lt. 0.0 )
          phi = phi + pi
        end do
        do while( phi .gt. 0.5*pi )
          phi = phi - pi      
        end do
      else if( alpha .gt. 0.0 .and. beta .gt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. 0.5*pi )
          phi = phi + pi
        end do
        do while( phi .gt. pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0 .and. beta .gt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. pi )
          phi = phi + pi
        end do
        do while( phi .gt. 1.5*pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0 .and. beta .lt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. 1.5*pi )
          phi = phi + pi
        end do
        do while( phi .gt. 2.0*pi )
          phi = phi - pi
        end do
      end if
      r   = sqrt(alpha**2+beta**2)
      r   = r / sqrt( sin(phi)**2 + cos0**2*cos(phi)**2 )
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine drandphi(alpha,beta,cos0,r,phi)
      implicit none
      double precision alpha,beta,cos0,r,phi
      double precision tanphi,pi
      pi = acos(-1.0)
      tanphi = -cos0*alpha/beta
      if( abs(tanphi) .ge. HUGE(tanphi) )then
        if( alpha .gt. 0.0 )then
          phi =  0.5 * pi
        else
          phi = -0.5 * pi
        end if
      else if( abs(tanphi) .lt. TINY(tanphi) )then
        if( beta .lt. 0.0 )then
          phi = 0.0
        else
          phi = pi
        end if
      else if( alpha .gt. 0.0 .and. beta .lt. 0.0 )then
        phi = atan( tanphi )         
        do while( phi .lt. 0.0 )
          phi = phi + pi
        end do
        do while( phi .gt. 0.5*pi )
          phi = phi - pi      
        end do
      else if( alpha .gt. 0.0 .and. beta .gt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. 0.5*pi )
          phi = phi + pi
        end do
        do while( phi .gt. pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0 .and. beta .gt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. pi )
          phi = phi + pi
        end do
        do while( phi .gt. 1.5*pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0 .and. beta .lt. 0.0 )then
        phi = atan( tanphi )
        do while( phi .lt. 1.5*pi )
          phi = phi + pi
        end do
        do while( phi .gt. 2.0*pi )
          phi = phi - pi
        end do
      end if
      r   = sqrt(alpha**2+beta**2)
      r   = r / sqrt( sin(phi)**2 + cos0**2*cos(phi)**2 )
      return
      end
!-----------------------------------------------------------------------
