!-----------------------------------------------------------------------
      subroutine randphi(alpha,beta,cos0,r,phi)
      use rtconstants, only: wp
      implicit none
      real(wp) alpha,beta,cos0,r,phi
      real(wp) tanphi,pi
      pi = acos(-1.0_wp)
      tanphi = -cos0*alpha/beta
      if( abs(tanphi) .ge. HUGE(tanphi) )then
        if( alpha .gt. 0.0_wp )then
          phi =  0.5_wp * pi
        else
          phi = -0.5_wp * pi
        end if
      else if( abs(tanphi) .lt. TINY(tanphi) )then
        if( beta .lt. 0.0_wp )then
          phi = 0.0_wp
        else
          phi = pi
        end if
      else if( alpha .gt. 0.0_wp .and. beta .lt. 0.0_wp )then
        phi = atan( tanphi )         
        do while( phi .lt. 0.0_wp )
          phi = phi + pi
        end do
        do while( phi .gt. 0.5_wp*pi )
          phi = phi - pi      
        end do
      else if( alpha .gt. 0.0_wp .and. beta .gt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. 0.5_wp*pi )
          phi = phi + pi
        end do
        do while( phi .gt. pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0_wp .and. beta .gt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. pi )
          phi = phi + pi
        end do
        do while( phi .gt. 1.5_wp*pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0_wp .and. beta .lt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. 1.5_wp*pi )
          phi = phi + pi
        end do
        do while( phi .gt. 2.0_wp*pi )
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
      use rtconstants, only: wp
      implicit none
      real(wp) alpha,beta,cos0,r,phi
      real(wp) tanphi,pi
      pi = acos(-1.0_wp)
      tanphi = -cos0*alpha/beta
      if( abs(tanphi) .ge. HUGE(tanphi) )then
        if( alpha .gt. 0.0_wp )then
          phi =  0.5_wp * pi
        else
          phi = -0.5_wp * pi
        end if
      else if( abs(tanphi) .lt. TINY(tanphi) )then
        if( beta .lt. 0.0_wp )then
          phi = 0.0_wp
        else
          phi = pi
        end if
      else if( alpha .gt. 0.0_wp .and. beta .lt. 0.0_wp )then
        phi = atan( tanphi )         
        do while( phi .lt. 0.0_wp )
          phi = phi + pi
        end do
        do while( phi .gt. 0.5_wp*pi )
          phi = phi - pi      
        end do
      else if( alpha .gt. 0.0_wp .and. beta .gt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. 0.5_wp*pi )
          phi = phi + pi
        end do
        do while( phi .gt. pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0_wp .and. beta .gt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. pi )
          phi = phi + pi
        end do
        do while( phi .gt. 1.5_wp*pi )
          phi = phi - pi      
        end do
      else if( alpha .lt. 0.0_wp .and. beta .lt. 0.0_wp )then
        phi = atan( tanphi )
        do while( phi .lt. 1.5_wp*pi )
          phi = phi + pi
        end do
        do while( phi .gt. 2.0_wp*pi )
          phi = phi - pi
        end do
      end if
      r   = sqrt(alpha**2+beta**2)
      r   = r / sqrt( sin(phi)**2 + cos0**2*cos(phi)**2 )
      return
      end
!-----------------------------------------------------------------------
