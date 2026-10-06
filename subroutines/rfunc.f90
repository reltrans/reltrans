!-----------------------------------------------------------------------
      function rfunc(a,mu0)
! Sets minimum rn to use for impact parameter grid depending on mu0
! This is just an analytic function based on empirical calculations:
! I simply set a=0.998, went through the full range of mu0, and then
! calculated the lowest rn value for which there was a disk crossing.
! The function used here makes sure the calculated rnmin is always
! slightly lower than the one required.
      use rtconstants, only: wp
      implicit none
      real(wp) rfunc,mu0,a
      if( a .gt. 0.8_wp )then
        rfunc = 1.5_wp + 0.5_wp * mu0**5.5_wp
        rfunc = min( rfunc , -0.1_wp + 5.6_wp*mu0 )
        rfunc = max( 0.1_wp , rfunc )
      else
        rfunc = 3.0_wp + 0.5_wp * mu0**5.5_wp
        rfunc = min( rfunc , -0.2_wp + 10.0_wp*mu0 )
        rfunc = max( 0.1_wp , rfunc )
      end if
      end
!-----------------------------------------------------------------------
