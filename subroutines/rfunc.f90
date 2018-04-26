!-----------------------------------------------------------------------
      function rfunc(a,mu0)
! Sets minimum rn to use for impact parameter grid depending on mu0
! This is just an analytic function based on empirical calculations:
! I simply set a=0.998, went through the full range of mu0, and then
! calculated the lowest rn value for which there was a disk crossing.
! The function used here makes sure the calculated rnmin is always
! slightly lower than the one required.
      implicit none
      double precision rfunc,mu0,a
      if( a .gt. 0.8 )then
        rfunc = 1.5d0 + 0.5d0 * mu0**5.5d0
        rfunc = min( rfunc , -0.1d0 + 5.6d0*mu0 )
        rfunc = max( 0.1d0 , rfunc )
      else
        rfunc = 3.0d0 + 0.5d0 * mu0**5.5d0
        rfunc = min( rfunc , -0.2d0 + 10.0d0*mu0 )
        rfunc = max( 0.1d0 , rfunc )
      end if
      end
!-----------------------------------------------------------------------
