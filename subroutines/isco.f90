!-----------------------------------------------------------------------
      function dISCO(a)
      !ISCO in Rg 
      use rtconstants, only: wp
      implicit none
      real(wp) a,dISCO,z1,z2
      z1 = ( 1.0_wp - a**2.0_wp )**(1.0_wp/3.0_wp)
      z1 = z1 * ( (1.0_wp+a)**(1.0_wp/3.0_wp)+(1.0_wp-a)**(1.0_wp/3.0_wp))+1.0_wp
      z2 = sqrt( 3.0_wp * a**2.0_wp + z1**2.0_wp )
      if(a.ge.0.0_wp)then
        dISCO = 3.0_wp + z2 - sqrt( (3.0_wp-z1) * (3.0_wp + z1 + 2.0_wp*z2) )
      else
        dISCO = 3.0_wp + z2 + sqrt( (3.0_wp-z1) * (3.0_wp + z1 + 2.0_wp*z2) )
      end if
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function ISCO(a)
      !ISCO in Rg 
      use rtconstants, only: wp
      implicit none
      real(wp) a,ISCO,z1,z2
      z1 = ( 1.0_wp - a**2.0_wp )**(1.0_wp/3.0_wp)
      z1 = z1 * ( (1.0_wp+a)**(1.0_wp/3.0_wp)+(1.0_wp-a)**(1.0_wp/3.0_wp))+1.0_wp
      z2 = sqrt( 3.0_wp * a**2.0_wp + z1**2.0_wp )
      if(a.ge.0.0_wp)then
        ISCO = 3.0_wp + z2 - sqrt( (3.0_wp-z1) * (3.0_wp + z1 + 2.0_wp*z2) )
      else
        ISCO = 3.0_wp + z2 + sqrt( (3.0_wp-z1) * (3.0_wp + z1 + 2.0_wp*z2) )
      end if
      return
      end
!-----------------------------------------------------------------------
