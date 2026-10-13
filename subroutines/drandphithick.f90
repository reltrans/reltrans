!-----------------------------------------------------------------------
      subroutine drandphithick(alpha,beta,cosi,costheta,r,phi)
!
! A disk with an arbitrary thickness
! The angle between the normal to the midplane and the disk surface is theta
! The inclination angle is i
      use rtconstants, only: wp
      implicit none
      real(wp) alpha,beta,cosi,sini,r,phi
      real(wp) pi,costheta,sintheta,x,a,b,c,det
      real(wp) mu,sinphi
!      double precision muplus,muminus,ra,rb,rab,xplus1,xminus1,xplus2,xminus2
      pi = acos(-1.0_wp)
      sintheta = sqrt( 1.0_wp - costheta**2 )
      sini     = sqrt( 1.0_wp - cosi**2 )
      x        = alpha / beta
      if( abs(alpha) .lt. abs(tiny(alpha)) .and. abs(beta) .lt. abs(tiny(beta))  )then
        mu = 0.0_wp
        r  = 0.0_wp
      else if( abs(beta) .lt. abs(tiny(beta)) )then
        mu     = sini*costheta/(cosi*sintheta)
        sinphi = sign( 1.0_wp , alpha ) * sqrt( 1.0_wp - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      else if( abs(alpha) .lt. abs(tiny(alpha)) )then
        mu     = 1.0_wp
        sinphi = 0.0_wp
        r      = beta / ( sini*costheta - cosi*sintheta )
      else
        a      = sintheta**2 + x**2*cosi**2*sintheta**2
        b      = -2*x**2*sini*cosi*sintheta*costheta
        c      = x**2*sini**2*costheta**2-sintheta**2
        det    = b**2 - 4.0_wp * a * c
        if( det .lt. 0.0_wp )then
           write(*,*)"Error in drandphithick"
           write(*,*)"determinant <0!!!"
           write(*,*)"costheta=",costheta
           write(*,*)"h/r=",costheta/sintheta
           write(*,*)"cosi=",cosi
           write(*,*)"alpha,beta=",alpha,beta
        end if
        if( beta .gt. 0.0_wp )then
          mu     = ( -b + sqrt( det ) ) / ( 2.0_wp * a )
        else
          mu     = ( -b - sqrt( det ) ) / ( 2.0_wp * a )
        end if
        sinphi = sign( 1.0_wp , alpha ) * sqrt( 1.0_wp - mu**2 )
        r      = alpha / ( sintheta * sinphi )
      end if
      phi = atan2( sinphi , mu )
      return
      end subroutine drandphithick
!-----------------------------------------------------------------------
