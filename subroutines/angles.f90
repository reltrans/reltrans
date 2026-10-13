!-----------------------------------------------------------------------
      function dinang(a,r,h,mus)
! mus = \cos\delta
! Calculates cosine of angle of incidence for a photon hitting the
! disk mid-plane from a lamppost source
! mui = pi_alpha n^alpha / ( pi_mu u^mu )
!     = pi_alpha n^alpha / ( -u^t )
!     = ( ql / r ) / ( u^t )
        use rtconstants, only: wp
        implicit none
      real(wp) dinang,a,r,h,mus
      real(wp) DeltaL,Delta,ql,Omega
      !Define useful combinations
      DeltaL  = h**2 - 2*h + a**2
      Delta   = r**2 - 2*r + a**2
      !Carter's constant of incoming photon      
      ql      = (1.0_wp-mus**2) * (h**2+a**2)**2/DeltaL - a**2
      ql      = sqrt( ql )
      !Orbital angular velocity of disk material
      Omega   = sign(a,1.0_wp) * 1.0_wp / ( r**1.5_wp + abs(a) )
      ! Calculate mui = (ql/r)/ut
      dinang = 1-2.0_wp/r + 4*a*Omega/r
      dinang = dinang - ( (r**2+a**2)**2 - a**2*Delta )/r**2 * Omega**2
      dinang = sqrt( dinang ) * ql / r
      return
      end
!-----------------------------------------------------------------------
      

!-----------------------------------------------------------------------
      function demang(a,mu0,r,alpha,beta)
! Calculates emission angle for the disk mid-plane
      use rtconstants, only: wp
      implicit none
      real(wp) demang,a,mu0,r,alpha,beta,mue
      real(wp) dlgfacthick
      !mue = dlgfac( a,mu0,alpha,r ) / r
      mue = dlgfacthick( a,mu0,alpha,r,0.0_wp ) / r
      mue = mue * sqrt( beta**2 + mu0**2 * (alpha**2-a**2) )
      demang = min( mue , 1.0_wp )
! g = - (1+z)^-1 / [ pe_\mu U^\mu ]
! \mue = - pe_\alpha n^\alpha / [ pe_\mu U^\mu ], therefore
! \mue = pe_\alpha n^\alpha g|_{z=0}
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function arg( z )
      use rtconstants, only: wp
      implicit none
      complex(wp) z
      real(wp) arg
      arg = atan2( aimag(z) , real(z, wp) )
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function darg( z )
      use rtconstants, only: wp
      implicit none
      complex(wp) z
      real(wp) darg
      darg = atan2( aimag(z) , real(z, wp) )
      return
      end
!-----------------------------------------------------------------------
