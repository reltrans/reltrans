

!-----------------------------------------------------------------------
      function dlgfac(a,mu0,alpha,r)
!c Calculates g-factor for a disk in the BH equatorial plane
      implicit none
      double precision dlgfac,a,mu0,alpha,r
      double precision sin0,omega,Delta,Sigma2,gtt,gtp,gpp
      sin0   = sqrt( 1.0 - mu0**2 )
      omega  = 1. / (r**1.5+a)
      Delta  = r**2 - 2*r + a**2
      Sigma2 = (r**2+a**2)**2 - a**2 * Delta
      gtt    = 4*a**2/Sigma2 - r**2*Delta/Sigma2
      gtp    = -2*a/r
      gpp    = Sigma2/r**2
      dlgfac = sqrt( -gtt - 2*omega*gtp - omega**2.*gpp )
      dlgfac = dlgfac / ( 1.+omega*alpha*sin0 )
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function areafac(r,a)
! Calculates dA/dr, where A is the surface area of a disk ring
! *ADJUSTED FOR ORBITAL MOTION* BY MULTIPLYING BY THE LORENTZ FACTOR
      implicit none
      real areafac,r,a
      real Dm,dArbydr,pi,lorfac
      pi = acos(-1.0)
      Dm     = r**2 - 2*r + a**2
      dArbydr = ( r**4+a**2*r**2+2*a**2*r ) / Dm
      dArbydr = 2*pi*sqrt(dArbydr)
      areafac = lorfac(r,a) * dArbydr
      return
      end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function dareafac(r,a)
! Calculates dA/dr, where A is the surface area of a disk ring
! *ADJUSTED FOR ORBITAL MOTION* BY MULTIPLYING BY THE LORENTZ FACTOR
      implicit none
      double precision dareafac,r,a
      double precision Dm,dArbydr,pi,dlorfac
      pi = acos(-1.0)
      Dm     = r**2 - 2*r + a**2
      dArbydr = ( r**4+a**2*r**2+2*a**2*r ) / Dm
      dArbydr = 2*pi*sqrt(dArbydr)
      dareafac = dlorfac(r,a) * dArbydr
      return
      end
!-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
      function glpfac(r,a,h)
! Calculates blue shift expreienced by a photon travelling from
! an on-axis point source to a point on a Keplerian disk
! Works for pro- and retrograde spins.
      implicit none
      real glpfac,r,a,h
      real angvel,Dh,gphiphi
      angvel = 1.0 / ( r**1.5 + abs(a) )
      Dh     = h**2 - 2*h + a**2
      gphiphi = r**2+a**2+2*a**2/r
      glpfac = 1-2./r+4*a*angvel/r-gphiphi*angvel**2
      glpfac = Dh/(h**2+a**2) / glpfac
      glpfac = sqrt( glpfac )
      return
      end
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      function dglpfac(r,a,h)
! Calculates blue shift expreienced by a photon travelling from
! an on-axis point source to a point on a Keplerian disk
! Works for pro- and retrograde spins.
      implicit none
      double precision dglpfac,r,a,h
      double precision angvel,Dh,gphiphi
      angvel = 1.0 / ( r**1.5 + abs(a) )
      Dh     = h**2 - 2*h + a**2
      gphiphi = r**2+a**2+2*a**2/r
      dglpfac = 1-2./r+4*a*angvel/r-gphiphi*angvel**2
      dglpfac = Dh/(h**2+a**2) / dglpfac
      dglpfac = sqrt( dglpfac )
      return
      end
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      function dgsofac(a,h)
! Calculates blue shift expreienced by a photon travelling from
! an on-axis point source to a distant, stationary observer.
! Works for pro- and retrograde spins.
      implicit none
      double precision dgsofac,a,h,Dh
      Dh      = h**2 - 2*h + a**2
      dgsofac = Dh / ( h**2 + a**2 )
      dgsofac = sqrt( dgsofac )
      return
      end
!-----------------------------------------------------------------------


    
!-----------------------------------------------------------------------
      function lorfac(r,a)
! Calculates Lorentz factor for rotating disk element
      implicit none
      real lorfac,r,a
      real Delta,BigA,Omega,v
      Delta  = r**2 - 2*r + a**2
      BigA   = (r**2+a**2)**2 - a**2*Delta
      Omega  = 1.0 / ( r**1.5 + abs(a) )
      v      = ( Omega * bigA - 2*a*r ) / ( r**2 * sqrt(Delta) )
      lorfac = ( 1 - v**2 )**(-0.5)
      return
      end
!-----------------------------------------------------------------------

      
!-----------------------------------------------------------------------
      function dlorfac(r,a)
! Calculates Lorentz factor for rotating disk element
      implicit none
      double precision dlorfac,r,a
      double precision Delta,BigA,Omega,v
      Delta   = r**2 - 2*r + a**2
      BigA    = (r**2+a**2)**2 - a**2*Delta
      Omega   = 1.0 / ( r**1.5 + abs(a) )
      v       = ( Omega * bigA - 2*a*r ) / ( r**2 * sqrt(Delta) )
      dlorfac = ( 1 - v**2 )**(-0.5)
      return
      end
!-----------------------------------------------------------------------
