

! !-----------------------------------------------------------------------
!       function dlgfac(a,mu0,alpha,r)
! !c Calculates g-factor for a disk in the BH equatorial plane
!       implicit none
!       double precision dlgfac,a,mu0,alpha,r
!       double precision sin0,omega,Delta,Sigma2,gtt,gtp,gpp
!       sin0   = sqrt( 1.0 - mu0**2 )
!       omega  = 1. / (r**1.5+a)
!       Delta  = r**2 - 2*r + a**2
!       Sigma2 = (r**2+a**2)**2 - a**2 * Delta
!       gtt    = 4*a**2/Sigma2 - r**2*Delta/Sigma2
!       gtp    = -2*a/r
!       gpp    = Sigma2/r**2
!       dlgfac = sqrt( -gtt - 2*omega*gtp - omega**2.*gpp )
!       dlgfac = dlgfac / ( 1.+omega*alpha*sin0 )
!       return
!       end
! !-----------------------------------------------------------------------

      
! !-----------------------------------------------------------------------
!       function areafac(r,a)
! ! Calculates dA/dr, where A is the surface area of a disk ring
! ! *ADJUSTED FOR ORBITAL MOTION* BY MULTIPLYING BY THE LORENTZ FACTOR
!       implicit none
!       real areafac,r,a
!       real Dm,dArbydr,pi,lorfac
!       pi = acos(-1.0)
!       Dm     = r**2 - 2*r + a**2
!       dArbydr = ( r**4+a**2*r**2+2*a**2*r ) / Dm
!       dArbydr = 2*pi*sqrt(dArbydr)
!       areafac = lorfac(r,a) * dArbydr
!       return
!       end
! !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      function dareafac(r,a)
! Calculates dA/dr, where A is the surface area of a disk ring
! *ADJUSTED FOR ORBITAL MOTION* BY MULTIPLYING BY THE LORENTZ FACTOR
! lor * sqrt(grr * g_phi_phi) * 2 pi check this 
        implicit none
      double precision               dareafac,r,a
      double precision               Dm,dArbydr,dlorfac
      double precision, parameter :: pi = acos(-1.0)
      Dm     = r**2 - 2*r + a**2
      dArbydr = ( r**4+a**2*r**2+2*a**2*r ) / Dm
      dArbydr = 2*pi*sqrt(dArbydr)
      dareafac = dlorfac(r,a) * dArbydr
      ! write(22,*) 1/dlorfac(r,a)
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


! !-----------------------------------------------------------------------
!       function dglpfac(r,a,h)
! ! Calculates blue shift expreienced by a photon travelling from
! ! an on-axis point source to a point on a Keplerian disk
! ! Works for pro- and retrograde spins.
!       implicit none
!       double precision dglpfac,r,a,h
!       double precision angvel,Dh,gphiphi
!       angvel = 1.0 / ( r**1.5 + abs(a) )
!       Dh     = h**2 - 2*h + a**2
!       gphiphi = r**2+a**2+2*a**2/r
!       dglpfac = 1-2./r+4*a*angvel/r-gphiphi*angvel**2
!       dglpfac = Dh/(h**2+a**2) / dglpfac
!       dglpfac = sqrt( dglpfac )
!       return
!       end
! !-----------------------------------------------------------------------

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
end function dlgfac
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function dlgfac_inside_isco(a, mu0, alpha, beta, r, t_r)
  !c Calculates g-factor for a disk in the BH equatorial plane
  ! imported from Andy's mummary code
  use isco 
  implicit none
  double precision dlgfac_inside_isco,a,mu0, alpha, r, beta
  double precision sin0,Delta
  double precision kr, vt_, vr_, vp_, vt, vr, vphi
  double precision eis, jis
  double precision lam, q
  integer          t_r

  Delta  = r**2. - 2.*r + a**2.0
  ! eis = (1. - 2./(3.*risco))**0.5
  ! jis = 2. * 3.**0.5 * (1 - 2.*a/(3.*risco**0.5))
  vt_ = vt(r,a) 
  vr_ = vr(r) 
  vp_ = vphi(r,a)

  sin0   = sqrt( 1.0 - mu0**2 )
  lam = -alpha * sin0

  q = (beta**2.0 - a**2.0 * mu0**2.0 + alpha**2.0 * mu0**2.0)
  
  kr=(r**4.-(q+lam**2.-a**2.)*r**2.+2.*r*(q+(lam-a)**2.) - a**2.0*q)/Delta**2.0

  dlgfac_inside_isco=1/(+vt_ - (-1.0)**t_r * sqrt(kr) * vr_ - vp_*lam)
  
  return
end function dlgfac_inside_isco
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function vt(r,a)
  !compute the time-like of the 4-velocity of the disk element
  use isco
  implicit none
  double precision, intent(IN) :: r, a
  double precision             :: vt
  vt = (e_isco*(r**3.0 + r * a**2.0 + 2.0 * a**2.0) -&
       2.0 * j_isco * a) / (r * (r**2.0 - 2. * r + a**2.0))
end function vt
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function vr(r)
  !compute the radial component of the 4-velocity of the disk element
  use isco
  implicit none
  double precision, intent(IN) :: r
  double precision             :: vr
  vr = -(2./(3.*risco))**0.5 * (risco/r - 1.0)**1.5
end function vr
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function vphi(r,a)
  !compute the phi component of the 4-velocity of the disk element
  use isco
  implicit none
  double precision, intent(IN) :: r, a
  double precision             :: vphi
  vphi = (2 * e_isco * a + j_isco * (r - 2) ) / (r * (r**2.0 - 2.* r + a**2.0))
end function vphi
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dlgfacthick(a,mu0,alpha,r,mu)
! Calculates g-factor for a photon travelling from disc to observer.
  ! Disc has constant mu.
  !NOTE: this unction is used only in the Newtonian loop 
  implicit none
  double precision dlgfacthick,a,mu0,alpha,r,mu
  double precision sin0,sindisk,angvel,mathcalA,Delta,Sig
  double precision gtt,gtphi,gphiphi,pt,pphi,num,den
! Useful factors
  sin0     = sqrt( 1.0 - mu0**2 )      
  sindisk  = sqrt( 1.d0 - mu**2 )
  angvel   = 1.0 / ( r**1.5 + abs(a) )
  mathcalA = (r**2+a**2)**2 - (r**2-2.0*r+a**2)*a**2*sindisk
  Delta    = r**2 - 2*r + a**2    
  Sig      = r**2 + a**2 * mu**2
! Metric components
  gtt      = -( 1.0 - 2.0*r / Sig )
  gtphi    = -2.0*a*r*sindisk**2 / Sig
  gphiphi  = mathcalA * sindisk**2 / Sig
! 4-momentum
  pt   = -a * ( alpha*sin0 + a*sindisk**2 )
  pt   = pt + ( r**2 + a**2 ) * ( r**2 + a**2 + a*alpha*sin0 ) / Delta
  pt   = pt / Sig
  pphi = -alpha*sin0/sindisk**2 -a
  pphi = pphi + a * ( r**2 + a**2 + a*alpha*sin0 ) / Delta
  pphi = pphi / Sig
! gdo calc
  num = sqrt( -gtt - 2.0*gtphi*angvel - gphiphi*angvel**2 )
  den = -gtt*pt - gtphi*(pt*angvel+pphi) - gphiphi*pphi*angvel
  dlgfacthick = num / den
  return
end function dlgfacthick
!-----------------------------------------------------------------------


! !-----------------------------------------------------------------------
! function dglpfacthick(r,a,h,mu)
! ! Calculates blue shift expreienced by a photon travelling from
! ! an on-axis point source to a point on a Keplerian disk with constant
! ! scaleheight (h/r=0 is mu=0)
! ! Works for pro- and retrograde spins.
!   implicit none
!   double precision dglpfacthick,r,a,h,mu,gsd
!   double precision angvel,Dh,gphiphi,sindisk,mathcalA,Sig
!   angvel   = 1.0 / ( r**1.5 + abs(a) )
!   Dh       = h**2 - 2*h + a**2
!   sindisk  = sqrt( 1.d0 - mu**2 )
!   mathcalA = (r**2+a**2)**2 - (r**2-2.0*r+a**2)*a**2*sindisk
!   Sig      = r**2 + a**2 * mu**2
!   gphiphi  = mathcalA * sindisk**2 / Sig
!   gsd = 1.0-2.0*r/Sig + 4.0*a*r*sindisk**2/Sig*angvel - gphiphi*angvel**2
!   gsd = Dh/(h**2+a**2) / gsd
!   gsd = sqrt( gsd )
!   dglpfacthick = gsd
! return
! end function dglpfacthick
! !-----------------------------------------------------------------------


!-----------------------------------------------------------------------
function dglpfacthick(r,a,h,mudisk,cosdelta)
! Calculates blue shift expreienced by a photon travelling from
! an on-axis point source to a point on a Keplerian disk with constant
! scaleheight (h/r=0 is mu=0)
! Works for pro- and retrograde spins.
  use blcoordinate
  use isco
  implicit none
  double precision, intent(in) :: r, a, h, mudisk, cosdelta
  double precision dglpfacthick,gsd
  double precision angvel,Dh,sindisk,mathcalA,Sig,Delta,gphiphi
  double precision kr,vt_,vr_,vt,vr,lambda,q,f1234(4),velocity(3),pcros
  double precision r1,mu1,phi1,t1,sigma1
  double precision r2,mu2,phi2,t2,sigma2
  double precision pr, pp, pt
  integer          t_r1, t_r2

  Dh = h**2 - 2*h + a**2
  if(r .ge. risco) then
     sindisk  = sqrt( 1.d0 - mudisk**2 )
     Sig      = r**2 + a**2 * mudisk**2
     mathcalA = (r**2+a**2)**2 - (r**2-2.0*r+a**2)*a**2*sindisk**2
     gphiphi  = mathcalA * sindisk**2 / Sig
     angvel   = 1.0 / ( r**1.5 + a )
     gsd = 1.0-2.0*r/Sig + 4.0*a*r*sindisk**2/Sig*angvel - gphiphi*angvel**2
     gsd = Dh/(h**2+a**2) / gsd
     gsd = sqrt( gsd )
     dglpfacthick = gsd
     !            write(156,*) r,gsd,"out"
  else
     Delta  = r**2. - 2.*r + a**2.0
     velocity = 0.0D0
     pr = cosdelta
     pt = 0.d0
     pp = sqrt(1 - cosdelta**2)
     call initialdirection(pr,pt,pp,0.d0,1.d0,a,h,velocity,lambda,q,f1234)
     pcros = Pemdisk(f1234,lambda,q,0.d0,1.d0,a,h,1.d0,mudisk,1.d15,1.d0+sqrt(1.d0-a**2))
     call YNOGK(pcros-1.d-5,f1234,lambda,q,0.d0,1.d0,a,h,1.d0,&
          r1,mu1,phi1,t1,sigma1, t_r1, t_r2) 
     call YNOGK(pcros,f1234,lambda,q,0.d0,1.d0,a,h,1.d0,&
          r2,mu2,phi2,t2,sigma2, t_r1, t_r2)
     if(pcros.lt.0.d0)then
        dglpfacthick=0.d0
     else
        vt_ = vt(r,a)
        vr_ = vr(r)
        ! kr = sqrt(r**4. - (q + lam**2. - a**2.) * r**2. + &
        !      2. * r * (q + (lam - a)**2.) - a**2.0 * q) / Delta
        !which is equivalent to the power of 2 of the following
        kr = sqrt( (r**2.0 + a**2.0)**2.0 - Delta * (q + a**2.0) )/ Delta
        if (kr .ne. kr) kr = 0.d0 ! this is needed since we used interpolated values of cosdelta for this patch of the disk (if you increase the ndelta, the negative value over the sqrt disappears

        ! write(22,*) r, cosdelta, kr, (r**2.0 + a**2.0)**2.0, Delta * (q + a**2.0), q, lambda, pcros, t_r1
        ! write(40,*) t_r1, (r2-r1) !t_r1 is 1 when (r2-r1) > 0, is 0 when (r2-r1) < 0
        !this is very strange I would have imagined the opposite
        if ((r2-r1).lt.0.d0)then
           dglpfac thick=-sqrt(Dh / (h**2 + a**2) )*(-vt_ - kr * vr_)
        else
           dglpfacthick=-sqrt(Dh / (h**2 + a**2) )*(-vt_ + kr * vr_)
        endif
        ! dglpfacthick=-sqrt(Dh / (h**2 + a**2) )*(-vt_ -(-1.0)**t_r1, kr * vr_)
        
     endif
  endif
  ! write(1003,*)r,dglpfacthick
  return
!  dlgfac_inside_isco=1/(+_vt - (-1.0)**t_r * sqrt(kr)*_vr - _vp*lam)

end function dglpfacthick
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

    
! !-----------------------------------------------------------------------
!       function lorfac(r,a)
! ! Calculates Lorentz factor for rotating disk element
!       implicit none
!       real lorfac,r,a
!       real Delta,BigA,Omega,v
!       Delta  = r**2 - 2*r + a**2
!       BigA   = (r**2+a**2)**2 - a**2*Delta
!       Omega  = 1.0 / ( r**1.5 + abs(a) )
!       v      = ( Omega * bigA - 2*a*r ) / ( r**2 * sqrt(Delta) )
!       lorfac = ( 1 - v**2 )**(-0.5)
!       return
!       end
! !-----------------------------------------------------------------------

      
! !-----------------------------------------------------------------------
! function Delta(r,a)
!   implicit none
!   double precision Delta, r, a
!   Delta    = r**2 - 2*r + a**2
! end function Delta
! !-----------------------------------------------------------------------


! !-----------------------------------------------------------------------
!       function dlorfac(r,a)
! ! Calculates Lorentz factor for rotating disk element
!       implicit none
!       double precision dlorfac,r,a
!       double precision Delta,BigA,Omega,v
!       Delta   = r**2 - 2*r + a**2
!       BigA    = (r**2+a**2)**2 - a**2*Delta
!       Omega   = 1.0 / ( r**1.5 + abs(a) )
!       v       = ( Omega * bigA - 2*a*r ) / ( r**2 * sqrt(Delta) )
!       dlorfac = ( 1 - v**2 )**(-0.5)
!       return
!       end
! !-----------------------------------------------------------------------
      
!-----------------------------------------------------------------------
      function dlorfac(r,a)
! Calculates Lorentz factor for rotating disk element
      use isco
      implicit none
      double precision dlorfac,r,a
      double precision Delta, BigA, vt_, vt
      Delta   = r**2 - 2*r + a**2
      BigA    = (r**2+a**2)**2 - a**2*Delta
      if (r .ge. risco) then
         dlorfac = (a+r**1.5d0 )/sqrt(2.0d0 * a * r**1.5d0 + &
              (-3.0d0 + r) * r**2.0d0) &
              * sqrt((r * Delta) / (r**3.d0 + (2.d0 + r)*a**2.d0))
      else
         vt_ = vt(r,a)
         dlorfac = vt_ * sqrt( (r * Delta) / (r**3.d0 + (2.d0 + r) * a**2.d0) )
      endif
      return
      end
!-----------------------------------------------------------------------
