!-----------------------------------------------------------------------
      subroutine sourcelag(a_spin,h,muobs,delt)
      use blcoordinate
      implicit none
      double precision a_spin,h,muobs
      integer nup,k,j
      double precision scal,mus,sins
      double precision velocity(3),cosdelta,sindelta,pp,pr,pt,lambda,q,f1234(4)
      double precision ptotal,x,y,z,xprev,yprev,zprev,delx,dely,delz
      double precision ra,mua,phya,timea,sigmaa,p,cosdum,d,diff,mudiff,par(3)
      double precision x1,x2,xacc,drtbis,dfar,tfar,rfar,sinOa,delt
      external mudiff
      scal = 1.d0
      !Initialize for a stationary, on-axis source
      mus  = 1.d0
      sins = 0.d0
      velocity = 0.d0
      !Find cosdelta that sends a geodesic to the observer at muobs
      par(1)=a_spin
      par(2)=h
      par(3)=muobs
      x1   = -1.d0  !0.998
      x2   = -muobs
      xacc = 1d-4
      cosdelta = drtbis(mudiff,x1,x2,xacc,par)
      !Now calculate time and distance along the alpha=beta=0 geodesic
      !for a point on the geodesic where it is already straight (i.e. a large p)
      sindelta = sqrt( 1.d0 - cosdelta**2 )
      !Calculate 4-momentum in source rest frame tetrad
      pp= sindelta
      pr= cosdelta
      pt= 0.d0
      call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
      ptotal = p_total(f1234(1),lambda,q,sins,mus,a_spin,h,scal)
      p = 0.9998 * ptotal
      call YNOGK(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
           ra,mua,phya,timea,sigmaa)
      !Now calculate the time taken by a straight ray
      sinOa = sqrt(1.d0-mua**2)
      x = sqrt(ra**2+a_spin**2)*sinOa*cos(phya)
      y = sqrt(ra**2+a_spin**2)*sinOa*sin(phya)
      z = ra*mua
      dfar = x*sinOa*cos(phya) + y*sinOa*sin(phya) + z*muobs
      !Output GR time - Newtonian time + d*c
      tfar = timea
      delt = tfar - dfar
      return
      end subroutine sourcelag
!-----------------------------------------------------------------------
