!*****************************************************************************************************
      subroutine getlens(a_spin,h,muobs,lens,delt)
! Routine to calculate the lensing factor l=d\cos\delta/d\cos(i)
! and the source to observer time lag.
! Both calculations need us to know the delta value for the geodesic
! that ends up at angle i at infinity.
! INPUTS
! a_spin       Dimensionless spin parameter
! h            Height of on-axis, isotropically emitting source
! muobs        Cosine of inclination angle
!
! OUTPUTS
! lens         Lensing factor
! delt         Source to observer time lag 
      use blcoordinate
      implicit none
      double precision sins,mus,a_spin,h,lambda,q,scal
      double precision velocity(3),f1234(4),pp,pr,pt
      double precision muobs,delt,drtbis,cosidel
      double precision ptotal,cosdelta1,dcosdelta,lens
      double precision mua,p,phya,ra,sigmaa,timea,mudiff
      double precision par(3),x1,x2,xacc,mu2
!      double precision cosdelta,cosdelta0,cosi,cosphi,costheta,deltamax,deltamin,mu1,newt,sinphi,sintheta
      external mudiff
      !Settings
      scal      = 1.d0   !Meaningless scaling factor
      mus       = 1.d0   !Position of source: mus=0 means on-axis
      sins      = 0.d0   !sin of same angle
      velocity  = 0.0D0  !3-velocity of source
      dcosdelta = 1d-2   !Step in cosdelta used for differentiation
      xacc      = 1d-6   !Accuracy of minimisation routine

      !First calculate the cosdelta corresponding to the input muobs
      
      !Set limits for minimisation routine
      call getlimits(sins,mus,a_spin,h,velocity,muobs,x1,x2)
      !Call minimisation routine
      par(1)=a_spin
      par(2)=h
      par(3)=muobs
      cosdelta1 = drtbis(mudiff,x1,x2,xacc,par)
      
      !Now calculate the lensing factor
      
      !Make cosdelta a little bit bigger and calculate the new cosi
      mu2 = cosidel(cosdelta1+dcosdelta,sins,mus,a_spin,h,velocity) 
      !Finally calculate the lensing factor
      lens = dcosdelta / ( muobs - mu2 )

      !Now calculate the source lag
      pr   = cosdelta1             !cosdelta
      pp   = sqrt( 1.d0 - pr**2 )  !sindelta
      pt   = 0.d0
      !Convert to LNRF (locally non-rotating reference frame)
      call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
      !Now calculate ptotal (value of p-coordinate at infinity)
      ptotal = p_total(f1234(1),lambda,q,sins,mus,a_spin,h,scal)
      p = 0.9999d0 * ptotal
      call YNOGK(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
           ra,mua,phya,timea,sigmaa)
      !Finally subtract the distance used - will do the same for
      !the disk to observer lags, meaning I don't need to use the
      !same distance for both calculations
      delt = timea - sqrt( ra**2 + h**2*(muobs**2-1.0) )      
      return
      end subroutine getlens
!*****************************************************************************************************


!-----------------------------------------------------------------------
      subroutine getlimits(sins,mus,a_spin,h,velocity,muobs,x1,x2)
! Minimisation routine will numerically calculate cosdelta for a given cosi.
! To do that, we need limits that bracket only one root. 
! This routine works out sensible limits
      implicit none
      double precision sins,mus,a_spin,h,velocity(3),muobs,x1,x2
      double precision cosdelta0,mua,cosidel,cosi,cosdelta
      !The first limit is always cosdelta=-1 (corresponding to cosi=1)
      !Can't take cosdelta too large because this will also braket
      !the ghost images solutions
      !Tactic: extrapolate the initially straight line function from
      !cosi = 1, to some well-chosen cosi value. The cosdelta resulting
      !From this extrapolation is my second limit.
      cosdelta0 = -0.98d0
      mua = cosidel(cosdelta0,sins,mus,a_spin,h,velocity)
      !Take the straight line from (cosi=1,cosdelta=-1) to (cosi=mua,cosdelta=cosdelta0)
      !and extrapolate down to cosi=-0.5
      cosi = -0.5
      cosdelta = (cosi-1.d0)*(cosdelta0+1.d0)/(mua-1.d0) - 1.0
      cosdelta = min( cosdelta , -muobs )  !-muobs is the Newtonian limit 
      !Use for limits
      x1 = -1.d0
      x2 = cosdelta
      return
      end subroutine getlimits
!-----------------------------------------------------------------------

      
!-----------------------------------------------------------------------
      function cosidel(cosdelta,sins,mus,a_spin,h,velocity)
! Inputs:
! cosdelta,sins,mus,a_spin,h,velocity
        
! Calculates cosi when given cosdelta and parameters
!
!        
      use blcoordinate
      implicit none
      double precision cosdelta,sins,mus,a_spin,h,velocity(3),cosidel
      double precision pr,pp,pt,lambda,q,f1234(4),ptotal
      double precision scal,p,ra,mua,phya,timea,sigmaa
      scal = 1.d0                  !Meaningless scaling factor
      pr   = cosdelta              !cosdelta
      pp   = sqrt( 1.d0 - pr**2 )  !sindelta
      pt   = 0.d0
      !Convert to LNRF (locally non-rotating reference frame)
      call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
      !Now calculate ptotal (value of p-coordinate at infinity)
      ptotal = p_total(f1234(1),lambda,q,sins,mus,a_spin,h,scal)
      p = 0.9999d0 * ptotal
      call YNOGK(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
           ra,mua,phya,timea,sigmaa)
      cosidel = mua
      return
      end function cosidel
!-----------------------------------------------------------------------
      
