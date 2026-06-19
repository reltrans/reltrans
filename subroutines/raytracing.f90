module raytracing
!> This module contains the shim subroutines for raytracing and lensing factor 
!> calculations. All the subroutines only serve as interfaces to the actual 
!> implementations but allows for easy substitution without changing the rest of
!> the code.
    use rtconstants, only: pi
    implicit none

contains

    subroutine raytrace_disk(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                        radi,mu,phi,time,sigma)
    !> Computes four Boyer-Lindquist coordinates (r,\mu,\phi,t) and affine 
    !> parameter \sigma as functions of parameter p, i.e. functions r(p), 
    !> \mu(p), \phi(p), t(p), \sigma as functions of parameter p, i.e. 
    !> functions r(p), \mu(p), \phi(p), t(p) and \sigma(p). Cf. discussions in 
    !> Yang & Wang (2012).
    !> Inputs:
    !>    p: independent variable, which must be nonnegative.
    !>    f1234: array of p_1, p_2, p_3, p_4, which are the components of 
    !>           four-momentum of a photon measured under the LNRF frame. This 
    !>           array can be computed by subroutine lambdaq(...), see below
    !>    lambda,q: motion constants, defined by lambda=L_z/E, q=Q/E^2.
    !>    sinobs,muobs: sinobs=sin(\theta_{obs}), muobs=cos(\theta_{obs}), where 
    !>                  \theta_{obs} is the inclination angle of the observer.
    !>    a_spin: spin of black hole, on interval (-1,1).
    !>    robs: radial coordinate of observer or initial position of photon.
    !>    scal: a dimensionless parameter to control the size of the images, 
    !>          which is usually set to 1.D0.
    !> Outputs:
    !>    radi: value of function r(p).
    !>    mu: value of function \mu(p).
    !>    phi: value of function \phi(p).
    !>    time: value of function t(p).
    !>    sigma: value of function \sigma(p).
    !>    tm1,tm2: number of times of photon meets turning points \mu_tp1 and 
    !>             \mu_tp2 respectively.
    !>    tr1,tr2: number of times of photon meets turning points r_tp1 and 
    !>             r_tp2 respectively.
        use blcoordinate, only: YNOGK
        implicit none
        double precision, intent(in) :: p
        double precision, intent(in) :: f1234(4), lambda, q, sinobs, muobs
        double precision, intent(in) :: a_spin, robs, scal
        double precision, intent(out) :: radi, mu, phi, time, sigma
        call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                   radi,mu,phi,time,sigma)
        return
    end subroutine raytrace_disk

    
    function p_disk_crossing(f1234,lambda,q,sins,mus,a_spin,h,scal,mudisk,     &
                 r_max,r_min)
        use blcoordinate, only: Pemdisk
        implicit none
        double precision, intent(in) :: f1234(4), lambda, q, sins, mus
        double precision, intent(in) :: a_spin, h, scal, mudisk, r_max, r_min
        double precision p_disk_crossing
        p_disk_crossing = Pemdisk(f1234,lambda,q,sins,mus,a_spin,h,            &
               scal,mudisk,r_max,r_min)
        return
    end function p_disk_crossing

    
    subroutine initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,    &
                                f1234)
    !> This subroutine is a shim function that allows for the substitution of 
    !> the initial_photon functionality without adjusting the rest of the code.
    !> Inputs:
    !>     pr, pt, pp: components of the initial photon momentum in the source
    !>                 rest frame.
    !>     sins, mus: sine and cosine of the source inclination angle.
    !>     a_spin: spin of the black hole.
    !>     h: height of the source above the black hole.
    !>     velocity: 3-velocity of the source.
    !> Outputs:
    !>     lambda, q: motion constants, defined by lambda=L_z/E, q=Q/E^2.
    !>     f1234: array of p_1, p_2, p_3, p_4, which are the components of 
    !>            four-momentum of a photon measured under the LNRF frame. This 
    !>            array can be computed by subroutine lambdaq(...), see below
        use blcoordinate, only: initialdirection
        implicit none
        double precision, intent(in) :: pr, pt, pp, sins, mus, a_spin, h
        double precision, intent(in) :: velocity(3)
        double precision, intent(out) :: lambda, q
        double precision, intent(out) :: f1234(4)
        call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,    &
                                f1234)
        return
    end subroutine initial_photon

    function p_coord_at_infinity(f1234,lambda,q,sins,mus,a_spin,h,scal)        &
                                result(p_coord_at_inf)
    !> This subroutine is a shim function that allows for the substitution of 
    !> the p_coord_at_infinity functionality without adjusting the rest of the code.
    !> Inputs:
    !>     f1234: array of p_1, p_2, p_3, p_4, which are the components of 
    !>            four-momentum of a photon measured under the LNRF frame.
    !>     lambda, q: motion constants, defined by lambda=L_z/E, q=Q/E^2.
    !>     sins, mus: sine and cosine of the source inclination angle.
    !>     a_spin: spin of the black hole.
    !>     h: height of the source above the black hole.
    !>     scal: a dimensionless parameter to control the size of the images,
    !>           which is usually set to 1.D0.
    !> Outputs:
    !>     p_coord_at_infinity: array of p_1, p_2, p_3, p_4, which are the 
    !>                         components of four-momentum of a photon measured 
    !>                         at infinity.
        use blcoordinate, only: p_total
        implicit none
        double precision, intent(in) :: f1234(4), lambda, q, sins, mus
        double precision, intent(in) :: a_spin, h, scal
        double precision :: p_coord_at_inf
        p_coord_at_inf = p_total(f1234(1),lambda,q,sins,mus,a_spin,h,scal)
    end function p_coord_at_infinity

    subroutine constants_of_motion(alpha,beta,robs,sinobs,muobs,a_spin,scal,   &
                                    velocity,f1234,lambda,q)
    !> This subroutine is a shim function that allows for the substitution of 
    !> the constants_of_motion functionality without adjusting the rest of the 
    !> code.
    !> Inputs:
    !>     alpha, beta: impact parameters of the photon.
    !>     robs: radial coordinate of observer or initial position of photon.
    !>     sinobs, muobs: sine and cosine of the observer inclination angle.
    !>     a_spin: spin of the black hole.
    !>     scal: a dimensionless parameter to control the size of the images,
    !>           which is usually set to 1.D0.
    !>     velocity: 3-velocity of the source.
    !> Outputs:
    !>     f1234: array of p_1, p_2, p_3, p_4, which are the components of 
    !>            four-momentum of a photon measured under the LNRF frame. This 
    !>            array can be computed by subroutine lambdaq(...), see below
    !>     lambda, q: motion constants, defined by lambda=L_z/E, q=Q/E^2.
        use blcoordinate, only: lambdaq
        implicit none
        double precision, intent(in) :: alpha, beta, robs, sinobs, muobs
        double precision, intent(in) :: a_spin, scal
        double precision, intent(in) :: velocity(3)
        double precision, intent(out) :: f1234(4), lambda, q
        call lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,&
                    lambda,q)
        return
    end subroutine constants_of_motion

    subroutine trace_disk_observer(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,&
                                    d)
    !> Traces rays in full GR for the camera defined by rn(nro), nro, nphi
    !> to convert alpha and beta to r and tau_do (don't care about phi)
    !> Used to be called GRtrace.
    !> Inputs:
    !>     nro: number of radial points in the disk.
    !>     nphi: number of azimuthal points in the disk.
    !>     rn: array of radial points in the disk.
    !>     mueff: effective cosine of the inclination angle of the disk.
    !>     mu0: cosine of the inclination angle of the observer.
    !>     spin: spin of the black hole.
    !>     rmin: minimum radius of the disk.
    !>     rout: maximum radius of the disk.
    !>     mudisk: cosine of the inclination angle of the disk.
    !>     d: distance from the black hole to the observer.
    !> Outputs:
    !>     pem1: array of p-coordinate at the disk for each ray.
    !>     taudo1: array of time coordinate at the disk for each ray.
    !>     re1: array of radial coordinate at the disk for each ray.
        !use blcoordinate, only: Pemdisk
        use dyn_gr
        implicit none
        integer, intent(in) :: nro,nphi
        double precision, intent(in) :: rn(nro),mueff,mu0,spin,rmin,rout
        double precision, intent(in) :: mudisk,d
        double precision phin,alpha,beta,cos0,sin0,scal
        double precision velocity(3),f1234(4)
        double precision lambda,q
        double precision pem,re,mucros,phie,taudo,sigmacros   
        integer i,j
        cos0  = mu0
        sin0  = sqrt(1.0-cos0**2)
        scal     = 1.d0
        velocity = 0.d0
        taudo1   = 0.0
        re1      = 0.0
        !TODO: kerrz optimisation here! 
        do i = 1,nro
            do j = 1,NPHI
                phin  = (j-0.5) * 2.d0 * pi / dble(nphi)
                alpha = rn(i) * sin(phin)
                beta  = -rn(i) * cos(phin) * mueff
                call constants_of_motion(-alpha,-beta,d,sin0,cos0,spin,scal,   &
                                        velocity,f1234,lambda,q)
                !Can try rin instead of rmin to save an if statement
                pem = p_disk_crossing(f1234,lambda,q,sin0,cos0,spin,d,scal,    &
                            mudisk,rout,rmin)  
                pem1(j,i) = pem
                !pem > 0 means there is a solution
                !pem < 0 means there is no solution
                if( pem .gt. 0.0d0 )then
                    call raytrace_disk(pem,f1234,lambda,q,sin0,cos0,spin,d,    &
                                        scal,re,mucros,phie,taudo,sigmacros)
                    taudo1(j,i) = taudo - d
                    re1(j,i)    = re
                 end if
              end do
        end do
        return
    end subroutine trace_disk_observer

    
    subroutine getdcos(a_spin,h,mudisk,n,nlp,rout,npts,r1,dcosdr,tc,cosd1,cosdout)
    !> For n values of the emission angle, delta, the code calculates the r and t 
    !> coordinates for the geodesic for mu=mudisk; i.e. the crossing points of a 
    !> thin disk.
    !> Note that mudisk = (h/r) / sqrt( (h/r)**2 + 1 )
    !> INPUTS
    !>    a_spin       Dimensionless spin parameter
    !>    h            Height of on-axis, isotropically emitting source
    !>    mudisk       cos(theta) of disk surface (mu=0 for h/r=0)
    !>    n            Number of values of emission angle delta (see Fig 1 Dauser 
    !>                 et al 2013) calculated
    !>    rout         Disk outer radius
    !> OUTPUTS
    !>    npts         Number of points recorded in arrays (leq n, since some 
    !>                 trial values will not hit the disk)
    !>    r1(n)        Radius of disk crossing
    !>    dcosdr(n)    Corresponding d\cos\delta/dr
    !>    tc(n)        Corresponding time coordinate
    !>    cosd1(n)     Corresponding \cos\delta
    !>    cosdout      cosd at the outer disk radius
        implicit none
        double precision, intent(in )   :: a_spin, h(2), mudisk, rout
        integer         , intent(in )   :: n, nlp
        integer         , intent(inout) :: npts(nlp)
        double precision, intent(inout) :: r1(n,nlp)
        double precision, intent(out)   :: dcosdr(n,nlp), tc(n,nlp), cosd1(n,nlp), cosdout(nlp)
        integer  m,j,k,counter,nout(nlp)
        double precision sins,mus,lambda,q,scal
        double precision rhorizon,velocity(3),f1234(4),pp,pr,pt
        double precision deltamin,deltamax, deltas,r_min,r_max,disco
        double precision rcros,mucros,phicros,tcros,sigmacros,pcros
        scal     = 1.d0   !Meaningless scaling factor
        mus      = 1.d0   !Position of source: mus=1 means on-axis
        sins     = 0.d0   !sin of same angle
        velocity = 0.0D0  !3-velocity of source
        rhorizon = 1.d0+sqrt(1.d0-a_spin**2)

        !loop over h here
        do m=1,nlp
            !Calculate smallest delta worth considering
            deltamin = acos( h(m) / sqrt( h(m)**2 + rhorizon**2 ) )
            !Consider arbitrarily large value of delta
            deltamax = pi
            !Set minimum and maximum disk radii
            r_min = disco( a_spin )
            r_max = 1d10
            !Go through n different values of the angle delta_s
            counter = 0
            nout(m) = 1
            do j = 1,n
            !Run through linear steps in the angle delta (see Fig 1; Dauser et al 2013)
                deltas   = deltamin + (j-1) * (deltamax-deltamin)/float(n-1)
                !Calculate 4-momentum in source rest frame tetrad
                pr = cos(deltas)           !cosdelta
                pp = sqrt( 1.d0 - pr**2 )  !sindelta
                pt= 0.d0
                !Convert to LNRF (locally non-rotating reference frame)
                call initial_photon(pr,pt,pp,sins,mus,a_spin,h(m),             &
                     velocity,lambda,q,f1234)
                !Calculate value of p-coordinate at mu=0
                pcros = p_disk_crossing(f1234,lambda,q,sins,mus,               &
                                    a_spin,h(m),scal,mudisk,r_max,r_min)
                !From that, calculate r, phi and t at mu=0
                call raytrace_disk(pcros,f1234,lambda,q,sins,mus,a_spin,       &
                              h(m),scal,rcros,mucros,phicros,tcros,sigmacros)
                if( pcros .gt. 0.0 )then
                    counter        = counter + 1
                    r1(counter,m)    = rcros
                    cosd1(counter,m) = pr    !cosdelta
                    tc(counter,m)    = tcros
                    if( rout .gt. r1(counter,m) ) nout(m) = counter
                end if
            end do 
            npts(m) = counter
        end do 
        
        !Calculate cosdout
        do m=1,nlp 
           if( nout(m) .eq. npts(m) )then
            !Extrapolate assuming Newtonian profile
                cosdout(m) = h(m)/sqrt(h(m)**2+rout**2)-h(m)/                  &
                            sqrt(h(m)**2+r1(npts(m),m)**2)+cosd1(npts(m),m)
            else
                !Inperpolate
                cosdout(m) = cosd1(nout(m)+1,m) - cosd1(nout(m),m)
                cosdout(m) = cosdout(m) * ( rout - r1( nout(m),m ) )
                cosdout(m) = cosdout(m) / ( r1(nout(m)+1,m) - r1(nout(m),m) )
                cosdout(m) = cosdout(m) + cosd1(nout(m),m)
            end if
        end do
        !Calculate d\delta/dr on the r-grid. Note that we need yet another loop 
        !over m because of how the counter npts is set up computationally, this 
        !costs no time whatsoever
        do m=1,nlp            
            npts(m) = npts(m) -1           
            do k = 1,npts(m)
                dcosdr(k,m) = abs((cosd1(k+1,m)-cosd1(k,m))/(r1(k+1,m)-r1(k,m)))
            end do
        end do
        !Discard the outer points as unreliable
        npts = npts - 7

        return
    end subroutine getdcos


    function cosidel(cosdelta,sins,mus,a_spin,h,velocity)
    !> Calculates cosi when given cosdelta and parameters
    !> Inputs:
    !> cosdelta,sins,mus,a_spin,h,velocity
        implicit none
        double precision cosdelta,sins,mus,a_spin,h,velocity(3),cosidel
        double precision pr,pp,pt,lambda,q,f1234(4),ptotal
        double precision scal,p,ra,mua,phya,timea,sigmaa
        double precision p_coord_at_inf
        scal = 1.d0                  !Meaningless scaling factor
        pr   = cosdelta              !cosdelta
        pp   = sqrt( 1.d0 - pr**2 )  !sindelta
        pt   = 0.d0
        !Convert to LNRF (locally non-rotating reference frame)
        call initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
        !Now calculate ptotal (value of p-coordinate at infinity)
        p_coord_at_inf = p_coord_at_infinity(f1234,lambda,q,sins,mus,a_spin,h, &
                                            scal)
        p = 0.9999d0 * p_coord_at_inf
        call raytrace_disk(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
                 ra,mua,phya,timea,sigmaa)
        cosidel = mua
        return
    end function cosidel


    subroutine getlimits(sins,mus,a_spin,h,velocity,muobs,x1,x2)
        !> Minimisation routine will numerically calculate cosdelta for a given cosi.
        !> To do that, we need limits that bracket only one root. 
        !> This routine works out sensible limits
        !> Inputs:
        !>     sins, mus: sine and cosine of the source inclination angle.
        !>     a_spin: spin of the black hole.
        !>     h: height of the source above the black hole.
        !>     velocity: 3-velocity of the source.
        !>     muobs: cosine of the observer inclination angle.
        !> Outputs:
        !>     x1, x2: limits for the minimisation routine.
        implicit none
        double precision sins,mus,a_spin,h,velocity(3),muobs,x1,x2
        double precision cosdelta0,mua,cosi,cosdelta
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
    

    subroutine getlens(a_spin,h,muobs,lens,delt,cosdelta1)
        !> Routine to calculate the lensing factor l=d\cos\delta/d\cos(i)
        !> and the source to observer time lag.
        !> Both calculations need us to know the delta value for the geodesic
        !> that ends up at angle i at infinity.
        !> INPUTS
        !> a_spin       Dimensionless spin parameter
        !> h            Height of on-axis, isotropically emitting source
        !> muobs        Cosine of inclination angle
        !
        !> OUTPUTS
        !> lens         Lensing factor
        !> delt         Source to observer time lag 
        implicit none
        double precision, intent(in)    :: a_spin,h, muobs
        double precision, intent(inout) :: cosdelta1
        double precision, intent(out)   :: lens, delt
        double precision sins,mus,lambda,q,scal
        double precision velocity(3),f1234(4),pp,pr,pt
        double precision ptotal,dcosdelta,drtbis
        double precision mua,p,phya,ra,sigmaa,timea
        double precision par(3),x1,x2,xacc,mu2
        double precision alpha,beta,b2,d
        double precision p_coord_at_inf

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

        !Set 4-momentum in the source frame
        pr   = cosdelta1             !cosdelta
        pp   = sqrt( 1.d0 - pr**2 )  !sindelta
        pt   = 0.d0
        !Convert to LNRF (locally non-rotating reference frame)
        call initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
        !Now calculate ptotal (value of p-coordinate at infinity)
        p_coord_at_inf = p_coord_at_infinity(f1234,lambda,q,sins,mus,a_spin,h, &
                                            scal)
        p = 0.9999d0 * p_coord_at_inf
        call raytrace_disk(p,f1234,lambda,q,sins,mus,a_spin,h,scal,&
             ra,mua,phya,timea,sigmaa)
        !Calcluate the distance from BH to centre of observer's camera
        !For an on-axis lamppost, alpha should always be 0, but the below is general
        if( muobs .eq. 1.d0 )then
            d = ra
        else
            alpha = -lambda / sqrt( 1.0 - muobs**2 )
            beta  = q - (alpha**2-a_spin**2)*muobs**2
            beta = sqrt(beta)
            b2   = alpha**2 + beta**2
            d    = sqrt( ra**2 - b2  )
        end if
        !Subtract the distance - will do the same for
        !the disk to observer lags, meaning I don't need to use the
        !same distance for both calculations
        delt = timea - d
        return
    end subroutine getlens


    function mudiff(cosdelta,par)
        !> Calculates muobs (cosine of distant inclination angle) when given
        !> cos(delta) (cosine of angle between initial photon trajectory and -z)
        !> Inputs:
        !>     cosdelta: cosine of angle between initial photon trajectory and -z
        !>     par: array of parameters, where par(1)=a_spin, par(2)=h, par(3)=muobs
        !> Outputs:
        !>     mudiff: difference between calculated muobs and input muobs
        implicit none
        double precision mudiff,cosdelta,par(3)
        double precision a_spin,h,muobs
        double precision scal,mus,sins
        double precision velocity(3),sindelta,pp,pr,pt,lambda,q,f1234(4)
        double precision ptotal,x,y,z,xprev,yprev,zprev,delx,dely,delz
        double precision ra,mua,phya,timea,sigmaa,p,cosdum
        double precision p_coord_at_inf
        a_spin = par(1)
        h      = par(2)
        muobs  = par(3)
        velocity = 0.0D0
        sindelta = sqrt( 1.d0 - cosdelta**2 )
        if ( sindelta .eq. 0.d0 )then
            cosdum = 1.d0
        else
            !Calculate 4-momentum in source rest frame tetrad
            pp= sindelta
            pr= cosdelta
            pt= 0.d0
            !Convert to LNRF (locally non-rotating reference frame)
            scal = 1.d0
            mus  = 1.d0
            sins = 0.d0
            call initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,           &
                                  lambda,q,f1234)
            p_coord_at_inf = p_coord_at_infinity(f1234,lambda,q,sins,mus,      &
                                  a_spin,h,scal)
            p = 0.9998d0 * p_coord_at_inf
            call raytrace_disk(p,f1234,lambda,q,sins,mus,a_spin,h,scal,        &
                       ra,mua,phya,timea,sigmaa)
            xprev = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*cos(phya)
            yprev = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*sin(phya)
            zprev = ra*mua
            p = 0.9999d0 * p_coord_at_inf
            call raytrace_disk(p,f1234,lambda,q,sins,mus,a_spin,h,scal,        &
                       ra,mua,phya,timea,sigmaa)
            x = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*cos(phya)
            y = sqrt(ra**2+a_spin**2)*sqrt(1.d0-mua**2)*sin(phya)
            z = ra*mua
            delx = x - xprev
            dely = y - yprev
            delz = z - zprev
            cosdum = delz / sqrt( delx**2 + dely**2 + delz**2 )
        end if
    
        mudiff = cosdum - muobs
        return
    end function mudiff

    
end module raytracing
