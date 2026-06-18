module raytracing
!> This module contains the shim subroutines for raytracing and lensing factor 
!> calculations. All the subroutines only serve as interfaces to the actual 
!> implementations but allows for easy substitution without changing the rest of
!> the code.
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
        double precision :: tm1, tm2
        double precision :: tr1, tr2
        call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                   radi,mu,phi,time,sigma)
        return
    end subroutine raytrace_disk

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
        implicit none
        integer, intent(in) :: nro,nphi
        double precision, intent(in) :: rn(nro),mueff,mu0,spin,rmin,rout
        double precision, intent(in) :: mudisk,d
        double precision, intent(in) :: phin,alpha,beta,cos0,sin0,scal
        double precision, intent(in) :: velocity(3),f1234(4)
        double precision, intent(out) :: lambda,q
        double precision pem,re,mucros,phie,taudo,sigmacros      
        integer i,j
        cos0  = mu0
        sin0  = sqrt(1.0-cos0**2)
        scal     = 1.d0
        velocity = 0.d0
        taudo1   = 0.0
        re1      = 0.0      
        do i = 1,nro
            do j = 1,NPHI
                phin  = (j-0.5) * 2.d0 * pi / dble(nphi)
                alpha = rn(i) * sin(phin)
                beta  = -rn(i) * cos(phin) * mueff
                call constants_of_motion(-alpha,-beta,d,sin0,cos0,spin,scal,   &
                                        velocity,f1234,lambda,q)
                !Can try rin instead of rmin to save an if statement
                pem = Pemdisk(f1234,lambda,q,sin0,cos0,spin,d,scal,mudisk,rout,&
                            rmin)  
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

end module raytracing