module raytracing
!> This module contains the shim subroutines for raytracing and lensing factor 
!> calculations. All the subroutines only serve as interfaces to the actual 
!> implementations but allows for easy substitution without changing the rest of
!> the code.
    implicit none

contains

    subroutine raytrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
    !> This subroutine is a shim function that allows for the substitution of the
    !> raytracing functionality without adjusting the rest of the code
    !> Inputs:
    !>   model_args: a real array containing the model parameters
    !> Outputs:
    !>   none
        use grtrace, only: GRtrace
        implicit none
        integer, intent(in) :: nro, nphi
        double precision, intent(in) :: rn, mueff, mu0, spin, rmin, rout
        double precision, intent(in) :: mudisk, d
        ! Call the GRtrace subroutine
        call GRtrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
        return
    end subroutine raytrace

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
        use YNOGK, only: YNOGK
        implicit none
        double precision, intent(in) :: p
        double precision, intent(in) :: f1234(4), lambda, q, sinobs, muobs
        double precision, intent(in) :: a_spin, robs, scal
        double precision, intent(out) :: radi, mu, phi, time, sigma
        double precision :: tm1, tm2
        double precision :: tr1, tr2
        call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                   radi,mu,phi,time,sigma,tm1,tm2,tr1,tr2)
        return
    end subroutine raytrace_disk

    subroutine initial_photon(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,    &
                                f1234)
    !> This subroutine is a shim function that allows for the substitution of the
    !> initialdirection functionality without adjusting the rest of the code.
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
        call initialdirection(pr,pt,pp,sins,mus,a_spin,h,velocity,lambda,q,f1234)
        return
    end subroutine initial_photon

end module raytracing