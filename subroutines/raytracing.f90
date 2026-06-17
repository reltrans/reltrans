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

    subroutine lensing_factor(a_spin,h,muobs,lens,delt,cosdelta1)
    !> This subroutine is a shim function that allows for the substitution of the
    !> lensing factor functionality without adjusting the rest of the code
    !> Inputs:
    !>   model_args: a real array containing the model parameters
    !> Outputs:
    !>   none
        use getlens, only: getlens
        implicit none
        double precision, intent(in) :: a_spin, h, muobs, lens, delt, cosdelta1
        ! Call the getlens subroutine
        call getlens(a_spin,h,muobs,lens,delt,cosdelta1)
        return
    end subroutine lensing_factor

    subroutine distant_inclination(cosdelta, par)
    !> This subroutine is a shim function that allows for the substitution of the
    !> distant inclination functionality without adjusting the rest of the code
    !> Inputs:
    !>   cosdelta: a real array containing the cosine of the inclination angles
    !>   par: a real array containing additional parameters
    !> Outputs:
    !>   none
        use mudiff, only: mudiff
        implicit none
        double precision, intent(in) :: cosdelta(:)
        double precision, intent(in) :: par(:)
        ! Call the mudiff subroutine
        call mudiff(cosdelta,par)
        return
    end subroutine distant_inclination

    subroutine raytrace_disk
    !> This subroutine is a shim function that allows for the substitution of the
    !> raytracing disk functionality without adjusting the rest of the code
    !> Inputs:
    !>   none
    !> Outputs:
    !>   radi-----------value of function r(p). 
    !>   mu-------------value of function \mu(p). 
    !>   phi------------value of function \phi(p).
    !>   time-----------value of function t(p).
    !>   sigma----------value of function \sigma(p).
    !>   tm1,tm2--------number of times of photon meets turning points \mu_tp1 and \mu_tp2
    !>                  respectively. 
    !>   tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
    !>                  respectively.
        use YNOGK, only: YNOGK
        call YNOGK

end module raytracing
!*                              respectively. 
!*               tr1,tr2--------number of times of photon meets turning points r_tp1 and r_tp2
!*                              respectively.            

end module raytracing