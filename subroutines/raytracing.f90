module raytracing
    implicit none

contains

    subroutine raytrace(model_args)
    !> This subroutine is a shim function that allows for the substitution of the
    !> raytracing functionality without adjusting the rest of the code
    !> Inputs:
    !>   model_args: a real array containing the model parameters
    !> Outputs:
    !>   none
        use grtrace, only: GRtrace
        implicit none
        real, intent(in) :: model_args(:)
        ! Call the GRtrace subroutine
        call GRtrace(nro,nphi,rn,mueff,mu0,spin,rmin,rout,mudisk,d)
        return
    end subroutine raytrace

    subroutine lensing_factor(model_args)
    !> This subroutine is a shim function that allows for the substitution of the
    !> lensing factor functionality without adjusting the rest of the code
    !> Inputs:
    !>   model_args: a real array containing the model parameters
    !> Outputs:
    !>   none
        use getlens, only: getlens
        implicit none
        real, intent(in) :: model_args(:)
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
        real, intent(in) :: cosdelta(:)
        real, intent(in) :: par(:)
        ! Call the mudiff subroutine
        call mudiff(cosdelta,par)
        return
    end subroutine distant_inclination

end module raytracing